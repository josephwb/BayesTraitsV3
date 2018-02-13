/*
*  BayesTriats 3.0
*
*  copyright 2017
*
*  Andrew Meade
*  School of Biological Sciences
*  University of Reading
*  Reading
*  Berkshire
*  RG6 6BX
*
* BayesTriats is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef NO_GSL
	#include <gsl/gsl_randist.h>
	#include <gsl/gsl_cdf.h>
#endif

#include "TypeDef.h"
#include "GenLib.h"
#include "Stones.h"
#include "Options.h"

void	PrintStones(FILE *Str, STONES *Stones)
{
	fprintf(Str, "Steppingstone sampler:\n");
	fprintf(Str, "        No Stones:                  %d\n", Stones->NoStones);
	fprintf(Str, "        Start It:                   %llu\n", Stones->ItStart);
	fprintf(Str, "        It Per Stone:               %d\n", Stones->ItPerStone);
	fprintf(Str, "        Sample Freq:                %d\n", Stones->SampleFreq);
	fprintf(Str, "        Dist:                       Beta(%f,%f)\n", Stones->Alpha, Stones->Beta);
}

void	OutputStoneHeadder(FILE *Out, STONES *Stones)
{
	
	PrintStones(Out, Stones);

	fprintf(Out, "Stone No\tPower\tN\tStone MLh\tRunning MLh\n");

	fflush(Out);
}

double	GetStoneHeat(STONES *Stones, long long Itter, double Heat)
{
	long long StoneIt;
	int CStone;

	if(Itter < Stones->ItStart)
		return Heat;
	
	StoneIt = Itter - Stones->ItStart;
	CStone = (int) StoneIt / Stones->ItPerStone;

	return Heat * Stones->Power[CStone];
}

int		StonesStarted(STONES *Stones, long long Itter)
{
	

	if(Stones == NULL)
		return FALSE;

	if(Itter > Stones->ItStart)
		return TRUE;

	return FALSE;
}

int		ChangeSample(STONES *Stones, int Itters)
{
	if(Stones == NULL)
		return TRUE;

	if(Itters < Stones->ItStart)
		return TRUE;

	return FALSE;
}

void	FreeStones(STONES *Stones)
{
	free(Stones->MLh);
	free(Stones->Power);
	free(Stones);
}

void	SetStonesP(STONES *Stones)
{
	int Index;
	double	X;

	for(Index=0;Index<Stones->NoStones;Index++)
	{
		X = (Index+1) / (double)Stones->NoStones;

#ifndef NO_GSL
		Stones->Power[Index] = gsl_cdf_beta_Qinv(X, Stones->Alpha, Stones->Beta);
#else
		Stones->Power[Index] = 1.0 - X;
#endif
	}
}

STONES*	CratesStones(int NoS, int Sample, double Alpha, double Beta)
{
	STONES *Ret;
	int Index;

	Ret = (STONES*)SMalloc(sizeof(STONES));

	Ret->NoStones = NoS;

	Ret->MLh = (double*)SMalloc(sizeof(double) * Ret->NoStones);
	Ret->Power = (double*)SMalloc(sizeof(double) * Ret->NoStones);

	for(Index=0;Index<Ret->NoStones;Index++)
		Ret->MLh[Index] = 0.0;

	Ret->Alpha = Alpha;
	Ret->Beta = Beta;
	

	Ret->ItPerStone = Sample;
	Ret->ItStart	= -1;
	Ret->SampleFreq = 1;
	Ret->Started	= FALSE;
	
	SetStonesP(Ret);

	return Ret;
}

double	GetMLhStoneSum(STONES *Stones, int Pos)
{
	int Index;
	double Ret;
	Ret = 0;

	for(Index=0;Index<Pos;Index++)
		Ret += Stones->MLh[Index];

	return Ret;
}

void	NewStone(STONES *Stones, long long Itter, double Lh, int CStone, FILE *Out)
{
	if(CStone != 0)
	{
		Stones->MLh[CStone-1] = log(Stones->Sum / Stones->N) + Stones->Scalar;
//		fprintf(Out, "\t\t\t\t\t\t\t\tSS:\t%d\tMLh\t%f\n", CStone, Stones->MLh[CStone-1]);

		Stones->Length = Stones->Power[CStone-1] - Stones->Power[CStone];

		fprintf(Out, "%d\t%f\t%d\t%f\t%f\n", CStone-1, Stones->Power[CStone-1], Stones->N, Stones->MLh[CStone-1], GetMLhStoneSum(Stones, CStone));

		if(CStone == Stones->NoStones)
			fprintf(Out, "Log marginal likelihood:\t%f\n", GetMLhStoneSum(Stones, CStone));
	}
	else
		Stones->Length = 1.0 - Stones->Power[CStone];
		
	Stones->Scalar = Lh * Stones->Length;

	Stones->LastLh = Lh;
	Stones->Sum = 0;
	Stones->N = 0;

//	fprintf(Out, "\t\t\t\t\t\t\t\tSS:\t%d\tStartLh\t%f\tPower\t%f\n", CStone, Lh, Stones->Power[CStone]);
	fflush(Out);
}

int		SampleStone(STONES *Stones, long long StoneIt, int CStone)
{
	long long CIt;

	if(StoneIt % Stones->SampleFreq != 0)
		return FALSE;

	CIt = CStone * Stones->ItPerStone;
	CIt = StoneIt - CIt;

	if(CIt < Stones->ItPerStone * 0.25)
		return FALSE;

	return TRUE;
}

void	StoneItter(STONES *Stones, long long Itter, double Lh, FILE *Out)
{
	long long StoneIt;
	int CStone;

	if(Itter < Stones->ItStart)
		return;
	
	StoneIt = Itter - Stones->ItStart;
	CStone = (int)StoneIt / Stones->ItPerStone;

	if(StoneIt % Stones->ItPerStone == 0)
	{
		NewStone(Stones, Itter, Lh, CStone, Out);
	}
	
	if(SampleStone(Stones, StoneIt, CStone) == TRUE)
	{
		Stones->Sum += exp((Lh * Stones->Length) - Stones->Scalar);
		Stones->N++;
	}

}

int		StoneExit(STONES *Stones, long long Itters)
{
	if(Stones == NULL)
		return FALSE;

	if(Itters > Stones->ItStart + (Stones->ItPerStone * Stones->NoStones))
		return TRUE;

	return FALSE;
}