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
#include <string.h>
#include <stdlib.h>
#include <math.h>

//#define MATHMAT

#include "TypeDef.h"
#include "Trees.h"
#include "GenLib.h"
#include "Data.h"
#include "Likelihood.h"
#include "Matrix.h"
#include "LinAlg.h"
#include "Rates.h"
#include "CKappa.h"
#include "Contrasts.h"
#include "RandLib.h"
#include "RandDists.h"
#include "Part.h"


#ifdef BTOCL
	#include "btocl_continuous.h"
#endif

#ifdef BTLAPACK
	#include "btlapack_interface.h"
#endif

double	MLFindAlphaMeanRegTC(TREES* Trees, TREE *Tree);


void	InitEstData(OPTIONS *Opt, TREES *Trees)
{
	int		*TempEst;
	int		TIndex;
	int		SIndex;
	TAXA	*Taxa;

	TempEst = (int*)SMalloc(sizeof(int) * (Trees->NoSites+1));


	for(SIndex=0;SIndex<Trees->NoSites+1;SIndex++)
		TempEst[SIndex] = FALSE;
	Opt->NoEstDataSite = 0;

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];

		for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
		{
			if(Taxa->EstDataP[SIndex] == TRUE)
			{
				if(TempEst[SIndex] == FALSE)
				{
					TempEst[SIndex] = TRUE;
					Opt->NoEstDataSite++;
				}
			}
		}

		if(Taxa->EstDepData == TRUE)
		{
			if(TempEst[Trees->NoSites] == FALSE)
			{
				TempEst[Trees->NoSites] = TRUE;
				Opt->NoEstDataSite++;
			}
		}
	}

	if(Opt->NoEstDataSite == 0)
	{
		free(TempEst);
		Opt->EstDataSites = NULL;
		return;
	}

	Opt->EstDataSites = (int*)SMalloc(sizeof(int) * Opt->NoEstDataSite);

	TIndex=0;
	for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
	{
		if(TempEst[SIndex] == TRUE)
			Opt->EstDataSites[TIndex++] = SIndex;
	}

	if(TempEst[Trees->NoSites] == TRUE)
		Opt->EstDataSites[TIndex] = -1;

	free(TempEst);
}

void	RemoveDependantData(OPTIONS *Opt, TREES *Trees)
{
	int		Index;
	int		TIndex;
	TAXA*	Taxa;
	int		DepNo;

	DepNo = 0;

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];

		if(Taxa->EstDataP[DepNo] == TRUE)
		{
			Taxa->EstDepData = TRUE;
			for(Index=DepNo+1;Index<Trees->NoSites;Index++)
				Taxa->EstDataP[Index-1] = Taxa->EstDataP[Index];
		}

		Taxa->Dependant = Taxa->ConData[DepNo];

		for(Index=DepNo+1;Index<Trees->NoSites;Index++)
			Taxa->ConData[Index-1] = Taxa->ConData[Index];
	}

	Trees->NoSites--;
}

NODE	TaxaToNode(TREES* Trees, TREE *Tree, TAXA *Taxa)
{
	int	NIndex;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		if(Tree->NodeList[NIndex]->Tip == TRUE)
		{
		//	if(strcmp(Tree->NodeList[NIndex]->Taxa->Name, Taxa->Name)==0)
				if(Tree->NodeList[NIndex]->Taxa->No == Taxa->No)
					return Tree->NodeList[NIndex];
		}
	}

	printf("Error %s::%d cannot match taxa %s %d\n", __FILE__, __LINE__, Taxa->Name, Taxa->No);
	exit(0);
	return NULL;
}

double	GetTinyBL(TREES *Trees)
{
	int		TIndex;
	int		NIndex;
	NODE	Node;
	double	Ret;

	Ret = 100000;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
	{
		for(NIndex=0;NIndex<Trees->Tree[TIndex]->NoNodes;NIndex++)
		{
			Node = Trees->Tree[TIndex]->NodeList[NIndex];
			if((Node != Trees->Tree[TIndex]->Root) &&
				(Node->Length != 0))
			{
				if(Node->Length < Ret)
					Ret = Node->Length;
			}
		}
	}

	return Ret;
}

void	CheckZeroTaxaBL(TREES *Trees)
{
	int		TIndex;
	int		NIndex;
	NODE	Node;
	double	TinyBL;

	/* Get the smallest branch in the tree */
	TinyBL = GetTinyBL(Trees);

	/* Make it even smaller */
	TinyBL = TinyBL * 0.001;

	for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
	{
		for(NIndex=0;NIndex<Trees->Tree[TIndex]->NoNodes;NIndex++)
		{
			Node = Trees->Tree[TIndex]->NodeList[NIndex];

			if((Node->Tip == TRUE) && (Node != Trees->Tree[TIndex]->Root))
			{
				if(Node->Length < TinyBL)
					Node->Length = TinyBL;
			}
		}
	}
}

double	DistToRoot(NODE N, NODE Root)
{
	if(N == Root)
		return 0;

	return N->Length  + (DistToRoot(N->Ans, Root));
}

double	FindCoVar(TREES* Trees, TREE *Tree, int T1, int T2)
{
	int		Index;
	double	Ret;
	NODE	N1;
	NODE	N2;

	Ret = 0;

    for(Index=0;Index<Tree->NoNodes;Index++)
		Tree->NodeList[Index]->Visited = FALSE;

	N1 = TaxaToNode(Trees, Tree, Trees->Taxa[T1]);
	N2 = TaxaToNode(Trees, Tree, Trees->Taxa[T2]);

	while(N1!=Tree->Root)
	{
		N1 = N1->Ans;
		N1->Visited = TRUE;
	}

	do
	{
		N2 = N2->Ans;
	} while(N2->Visited != TRUE);

	Ret = DistToRoot(N2, Tree->Root);

	return Ret;
}

int		GetMapID(TREES *Trees, int ID)
{
	int Size, Index;

	Size = Trees->NoTaxa;
	for(Index=0;Index<Size;Index++)
	{
		if(Trees->Taxa[Index]->No == ID)
			return Index;
	}

	printf("%s::%d\tGetMapID ID = %d\n", __FILE__, __LINE__, ID);
	exit(0);
	return -1;
}

void	MapPartID(TREES* Trees, PART *Part, int *Map)
{
	int Index;

	for(Index=0;Index<Part->NoTaxa;Index++)
		Map[Index] = GetMapID(Trees, Part->Taxa[Index]);
}

void	RecCalcV(TREES* Trees, double **Mat, NODE N, PART *DiffPart)
{
	int x,y, XPos, YPos;

	double Dist;


	Dist = N->Ans->DistToRoot;

	GetPartDiff(N->Ans->Part, N->Part, DiffPart);

	for(x=0;x<N->Part->NoTaxa;x++)
	{
		XPos = N->Part->Taxa[x];
		for(y=0;y<DiffPart->NoTaxa;y++)
		{
			YPos = DiffPart->Taxa[y];
			Mat[YPos][XPos] = Dist;
			Mat[XPos][YPos] = Dist;
		}
	}

	if(N->Tip == TRUE)
	{
		XPos = GetMapID(Trees, N->Taxa->No);
		Mat[XPos][XPos] = N->DistToRoot;
		return;
	}

	for(x=0;x<N->NoNodes;x++)
		RecCalcV(Trees, Mat, N->NodeList[x], DiffPart);
}

void	CaclPVarCoVarRec(TREES* Trees, TREE *Tree, MATRIX *V)
{
	int		Index;
	NODE	N;
	PART	*CPart;

	CPart = CreatPart(Trees->NoTaxa);

	N = Tree->Root;
	for(Index=0;Index<N->NoNodes;Index++)
		RecCalcV(Trees, V->me, N->NodeList[Index], CPart);

	FreePart(CPart);
}

int		GetVPosNode(NODE N, TREES *Trees)
{
	if(N->Tip == TRUE)
		return GetMapID(Trees, N->Taxa->No);

	return GetVPosNode(N->NodeList[0], Trees);
}

void	SetNodeVPos(NODE N, TREE *Tree, TREES *Trees)
{
	int X,Y;

	if(Tree->Root == N)
		return;

	if(N->Tip == TRUE)
	{
		X = GetMapID(Trees, N->Taxa->No);
		N->VPosX = X;
		N->VPosY = X;
		return;
	}

	X = GetVPosNode(N->NodeList[0], Trees);
	Y = GetVPosNode(N->NodeList[1], Trees);

	N->VPosX = X;
	N->VPosY = Y;
}

void	SetNodesVPos(TREES *Trees, TREE *Tree)
{
	int NIndex;
	NODE N;

	for(NIndex=0;NIndex<Tree->NoNodes;NIndex++)
	{
		N = Tree->NodeList[NIndex];
		SetNodeVPos(N, Tree, Trees);
//		printf("%f\t%d\t", N->Length, N->NoNodes);
//		PrintPart(stdout, Trees, N->Part);
//		printf("\n");
	}

//	exit(0);
}

void	VToTree(MATRIX *V, TREE *Tree)
{
		int Index;
	NODE N;
	double AnsLen;


	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N != Tree->Root)
		{
			if(N->Ans == Tree->Root)
				N->Length = V->me[N->VPosX][N->VPosY];
			else
			{
				AnsLen = V->me[N->Ans->VPosX][N->Ans->VPosY];
				N->Length = V->me[N->VPosX][N->VPosY] - AnsLen;
			}
		}
	}
}

void	TreeToV(TREES* Trees, TREE *Tree, MATRIX *V)
{
	#ifdef ID_MATRIX
		SetIdentityMatrix(V);
		return;
	#endif

// Fast method to find V.
	CaclPVarCoVarRec(Trees, Tree, V);
//	PrintMatrix(Tree->ConVars->V, "V=", stdout);exit(0);
	return;

// Slow method to find V
/*	int		x,y;
	double	CoVar;
	NODE	N;
	double	*WV1, *WV2;
	MATRIX	*TMat;
	int	Ret;


	for(x=0;x<Trees->NoTaxa;x++)
	{
		for(y=x;y<Trees->NoTaxa;y++)
		{
			if(x!=y)
			{
				CoVar = FindCoVar(Trees, Tree, x,y);
				V->me[x][y] = CoVar;
				V->me[y][x] = CoVar;
			}
			else
			{
				N = TaxaToNode(Trees, Tree, Trees->Taxa[x]);
				CoVar = DistToRoot(N, Tree->Root);
				V->me[x][x] = CoVar;
			}
		}
	}
	*/
/*
	PrintMathematicaMatrix(Tree->VarCoVar, "Var Co var for a tree", stdout);
	PrintMatrix(Tree->VarCoVar, "Var Co var for a tree", stdout);

	WV1 = (double*)malloc(sizeof(double)*Tree->VarCoVar->NoOfCols);
	WV2 = (double*)malloc(sizeof(double)*Tree->VarCoVar->NoOfCols);

	TMat = AllocMatrix(Tree->VarCoVar->NoOfRows, Tree->VarCoVar->NoOfCols);

	Ret = InvertMatrix(Tree->VarCoVar->me, Tree->VarCoVar->NoOfCols, WV1, WV2, TMat->me);
	Ret = InvertMatrix(Tree->VarCoVar->me, Tree->VarCoVar->NoOfCols, WV1, (int*)WV2, TMat->me);

	printf("Ret\t%d\n", Ret);
	fflush(stdout);

	free(WV1);
	free(WV2);
	FreeMatrix(TMat);
	PrintMathematicaMatrix(TMat, "Inv vec", stdout);
*/
}

void	CalcZ(TREES* Trees, TREE *Tree, OPTIONS *Opt)
{
	int	SIndex, TIndex;
	int	ZPos;

	ZPos = 0;

	if(Opt->Model == M_CONTINUOUS_REG)
	{
		for(SIndex=0;SIndex<Trees->NoTaxa;SIndex++)
			Tree->ConVars->Z[SIndex] = Trees->Taxa[SIndex]->Dependant;
	}
	else
	{
		for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
		{
			for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++,ZPos++)
				Tree->ConVars->Z[ZPos] = Trees->Taxa[TIndex]->ConData[SIndex];
		}
	}
}

double	FindSum(TAXA** Taxa, int NoTaxa, int TraitNo)
{
	double	Ret=0;
	int		Index;

	for(Index=0;Index<NoTaxa;Index++)
		Ret += Taxa[Index]->ConData[TraitNo];

	return Ret;
}

double	FindSumSqu(TAXA** Taxa, int NoTaxa, int TraitNo)
{
	double	Ret;
	int		Index;

	Ret = 0;
	for(Index=0;Index<NoTaxa;Index++)
		Ret += (Taxa[Index]->ConData[TraitNo] * Taxa[Index]->ConData[TraitNo]);

	return Ret;
}

double	CalcVar(TREES* Trees, int TraitNo)
{
	double	Ret=0;
	double	Sum=0;
	double	SumSqu=0;

	Sum		= FindSum(Trees->Taxa, Trees->NoTaxa, TraitNo);
	SumSqu	= FindSumSqu(Trees->Taxa, Trees->NoTaxa, TraitNo);

	Ret = (Sum * Sum) / (double)Trees->NoTaxa;
	Ret = SumSqu - Ret;
	Ret = Ret / (double)(Trees->NoTaxa - 1);

	return Ret;
}

double	CalcCoVar(TREES* Trees, int T1, int T2)
{
	double	Ret=0;
	double	Sum1, Sum2;
	int		Index;

	Sum1 = FindSum(Trees->Taxa, Trees->NoTaxa, T1);
	Sum2 = FindSum(Trees->Taxa, Trees->NoTaxa, T2);

	for(Index=0;Index<Trees->NoTaxa;Index++)
		Ret = Ret + (Trees->Taxa[Index]->ConData[T1] * Trees->Taxa[Index]->ConData[T2]);

	Ret = Ret - ((Sum1 * Sum2) / (double)Trees->NoTaxa);
/*	Ret = Ret / (double)(Trees->NoTaxa - 1); */
	Ret = Ret / (double)(Trees->NoTaxa);

	return Ret;
}

void	CalcDVarCoVar(TREES* Trees, TREE *Tree)
{
	int		x,y;
	double	Val;

	for(x=0;x<Trees->NoSites;x++)
	{
		for(y=x;y<Trees->NoSites;y++)
		{
			if(x!=y)
			{
				Val = CalcCoVar(Trees, x, y);
				Tree->ConVars->Sigma->me[x][y] = Val;
				Tree->ConVars->Sigma->me[y][x] = Val;
			}
			else
			{
				Val = CalcVar(Trees, x);
				Tree->ConVars->Sigma->me[x][x] = Val;
			}
		}
	}
}

void	FreeConVar(CONVAR* ConVar, int NoTaxa)
{
	FreeMatrix(ConVar->V);
	FreeMatrix(ConVar->InvV);
	FreeMatrix(ConVar->Sigma);
	FreeMatrix(ConVar->InvSigma);
//	FreeMatrix(ConVar->KProd);
//	FreeMatrix(ConVar->InvKProd);

	if(ConVar->TrueV != NULL)
		FreeMatrix(ConVar->TrueV);

	if(ConVar->DepVect != NULL)
		free(ConVar->DepVect);

	free(ConVar->Alpha);
	free(ConVar->Z);
	free(ConVar->ZA);
	free(ConVar->ZATemp);

	free(ConVar->TVect1);
	free(ConVar->TVect2);
	free(ConVar->TVect3);
	free(ConVar->SVect);

	if(ConVar->TVT != NULL)
		FreeMatrix(ConVar->TVT);

	if(ConVar->TVTTemp != NULL)
		FreeMatrix(ConVar->TVTTemp);

	if(ConVar->InvXVX != NULL)
		FreeMatrix(ConVar->InvXVX);

	if(ConVar->Beta	!= NULL)
		free(ConVar->Beta);

	if(ConVar->MultiVarNormState != NULL)
		free(ConVar->MultiVarNormState);

	if(ConVar->MultiVarNormTemp != NULL)
		free(ConVar->MultiVarNormTemp);

	#ifdef BTOCL
	btocl_FreeConVar(ConVar);
	#endif

	free(ConVar);

	ConVar = NULL;
}

CONVAR*	AllocConVar(OPTIONS *Opt, TREES* Trees)
{
	CONVAR* Ret;
	int		Lager;

	Ret = (CONVAR*)SMalloc(sizeof(CONVAR));

	Ret->TVT		=	NULL;

	Ret->TVTTemp	=	NULL;
	Ret->InvXVX		=	NULL;
	Ret->Beta		=	NULL;

	if(Opt->Model == M_CONTINUOUS_DIR)
	{
		Ret->TVT	=	AllocMatrix(2, 2);
		Ret->TVTTemp=	AllocMatrix(2, Trees->NoTaxa);
	}

	if(Opt->Model == M_CONTINUOUS_REG)
	{
		if(Opt->AlphaZero == FALSE)
		{
			Ret->TVT	=	AllocMatrix(Trees->NoTaxa, Trees->NoSites + 1);
			Ret->TVTTemp=	AllocMatrix(Trees->NoTaxa, Trees->NoSites + 1);
		}
		else
		{
			Ret->TVT	=	AllocMatrix(Trees->NoTaxa, Trees->NoSites);
			Ret->TVTTemp=	AllocMatrix(Trees->NoTaxa, Trees->NoSites);
		}

		if(Opt->Analsis == ANALML)
			Ret->InvXVX	= AllocMatrix(Trees->NoSites+1, Trees->NoSites+1);
		else
			Ret->InvXVX = NULL;
	}

	Ret->V		= AllocMatrix(Trees->NoTaxa, Trees->NoTaxa);
	Ret->InvV	= AllocMatrix(Trees->NoTaxa, Trees->NoTaxa);

	if(Opt->InvertV == TRUE)
		Ret->TrueV = AllocMatrix(Trees->NoTaxa, Trees->NoTaxa);
	else
		Ret->TrueV = NULL;

	if(Opt->Model == M_CONTINUOUS_REG)
	{
//		Ret->KProd		= AllocMatrix(Trees->NoTaxa, Trees->NoTaxa);
//		Ret->InvKProd	= AllocMatrix(Trees->NoTaxa, Trees->NoTaxa);
		Ret->Sigma		= AllocMatrix(1, 1);
		Ret->InvSigma	= AllocMatrix(1, 1);
	}
	else
	{
//		Ret->KProd		= AllocMatrix(Trees->NoOfSites * Trees->NoTaxa, Trees->NoOfSites * Trees->NoTaxa);
//		Ret->InvKProd	= AllocMatrix(Trees->NoOfSites * Trees->NoTaxa, Trees->NoOfSites * Trees->NoTaxa);
		Ret->Sigma		= AllocMatrix(Trees->NoSites, Trees->NoSites);
		Ret->InvSigma	= AllocMatrix(Trees->NoSites, Trees->NoSites);
	}

	Ret->Beta	=	NULL;

	if((Opt->Model == M_CONTINUOUS_DIR) || (Opt->Model == M_CONTINUOUS_RR))
		Ret->Alpha	=	(double*)malloc(sizeof(double) * Trees->NoSites);
	else
		Ret->Alpha	=	(double*)malloc(sizeof(double) * 1);

	if(Ret->Alpha == NULL)
		MallocErr();

	if((Opt->Model == M_CONTINUOUS_DIR) || (Opt->Model == M_CONTINUOUS_REG))
	{
		Ret->Beta = (double*)malloc(sizeof(double) * Trees->NoSites);
		if(Ret->Beta == NULL)
			MallocErr();
	}

	if(Opt->Model == M_CONTINUOUS_REG)
	{
		Ret->Z		=	(double*)SMalloc(sizeof(double) * Trees->NoTaxa);
		Ret->ZA		=	(double*)SMalloc(sizeof(double) * Trees->NoTaxa);
		Ret->ZATemp	=	(double*)SMalloc(sizeof(double) * Trees->NoTaxa);
		Ret->DepVect=	(double*)SMalloc(sizeof(double) * Trees->NoTaxa);
	}
	else
	{
		Ret->Z		=	(double*)SMalloc(sizeof(double) * Trees->NoSites * Trees->NoTaxa);
		Ret->ZA		=	(double*)SMalloc(sizeof(double) * Trees->NoSites * Trees->NoTaxa);
		Ret->ZATemp	=	(double*)SMalloc(sizeof(double) * Trees->NoSites * Trees->NoTaxa);
		Ret->DepVect=	NULL;
	}

	if(Trees->NoTaxa > Trees->NoSites)
		Lager = Trees->NoTaxa;
	else
		Lager = Trees->NoSites;

	Ret->TVect1	=	(double*)SMalloc(sizeof(double) * Lager);
	Ret->TVect2	=	(double*)SMalloc(sizeof(double) * Lager);
	Ret->TVect3	=	(double*)SMalloc(sizeof(double) * Lager);
	Ret->SVect	=	(double*)SMalloc(sizeof(double) * Lager);

	Ret->MultiVarNormState = NULL;
	Ret->MultiVarNormTemp = NULL;

	Ret->LogDetOfSigma	= 1;
	Ret->LogDetOfV		= 1;

	#ifdef BTOCL
	btocl_AllocConVar(Ret,Trees);
	#endif

	return Ret;
}


void	FindTVT(TREES* Trees, TREE *Tree, int AlphaZero)
{
	int			Row, Col, Index;
	double		Temp;
	MATRIX*		TMat;
	CONVAR*		ConVars;

	ConVars = Tree->ConVars;

	TMat = AllocMatrix(2, Trees->NoTaxa);

	for(Col=0;Col<Trees->NoTaxa;Col++)
	{
		if(AlphaZero == FALSE)
			ConVars->TVTTemp->me[0][Col] = 1;
		else
			ConVars->TVTTemp->me[0][Col] = 0;

		ConVars->TVTTemp->me[1][Col] = Tree->ConVars->V->me[Col][Col];
	}

	for(Row=0;Row<2;Row++)
	{
		for(Col=0;Col<Trees->NoTaxa;Col++)
		{
			Temp = 0;
			for(Index=0;Index<Trees->NoTaxa;Index++)
			{
				Temp += ConVars->TVTTemp->me[Row][Index] * ConVars->InvV->me[Index][Col];
			}
			TMat->me[Row][Col] = Temp;
		}
	}

	Temp = 0;
	for(Index=0;Index<Trees->NoTaxa;Index++)
		Temp += ConVars->TVTTemp->me[0][Index] * TMat->me[0][Index];

	ConVars->TVT->me[0][0] = Temp;

	Temp = 0;
	for(Index=0;Index<Trees->NoTaxa;Index++)
		Temp += ConVars->TVTTemp->me[1][Index] * TMat->me[0][Index];

	ConVars->TVT->me[0][1] = Temp;

	Temp = 0;
	for(Index=0;Index<Trees->NoTaxa;Index++)
		Temp += ConVars->TVTTemp->me[0][Index] * TMat->me[1][Index];

	ConVars->TVT->me[1][0] = Temp;

	Temp = 0;
	for(Index=0;Index<Trees->NoTaxa;Index++)
		Temp += ConVars->TVTTemp->me[1][Index] * TMat->me[1][Index];

	ConVars->TVT->me[1][1] = Temp;

	FreeMatrix(TMat);

	if(AlphaZero == TRUE)
	{
		ConVars->TVT->me[1][1] = 1.0 / ConVars->TVT->me[1][1];
		return;
	}

	Temp = ConVars->TVT->me[0][0];
	ConVars->TVT->me[0][0] = ConVars->TVT->me[1][1];
	ConVars->TVT->me[1][1] = Temp;

	Temp = 1.0 / (ConVars->TVT->me[0][0] * ConVars->TVT->me[1][1] - ConVars->TVT->me[1][0] * ConVars->TVT->me[0][1]);

	ConVars->TVT->me[0][0] = Temp * ConVars->TVT->me[0][0];
	ConVars->TVT->me[0][1] = Temp * -ConVars->TVT->me[0][1];
	ConVars->TVT->me[1][0] = Temp * -ConVars->TVT->me[1][0];
	ConVars->TVT->me[1][1] = Temp * ConVars->TVT->me[1][1];
}

void	FindMLRagVals(TREES* Trees, TREE *Tree, OPTIONS *Opt)
{
	TEMPCONVAR*	TempCon;
	int		x,y, Err;
	CONVAR	*CV;


	CV = Tree->ConVars;
	TempCon = Trees->TempConVars;

	for(x=0;x<Trees->NoTaxa;x++)
		TempCon->Y[x] = Trees->Taxa[x]->Dependant;


	if(Opt->AlphaZero == FALSE)
	{
		for(x=0;x<Trees->NoTaxa;x++)
			TempCon->X->me[x][0] = 1;

		for(x=1;x<Trees->NoSites+1;x++)
		{
			for(y=0;y<Trees->NoTaxa;y++)
				TempCon->X->me[y][x] = Trees->Taxa[y]->ConData[x-1];
		}
	}
	else
	{
		for(x=0;x<Trees->NoSites;x++)
		{
			for(y=0;y<Trees->NoTaxa;y++)
				TempCon->X->me[y][x] = Trees->Taxa[y]->ConData[x];
		}
	}



//	PrintMatrix(CV->InvV, "M=", stdout);exit(0);

	/* Calc X'.InvV.Y */
	MatrixByVectMult(CV->InvV, TempCon->Y, CV->TVect1);
	VectByMatrixMult(CV->TVect1, TempCon->X, CV->TVect2);
	/* Now X'.InvV.Y is in TVect2 */

	/* Calc d */
	MatrixMult(CV->InvV, TempCon->X, CV->TVTTemp);

	Transpose(TempCon->X, TempCon->TranX);
	MatrixMult(TempCon->TranX, CV->TVTTemp, TempCon->NX);

//	PrintMatrix(TempCon->NX, "NX=",stdout);

	Err = InvertMatrix(TempCon->NX->me, TempCon->NX->NoOfCols, CV->TVect1, (int*)CV->TVect3, CV->InvXVX->me);

	if(Err != FALSE)
	{
		printf("Error inverting sigma matrix, matrix is singular\n");
		exit(0);
	}

	MatrixByVectMult(CV->InvXVX, CV->TVect2, CV->TVect1);

	if(Opt->AlphaZero == FALSE)
	{
		CV->Alpha[0] = CV->TVect1[0];
		for(x=1;x<Trees->NoSites+1;x++)
			CV->Beta[x-1] = CV->TVect1[x];
	}
	else
	{
		CV->Alpha[0] = 0;
		for(x=0;x<Trees->NoSites;x++)
			CV->Beta[x] = CV->TVect1[x];
	}

	if(Opt->TestCorrel == FALSE)
	{
		CV->Alpha[0] = MLFindAlphaMeanRegTC(Trees, Tree);

		for(x=0;x<Trees->NoSites;x++)
			CV->Beta[x] = 0;
	}

/*
	printf("ML Vals\n");
	for(x=0;x<Trees->NoOfSites+1;x++)
		printf("%f\t", CV->TVect1[x]);
	printf("\nML Done\n");
*/

/*	To Check With Mathematica.
	PrintMathematicaMatrix(X, "X = ", stdout);
	PrintMathematicaMatrix(CV->InvV, "InvV = ", stdout);
	PrintMathematicaVect(Y, Trees->NoTaxa, "Y = ", stdout);
	printf("B1 = Inverse[Transpose[X].InvV.X]\n");
	printf("B2 = Transpose[X].InvV.Y\n");
	printf("B = B1.B2\n");

	exit(0);
*/
}

void	MLFindAlphaBeta(TREES* Trees, TREE *Tree, int Site, int AlphaZero)
{
	int			Row, Col, Index;
	double		Temp;
	double		TempTVX[2];
	MATRIX*		TMat;
	CONVAR*		ConVars;

	ConVars = Tree->ConVars;


	TMat = AllocMatrix(2, Trees->NoTaxa);

	for(Row=0;Row<2;Row++)
	{
		for(Col=0;Col<Trees->NoTaxa;Col++)
		{
			Temp = 0;
			for(Index=0;Index<Trees->NoTaxa;Index++)
			{
				Temp += ConVars->TVTTemp->me[Row][Index] * ConVars->InvV->me[Index][Col];
			}
			TMat->me[Row][Col] = Temp;
		}
	}

	Temp = 0;
	for(Index=0;Index<Trees->NoTaxa;Index++)
		Temp += Trees->Taxa[Index]->ConData[Site] * TMat->me[0][Index];

	TempTVX[0] = Temp;

	Temp = 0;
	for(Index=0;Index<Trees->NoTaxa;Index++)
		Temp += Trees->Taxa[Index]->ConData[Site] * TMat->me[1][Index];

	TempTVX[1] = Temp;

	if(AlphaZero == TRUE)
	{
		Tree->ConVars->Alpha[Site] = 0;
		Tree->ConVars->Beta[Site] = ConVars->TVT->me[1][1] * TempTVX[1];

	}
	else
	{
		Tree->ConVars->Alpha[Site]	= (ConVars->TVT->me[0][0] * TempTVX[0]) + (ConVars->TVT->me[0][1] * TempTVX[1]);
		Tree->ConVars->Beta[Site]	= (ConVars->TVT->me[1][0] * TempTVX[0]) + (ConVars->TVT->me[1][1] * TempTVX[1]);
	}

	FreeMatrix(TMat);
}

double	MLFindAlphaReg(TREES* Trees, TREE *Tree, double *Data)
{
	double	P1=0;
	double	P2;
	double	ColTemp;
	int		x,y;

	for(x=0;x<Trees->NoTaxa;x++)
	{
		for(y=x+1;y<Trees->NoTaxa;y++)
			P1 = P1 + (2 * Tree->ConVars->InvV->me[x][y]);
	}

	for(x=0;x<Trees->NoTaxa;x++)
		P1 = P1 + Tree->ConVars->InvV->me[x][x];

	P1 = 1 / P1;

	P2 = 0;
	for(y=0;y<Trees->NoTaxa;y++)
	{
		ColTemp = 0;
		for(x=0;x<Trees->NoTaxa;x++)
			ColTemp = ColTemp + Tree->ConVars->InvV->me[x][y];
		P2 += ColTemp * Data[y];
	}

	return P1 * P2;
}

/*

double	MLFindAlphaReg(TREES* Trees, TREE *Tree)
{
	double	P1=0;
	double	P2;
	double	ColTemp;
	int		x,y;

	for(x=0;x<Trees->NoTaxa;x++)
	{
		for(y=x+1;y<Trees->NoTaxa;y++)
			P1 = P1 + (2 * Tree->ConVars->InvV->me[x][y]);
	}

	for(x=0;x<Trees->NoTaxa;x++)
		P1 = P1 + Tree->ConVars->InvV->me[x][x];

	P1 = 1 / P1;

	P2 = 0;
	for(y=0;y<Trees->NoTaxa;y++)
	{
		ColTemp = 0;
		for(x=0;x<Trees->NoTaxa;x++)
			ColTemp = ColTemp + Tree->ConVars->InvV->me[x][y];
		P2 += ColTemp * Trees->Taxa[y].Dependant;
	}

	return P1 * P2;
}
*/
double	MLFindAlphaMean(TREES* Trees, TREE *Tree, int Site)
{
	double	P1;
	double	P2;
	double	ColTemp;
	int		x,y;

	P1 = 0;
	for(x=0;x<Trees->NoTaxa;x++)
	{
		for(y=x+1;y<Trees->NoTaxa;y++)
			P1 = P1 + (2 * Tree->ConVars->InvV->me[x][y]);
	}

	for(x=0;x<Trees->NoTaxa;x++)
		P1 = P1 + Tree->ConVars->InvV->me[x][x];

	P1 = 1 / P1;

	P2 = 0;
	for(y=0;y<Trees->NoTaxa;y++)
	{
		ColTemp = 0;
		for(x=0;x<Trees->NoTaxa;x++)
			ColTemp = ColTemp + Tree->ConVars->InvV->me[x][y];
		P2 += ColTemp * Trees->Taxa[y]->ConData[Site];
	//	P2 += ColTemp * Trees->Taxa[y]->Dependant;
	}

	return P1 * P2;
}

double	MLFindAlphaMeanRegTC(TREES* Trees, TREE *Tree)
{
	double	P1;
	double	P2;
	double	ColTemp;
	int		x,y;

	P1 = 0;
	for(x=0;x<Trees->NoTaxa;x++)
	{
		for(y=x+1;y<Trees->NoTaxa;y++)
			P1 = P1 + (2 * Tree->ConVars->InvV->me[x][y]);
	}

	for(x=0;x<Trees->NoTaxa;x++)
		P1 = P1 + Tree->ConVars->InvV->me[x][x];

	P1 = 1 / P1;

	P2 = 0;
	for(y=0;y<Trees->NoTaxa;y++)
	{
		ColTemp = 0;
		for(x=0;x<Trees->NoTaxa;x++)
			ColTemp = ColTemp + Tree->ConVars->InvV->me[x][y];
		P2 += ColTemp * Trees->Taxa[y]->Dependant;
	}

	return P1 * P2;
}


double	FindMLVarMatic(TREES* Trees, TREE *Tree, int SLS)
{
	double	Ret;
	int		x,y;
	CONVAR	*CV;

	CV = Tree->ConVars;

#ifdef OPENMP_THR
	#pragma omp parallel for private(x, Ret) num_threads(4)
#endif
	for(y=0;y<Trees->NoTaxa;y++)
	{
		Ret = 0;
		for(x=0;x<Trees->NoTaxa;x++)
			Ret = Ret + (CV->TVect1[x] * CV->InvV->me[x][y]);

		CV->TVect2[y] = Ret;
	}

	Ret = 0;
	for(x=0;x<Trees->NoTaxa;x++)
		Ret = Ret + (CV->TVect2[x] * CV->TVect3[x]);

	// Use sum of least squares, should not be used.
//	if(SLS == TRUE)
//		return Ret * (1.0/(Trees->NoTaxa - (Trees->NoOfSites+1)));

	Ret = Ret * (1.0/Trees->NoTaxa);

	return Ret;
}

double	FindMLVar(TREES* Trees, TREE *Tree, int Site1, double Alpha1, double Beta1, int Site2, double Alpha2, double Beta2)
{
	int		x;
	CONVAR	*CV;

	CV = Tree->ConVars;

	for(x=0;x<Trees->NoTaxa;x++)
	{
		CV->TVect1[x] = Trees->Taxa[x]->ConData[Site1]  - (Alpha1 + (Beta1 * CV->V->me[x][x]));
		CV->TVect3[x] = Trees->Taxa[x]->ConData[Site2]  - (Alpha2 + (Beta2 * CV->V->me[x][x]));

		CV->TVect2[x] = 0;
	}

	return FindMLVarMatic(Trees, Tree, FALSE);
}


double	FindMLRegVar(TREES* Trees, TREE *Tree)
{
	int		TIndex;
	int		x;
	CONVAR	*CV;
	double	Reg;

	CV = Tree->ConVars;


	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Reg = CV->Alpha[0];
		for(x=0;x<Trees->NoSites;x++)
			Reg += CV->Beta[x] * Trees->Taxa[TIndex]->ConData[x];

		CV->TVect1[TIndex] = Trees->Taxa[TIndex]->Dependant  - Reg;
		CV->TVect3[TIndex] = CV->TVect1[TIndex];
		CV->TVect2[TIndex] = 0;

	}

	Reg = FindMLVarMatic(Trees, Tree, TRUE);

	return Reg;
}

void	CalcSigma(OPTIONS *Opt, TREES* Trees, TREE *Tree, double* Means, double* Beta)
{
	int		x,y;
	double	Val;

	if((Opt->Analsis == ANALML) && (Opt->Model == M_CONTINUOUS_DIR))
		FindTVT(Trees, Tree, Opt->AlphaZero);

	if((Opt->Analsis == ANALML) && (Opt->Model == M_CONTINUOUS_REG))
		FindMLRagVals(Trees, Tree, Opt);

	for(x=0;x<Trees->NoSites;x++)
	{
		if(Opt->Analsis == ANALML)
		{
			if(Opt->Model == M_CONTINUOUS_DIR)
				MLFindAlphaBeta(Trees, Tree, x, Opt->AlphaZero);

			if(Opt->Model == M_CONTINUOUS_RR)
			{
				if(Opt->AlphaZero == FALSE)
					Tree->ConVars->Alpha[x] = MLFindAlphaMean(Trees, Tree, x);
				else
					Tree->ConVars->Alpha[x] = 0.0;
			}
		}
		else
		{
			if((Opt->Model == M_CONTINUOUS_RR) || (Opt->Model == M_CONTINUOUS_DIR))
				Tree->ConVars->Alpha[x] = Means[x];
			else
				Tree->ConVars->Alpha[0] = Means[0];

			if((Opt->Model == M_CONTINUOUS_DIR) || (Opt->Model == M_CONTINUOUS_REG))
				Tree->ConVars->Beta[x] = Beta[x];
		}
	}

	if(Opt->Model == M_CONTINUOUS_REG)
	{
		Tree->ConVars->Sigma->me[0][0] = FindMLRegVar(Trees, Tree);
		return;
	}

	for(x=0;x<Trees->NoSites;x++)
	{
		for(y=x;y<Trees->NoSites;y++)
		{
			if(Opt->TestCorrel == TRUE)
			{
				if(Opt->Model == M_CONTINUOUS_RR)
					Val = FindMLVar(Trees, Tree, x, Tree->ConVars->Alpha[x], 0, y, Tree->ConVars->Alpha[y], 0);
				else
					Val = FindMLVar(Trees, Tree, x, Tree->ConVars->Alpha[x], Tree->ConVars->Beta[x], y, Tree->ConVars->Alpha[y], Tree->ConVars->Beta[y]);

				Tree->ConVars->Sigma->me[x][y] = Val;
				Tree->ConVars->Sigma->me[y][x] = Val;
			}
			else
			{
				if(x==y)
				{
					if(Opt->Model == M_CONTINUOUS_RR)
						Tree->ConVars->Sigma->me[x][x] = FindMLVar(Trees, Tree, x, Tree->ConVars->Alpha[x], 0, y, Tree->ConVars->Alpha[y], 0);
					else
						Tree->ConVars->Sigma->me[x][x] = FindMLVar(Trees, Tree, x, Tree->ConVars->Alpha[x], Tree->ConVars->Beta[x], y, Tree->ConVars->Alpha[y], Tree->ConVars->Beta[y]);
				}
				else
				{
					Tree->ConVars->Sigma->me[x][y] = 0;
					Tree->ConVars->Sigma->me[y][x] = 0;
				}
			}
		}
	}
}

double	CalcZAReg(TREES* Trees, TREE* Tree, int TaxaNo)
{
	double	Ret;
	int		Index;
	CONVAR	*CV;
	TAXA	*Taxa;

	CV	= Tree->ConVars;
	Taxa= Trees->Taxa[TaxaNo];

	Ret = CV->Alpha[0];
	for(Index=0;Index<Trees->NoSites;Index++)
		Ret += CV->Beta[Index] * Taxa->ConData[Index];

	return Ret;
}

void	RegCalcZAlpha(TREES* Trees, TREE *Tree)
{
	int	TIndex;
	CONVAR	*CV;

	CV = Tree->ConVars;

	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
		CV->ZA[TIndex] = CV->Z[TIndex] - CalcZAReg(Trees, Tree, TIndex);

}

void	CalcZAlpha(TREES* Trees, TREE *Tree, MODEL Model)
{
	int	SIndex, TIndex;
	int	ZPos;

	ZPos = 0;

	if(Model == M_CONTINUOUS_REG)
	{
		RegCalcZAlpha(Trees, Tree);
		return;
	}

	for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
	{
		for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++,ZPos++)
		{
			switch(Model)
			{
				case M_CONTINUOUS_RR:
					Tree->ConVars->ZA[ZPos] = Tree->ConVars->Z[ZPos] - Tree->ConVars->Alpha[SIndex];
				break;

				case M_CONTINUOUS_DIR:
					Tree->ConVars->ZA[ZPos] = Tree->ConVars->Z[ZPos] - (Tree->ConVars->Alpha[SIndex] + (Tree->ConVars->Beta[SIndex] * Tree->ConVars->V->me[TIndex][TIndex]));
				break;


				// Keep CLang happy
				default:
					break;
			}
		}
	}
}


double	FindDet(TREES* Trees, TREE *Tree, OPTIONS *Opt)
{
	if(Opt->Model == M_CONTINUOUS_REG)
	{
		return (Tree->ConVars->LogDetOfV) + ((double)Trees->NoTaxa * Tree->ConVars->LogDetOfSigma);
	}

	return ((double)Trees->NoSites * Tree->ConVars->LogDetOfV) + ((double)Trees->NoTaxa * Tree->ConVars->LogDetOfSigma);
}

/*
<< LinearAlgebra`MatrixManipulation`;
BlockMatrix[Outer[Times, S, V]]
*/
/*
void	CalcOU(MATRIX *V,  double Alpha)
{
	int x,y;
	double T, T1, T2;
	double Scale;

	if(Alpha <= MIN_OU)
		Alpha = MIN_OU;

	Scale = 1.0 / (2.0 * Alpha);

	for(x=0;x<V->NoOfCols;x++)
	{
		for(y=x+1;y<V->NoOfRows;y++)
		{
			T = 0.5 * (V->me[x][x] + V->me[y][y]);

//			if(V->me[x][x] != V->me[y][y])
//				printf("ee\n");

			if(T < V->me[x][y])
				T = V->me[x][y];

			// T may be < V->me[x][y], for non ultrametic trees. Possibly an issue, don't know.
			T1 = exp(-2.0 * Alpha * (T - V->me[x][y]));
			T2 = 1.0 - exp(-2.0 * Alpha * V->me[x][y]);

			T1 = T1 * T2;

	//		for the sig^2 / 2 alpha multiplyer
			T1 = T1 * Scale;

			V->me[x][y] = T1;
			V->me[y][x] = T1;
		}
	}

	for(x=0;x<V->NoOfCols;x++)
	{
		T = 1.0 - exp(-2 * Alpha * V->me[x][x]);

//		for the sig^2 / 2 alpha multiplyer
		T = T * Scale;
		V->me[x][x] = T;
	}
}
*/

double FindT(MATRIX *V)
{
	int Index;
	double Ret;

	Ret = V->me[0][0];
	for(Index=1;Index<V->NoOfRows;Index++)
	{
		if(V->me[Index][Index] > Ret)
			Ret = V->me[Index][Index];
	}

	return Ret;
}
/*
void	CalcOU(TREES *Trees, TREE *Tree, MATRIX *V,  double Alpha)
{
	int x,y;
	double T, T1, T2;
	double Scale;

	if(Alpha <= MIN_OU)
		Alpha = MIN_OU;

	Scale = 1.0 / (2.0 * Alpha);
	T = FindT(V);

	for(x=0;x<V->NoOfCols;x++)
	{
		for(y=x+1;y<V->NoOfRows;y++)
		{
			// T may be < V->me[x][y], for non ultrametic trees. Possibly an issue, don't know.
			T1 = exp(-2.0 * Alpha * (T - V->me[x][y]));
			T2 = 1.0 - exp(-2.0 * Alpha * V->me[x][y]);

			T1 = T1 * T2;

	//		for the sig^2 / 2 alpha multiplyer
			T1 = T1 * Scale;

			V->me[x][y] = T1;
			V->me[y][x] = T1;
		}
	}

	for(x=0;x<V->NoOfCols;x++)
	{
		T1 = 1.0 - exp(-2 * Alpha * V->me[x][x]);
//		T1 = 1.0 - exp(-2 * Alpha * T);

//		for the sig^2 / 2 alpha multiplyer
		T1 = T1 * Scale;
		V->me[x][x] = T1;
	}
}
*/


double	GetOUSharedPath(TREES *Trees, TREE *Tree, double Shared, int X, int Y)
{
	int Index, IDX, IDY;
	double Ret, DX, DY;
	NODE N;

	IDX = Trees->Taxa[X]->No;
	IDY = Trees->Taxa[Y]->No;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N  = Tree->NodeList[Index];
		if(N->TipID == IDX)
			DX = N->DistToRoot;

		if(N->TipID == IDY)
			DY = N->DistToRoot;
	}

	Ret = (DX - Shared) + (DY - Shared);

	return Ret;
}
/*
// Have to be careful.
// Function uses both tree BL and V, nothing can be applied beforehand.
void	CalcOU(TREES *Trees, TREE *Tree, MATRIX *V,  double Alpha)
{
	int x,y;
	double RPath, T, T1, T2;
	double Scale;

//		PrintMatrix(V, "V=", stdout);exit(0);
	//	Alpha = 0.1;

	if(Alpha <= MIN_OU)
		Alpha = MIN_OU;

	Scale = 1.0 / (2.0 * Alpha);
	T = FindT(V);

	// T may be < V->me[x][y], for non ultrametic trees. Possibly an issue, don't know.
	for(x=0;x<V->NoOfCols;x++)
	{
		for(y=x;y<V->NoOfRows;y++)
		{
			RPath = GetOUSharedPath(Trees, Tree, V->me[x][y], x, y);
			//T1 = exp(-2.0 * Alpha * (T - V->me[x][y]));
			T1 = exp(-Alpha * RPath);
			T2 = 1.0 - exp(-2.0 * Alpha * V->me[x][y]);

			T1 = T1 * T2;

			// For the sig^2 / 2 alpha multiplyer
			T1 = T1 * Scale;

			V->me[x][y] = T1;
			V->me[y][x] = T1;
		}
	}

//	VToTree(V, Tree);
//	SaveTrees("vtree.trees", Trees);
	PrintMatrix(V, "V=", stdout);exit(0);

}
*/


double	CalcOUCell(double Alpha, double Scale, double XPath, double YPath, double Shared)
{
	double RPath, T1, T2;

	RPath = (XPath - Shared) + (YPath - Shared);

	T1 = exp(-Alpha * RPath);
	T2 = 1.0 - exp(-2.0 * Alpha * Shared);

	T1 = T1 * T2;

	return T1 * Scale;
}

void	CalcOU(TREES *Trees, TREE *Tree, MATRIX *V,  double Alpha)
{
	int x,y;
	double Scale;

	//	PrintMatrix(V, "V=", stdout);exit(0);
	//	Alpha = 0.1;

	if(Alpha <= MIN_OU)
		Alpha = MIN_OU;

	Scale = 1.0 / (2.0 * Alpha);

	// T may be < V->me[x][y], for non ultrametic trees. Possibly an issue, don't know.
	// Off diag must be done first.
	for(x=0;x<V->NoOfCols;x++)
	{
		for(y=x+1;y<V->NoOfRows;y++)
		{
			V->me[x][y] = CalcOUCell(Alpha, Scale, V->me[x][x], V->me[y][y], V->me[x][y]);
			V->me[y][x] = V->me[x][y];
		}
	}

	// do diag
	for(x=0;x<V->NoOfCols;x++)
		V->me[x][x] = CalcOUCell(Alpha, Scale, V->me[x][x], V->me[x][x], V->me[x][x]);

	//	VToTree(V, Tree);
	//	SaveTrees("vtree.trees", Trees);
//	PrintMatrix(V, "V=", stdout);exit(0);
}


/* Pre Non Ultramtrix
void	CalcOU(MATRIX *V,  double Alpha)
{
	int x,y;
	double T, T1, T2;
	double Scale;

//	PrintMatrix(V, "V=", stdout);exit(0);

//	Alpha = 0.1;

	if(Alpha <= MIN_OU)
		Alpha = MIN_OU;

	Scale = 1.0 / (2.0 * Alpha);
	T = FindT(V);


	// T may be < V->me[x][y], for non ultrametic trees. Possibly an issue, don't know.
	for(x=0;x<V->NoOfCols;x++)
	{
		for(y=x;y<V->NoOfRows;y++)
		{
			T1 = exp(-2.0 * Alpha * (T - V->me[x][y]));
			T2 = 1.0 - exp(-2.0 * Alpha * V->me[x][y]);

			T1 = T1 * T2;

	//		for the sig^2 / 2 alpha multiplyer
			T1 = T1 * Scale;
			V->me[x][y] = T1;
			V->me[y][x] = T1;
		}
	}
}
*/


/*
void	CalcDelta(MATRIX *V, double Delta)
{
	int x,y;

	for(x=0;x<V->NoOfCols;x++)
	{
		for(y=0;y<V->NoOfRows;y++)
		{
			V->me[x][y]= pow(V->me[x][y], Delta);
		}
	}
}
*/

void	CalcDelta(MATRIX *V, double Delta)
{
	int x,y;

	for(x=0;x<V->NoOfCols;x++)
	{
		V->me[x][x] = pow(V->me[x][x], Delta);
		for(y=x+1;y<V->NoOfRows;y++)
		{
			V->me[x][y] = pow(V->me[x][y], Delta);
			V->me[y][x] = V->me[x][y];
		}
	}
}

void	CalcLambda(MATRIX *V, double Lambda)
{
	int x,y;

	for(x=0;x<V->NoOfCols;x++)
	{
		for(y=x;y<V->NoOfRows;y++)
		{
			if(x!=y)
			{
				V->me[x][y] = V->me[x][y] * Lambda;
				V->me[y][x] = V->me[x][y];
			}
		}
	}
}
#ifdef akdlkljkajlk
MATRIX*	FindRegVar(TREES *Trees, RATES* Rates)
{
	MATRIX	*Ret;
	static MATRIX	*XT=NULL;
	static MATRIX	*X=NULL;
	static MATRIX	*TempV1;
	static MATRIX	*TempV2;
	TAXA	*Taxa;
	TREE	*Tree;
	CONVAR	*CV;
	int		x,y;

	Tree = &Trees->Tree[Rates->TreeNo];
	CV = Tree->ConVars;

	if(X == NULL)
	{
		X		= AllocMatrix(Trees->NoTaxa, Trees->NoOfSites+1);
		XT		= AllocMatrix(Trees->NoOfSites+1, Trees->NoTaxa);
		TempV1	= AllocMatrix(Trees->NoOfSites+1, Trees->NoTaxa);
		TempV2	= AllocMatrix(Trees->NoOfSites+1, Trees->NoOfSites+1);

		for(x=0;x<Trees->NoTaxa;x++)
		{
			Taxa = &Trees->Taxa[x];
			X->me[x][0] = 1;
			for(y=0;y<Trees->NoOfSites;y++)
				X->me[x][y+1] = Taxa->ConData[y];
		}

		Transpose(X, XT);
	}

	Ret = AllocMatrix(Trees->NoOfSites+1, Trees->NoOfSites+1);

	/* Do
		Sig*Inverse[Transpose[X].InvV.X]
	*/
	MatrixMult(XT, CV->InvV, TempV1);
	MatrixMult(TempV1, X, TempV2);

	InvertMatrix(TempV2->me, Trees->NoOfSites+1, CV->TVect1,(int*)CV->TVect2, Ret->me);

	ScaleMatrix(Ret, CV->Sigma->me[0][0]);

/*	printf("\n");
	PrintMathematicaMatrix(X, "X=", stdout);
	PrintMathematicaMatrix(CV->InvV, "InvV=", stdout);
	printf("Sig=%f;\n", CV->Sigma->me[0][0]);

	printf("Sqrt[Sig*Inverse[Transpose[X].InvV.X]]\n");
	exit(0);
*/

	return Ret;
}
#endif

MATRIX*	FindRegVar(TREES *Trees, RATES* Rates, int AlphaZero)
{
	MATRIX	*Ret;
	TEMPCONVAR	*TempCon;
	int		Err;
	TREE	*Tree;
	CONVAR	*CV;
	int		Size;

	TempCon = Trees->TempConVars;

	Tree = Trees->Tree[Rates->TreeNo];
	CV = Tree->ConVars;

	if(AlphaZero == FALSE)
		Size = Trees->NoSites+1;
	else
		Size = Trees->NoSites;

	Ret = AllocMatrix(Size, Size);

	MatrixMult(TempCon->XT, CV->InvV, TempCon->TempV1);
	MatrixMult(TempCon->TempV1, TempCon->RVX, TempCon->TempV2);

	Err = InvertMatrix(TempCon->TempV2->me, Size, CV->TVect1,(int*)CV->TVect2, Ret->me);
	if(Err != FALSE)
	{
		printf("Error inverting matrix, matrix is singular, %s::%d\n", __FILE__, __LINE__);
		exit(0);
	}

	ScaleMatrix(Ret, CV->Sigma->me[0][0]);

/*	printf("\n");
	PrintMathematicaMatrix(X, "X=", stdout);
	PrintMathematicaMatrix(CV->InvV, "InvV=", stdout);
	printf("Sig=%f;\n", CV->Sigma->me[0][0]);

	printf("Sqrt[Sig*Inverse[Transpose[X].InvV.X]]\n");
	exit(0);
*/

	return Ret;
}

void	SetEstData(TREES *Trees, RATES* Rates)
{
	int		TIndex;
	TAXA*	Taxa;
	int		SIndex;
	int		RIndex;

	RIndex=0;
	for(TIndex=0;TIndex<Trees->NoTaxa;TIndex++)
	{
		Taxa = Trees->Taxa[TIndex];

		if(Taxa->EstData == TRUE)
		{
			for(SIndex=0;SIndex<Trees->NoSites;SIndex++)
			{
				if(Taxa->EstDataP[SIndex] == TRUE)
					Taxa->ConData[SIndex] = Rates->EstData[RIndex++];
			}

			if(Taxa->EstDepData == TRUE)
				Taxa->Dependant = Rates->EstData[RIndex++];
		}
	}
}

void	PrintMathmatCode(void)
{
/*	printf("Part3 = -0.5*ZA.Inverse[BlockMatrix[Outer[Times, Sigma, (V^Delta)]]].ZA\n");
	printf("Part2 = Log[Det[BlockMatrix[Outer[Times, (Sigma), (V^Delta)]]]^-0.5]\n");
	printf("Part1 = 1.8378770664093453 ((-(NoTaxa*NoOfSites))/2)\n");
	printf("Lh = Part1 + Part2 + Part3\n");
	printf("Lh = (Log[2 Pi] ((-(NoTaxa*NoOfSites))/2)) + (Log[Det[BlockMatrix[Outer[Times, Sigma, (V^Delta)]]]^-0.5]) + (-0.5*ZA.Inverse[BlockMatrix[Outer[Times, Sigma, (V^Delta)]]].ZA)\n");
*/
//return;
	printf("Part3 = -0.5*ZA.Inverse[ArrayFlatten[Outer[Times, Sigma, V]]].ZA;\n");
//	printf("Part2 = Log[Det[ArrayFlatten[Outer[Times, (Sigma), V]]]^-0.5];\n");
	printf("Part2 = ((NoOfSites * Log[Det[V]]) + (NoTaxa * Log[Det[Sigma]]) * -0.5;\n");
	printf("Part1 = 1.8378770664093453 ((-(NoTaxa*NoOfSites))/2);\n");
	printf("Lh = Part1 + Part2 + Part3;\n");
	printf("Lh = (Log[2 Pi] ((-(NoTaxa*NoOfSites))/2)) + (Log[Det[ArrayFlatten[Outer[Times, Sigma, V]]]^-0.5]) + (-0.5*ZA.Inverse[ArrayFlatten[Outer[Times, Sigma, V]]].ZA)\n");

}


int	CalcInvV(TREES *Trees, TREE *Tree)
{
#ifdef BTOCL
	return btocl_FindInvV(Trees, Tree);
#endif

#ifdef BTLAPACK
	return btlapack_FindInvV(Trees, Tree);
#endif
	return FALSE;
}

double	LHRandWalk(OPTIONS *Opt, TREES* Trees, RATES* Rates)
{
	double	Val;
	MATRIX	*TMat;
	int		Index, Err;
	int		Len;
	double	Det;
	double	Ret;
	TREE	*Tree;
	CONVAR	*CV;


	Err = FALSE;

	Tree = Trees->Tree[Rates->TreeNo];
	CV = Tree->ConVars;

	if(Rates->UseEstData == TRUE)
		SetEstData(Trees, Rates);

	if(Rates->UseEstData == TRUE)
		CalcZ(Trees, Tree, Opt);

#ifdef MATHMAT

	printf("(* **************** Start Tree %d *********************** *)\n", Rates->TreeNo);
	printf("NoOfSites = %d;\n", Trees->NoOfSites);
	printf("NoTaxa = %d;\n", Trees->NoTaxa);

	if(Tree->ConVars->TrueV != NULL)
		PrintMathematicaTFMatrix(Tree->ConVars->TrueV, "V = ", stdout);
	else
		PrintMathematicaTFMatrix(Tree->ConVars->V, "V = ", stdout);

//	PrintMatrix(Tree->ConVars->V, "V = ", stdout);
#endif



	if(Opt->InvertV	== TRUE)
	{
	//	PrintMatrix(Tree->ConVars->V, "V = ", stdout);

		if(Opt->EstKappa == TRUE)
			MakeKappaV(Trees, Tree, Rates->Kappa);
		else
			CopyMatrix(Tree->ConVars->V, Tree->ConVars->TrueV);


		if(Opt->EstOU == TRUE)
			CalcOU(Trees, Tree, Tree->ConVars->V, Rates->OU);

		if(Opt->FixOU != -1)
			CalcOU(Trees, Tree, Tree->ConVars->V, Opt->FixOU);

//		PrintMatrix(Tree->ConVars->V, "V = ", stdout);


		if(Opt->EstDelta == TRUE)
			CalcDelta(Tree->ConVars->V, Rates->Delta);

		if(Opt->FixDelta != -1)
			CalcDelta(Tree->ConVars->V, Opt->FixDelta);

		if(Opt->EstLambda == TRUE)
			CalcLambda(Tree->ConVars->V, Rates->Lambda);

		if(Opt->FixLambda != -1)
			CalcLambda(Tree->ConVars->V, Opt->FixLambda);

	//	PrintMatrix(Tree->ConVars->V, "V = ", stdout);
	//	PrintMathematicaTFMatrix(Tree->ConVars->V, "V = ", stdout);

		if(CalcInvV(Trees, Tree) == FALSE)
			return ERRLH;
	}

//	PrintMatrix(Tree->ConVars->V, "V", stdout);	exit(0);

	CalcSigma(Opt, Trees, Tree, Rates->Means, Rates->Beta);

//	CV->Sigma->me[0][0]  = CV->Sigma->me[0][0]  / (2.0 * Rates->OU);


#ifdef MATHMAT
	PrintMathematicaMatrix(Tree->ConVars->Sigma, "Sigma = ", stdout);
#endif

	if(Opt->Model == M_CONTINUOUS_REG)
	{
		CV->InvSigma->me[0][0] = 1.0 / CV->Sigma->me[0][0];

		CV->LogDetOfSigma = log(CV->Sigma->me[0][0]);

	}
	else
	{
		TMat = AllocMatrix(Trees->NoSites, Trees->NoSites);
		CopyMatrix(TMat, CV->Sigma);
		Err = InvertMatrixAndDet(TMat->me, Trees->NoSites, CV->TVect1, (int*)CV->TVect2, CV->InvSigma->me, &CV->LogDetOfSigma);
		FreeMatrix(TMat);

		if(Err != FALSE)
		{
			return ERRLH;
			printf("Error inverting matrix, matrix is singular, %s::%d\n", __FILE__, __LINE__);
			exit(0);
		}
	}

	CalcZAlpha(Trees, Tree, Opt->Model);

#ifdef MATHMAT
	PrintMathematicaVect(Tree->ConVars->ZA, Trees->NoOfSites * Trees->NoTaxa, "ZA = ", stdout);
#endif


//	KroneckerProduct(Tree->ConVars->InvSigma, Tree->ConVars->InvV, Tree->ConVars->InvKProd);
//	VectByMatrixMult(Tree->ConVars->ZA, Tree->ConVars->InvKProd, Tree->ConVars->ZATemp);
//	btdebug_enter("kronecker");
#ifdef BTOCL
    // only helps for dim > 2000 or sigma > 1
	btocl_VectByKroneckerMult(Tree);
#else
	VectByKroneckerMult(Tree->ConVars->ZA, Tree->ConVars->InvSigma, Tree->ConVars->InvV,Tree->ConVars->ZATemp);
#endif
//	btdebug_exit("kronecker");

	if(Opt->Model == M_CONTINUOUS_REG)
		Len = Trees->NoTaxa;
	else
		Len = Trees->NoSites * Trees->NoTaxa;

	Val = 0;
	for(Index=0;Index<Len;Index++)
		Val += Tree->ConVars->ZATemp[Index] * Tree->ConVars->ZA[Index];

	Val = -0.5 * Val;

	Det = FindDet(Trees, Tree, Opt);

	Det = -0.5 * Det;

	if(Opt->Model == M_CONTINUOUS_REG)
		Ret = -(double)(Trees->NoTaxa);
	else
		Ret = -(double)(Trees->NoSites * Trees->NoTaxa);

	Ret = Ret / 2.0;

	Ret = Ret * 1.837877066409345483560659472811;

#ifdef MATHMAT

	PrintMathmatCode();

	printf("(* Part1 = %f *)\n", Ret);
	printf("(* Part2 = %f *)\n", Det);
	printf("(* Part3 = %f *)\n", Val);

	printf("(* LH = %f *)\n", Ret + Det + Val);

	printf("(* **************** End Tree %d *********************** *)\n", Rates->TreeNo);
	exit(0);
#endif

	Ret = Ret + Det + Val;

	if(ValidLh(Ret, Opt->ModelType) == FALSE)
		Ret = ERRLH;

	Rates->Lh = Ret;

	return Ret;
}

void	TreeBLToPower(TREES *Trees, TREE *Tree, double Power)
{
	int	Index;
	NODE N;

	for(Index=0;Index<Tree->NoNodes;Index++)
	{
		N = Tree->NodeList[Index];
		if(N != Tree->Root)
			N->Length = pow(N->Length, Power);
	}
}

void	InitMultiVarCV(OPTIONS *Opt, TREES* Trees, int TreeNo)
{
	TREE *Tree;
	CONVAR *CV;
	int Index, Size;
	double *FV, *Mean;

	Tree = Trees->Tree[TreeNo];
	CV = Tree->ConVars;

	Mean = (double*)SMalloc(sizeof(double) * Trees->NoTaxa);

	for(Index=0;Index<Trees->NoTaxa;Index++)
		Mean[Index] = 0.0;

	FV = FlattenMatix(CV->V);

	Size = (Trees->NoTaxa * Trees->NoTaxa) + (Trees->NoTaxa + 1);
	CV->MultiVarNormState = (double*)SMalloc(sizeof(double) * Size);
	CV->MultiVarNormTemp = (double*)SMalloc(sizeof(double) * Size);

	setgmn(Mean, FV, Trees->NoTaxa, CV->MultiVarNormState);

	free(Mean);
	free(FV);
}

void	VTest(TREES* Trees, TREE *Tree)
{
	CONVAR	*CV;
//	return;
	CV = Tree->ConVars;
//	PrintMatrix(CV->V, "V=", stdout);


//	MakeKappaV(Trees, Tree, 0.1);

//	CalcOU(CV->V, 0);
//	CaclKappa(CV->V, 2);
//	CalcDelta(CV->V, 2);
//	CalcLambda(CV->V, .5);

	VToTree(CV->V, Tree);

	SaveTrees("treeout3.trees", Trees);
	exit(0);
}


void	InitContinusTree(OPTIONS *Opt, TREES* Trees, int TreeNo)
{
	int		Index;
	CONVAR	*CV;
	TREE	*Tree;

	Trees->Tree[TreeNo]->ConVars = AllocConVar(Opt, Trees);
	Tree = Trees->Tree[TreeNo];
	CV = Tree->ConVars;

	if(Opt->NodeData == TRUE || Opt->NodeBLData == TRUE)
		SetTreeAsData(Opt, Trees, TreeNo);

	if(Opt->FixKappa != -1)
	{
		TreeBLToPower(Trees, Tree, Opt->FixKappa);
		SetTreeDistToRoot(Tree);
	}

	TreeToV(Trees, Tree, CV->V);

	SetNodesVPos(Trees, Tree);

	if(EstData(Trees) == TRUE)
		InitMultiVarCV(Opt, Trees, TreeNo);

	if(Opt->InvertV == TRUE)
		CopyMatrix(CV->TrueV, CV->V);

	if(Opt->FixOU != -1)
		CalcOU(Trees, Tree, CV->V, Opt->FixOU);

	if(Opt->FixDelta != -1)
		CalcDelta(CV->V, Opt->FixDelta);

	if(Opt->FixLambda != -1)
		CalcLambda(CV->V, Opt->FixLambda);

//	VTest(Trees, Trees->Tree[TreeNo]);

	CalcZ(Trees, Tree, Opt);



	CalcInvV(Trees, Trees->Tree[TreeNo]);

	if(Opt->Model == M_CONTINUOUS_REG)
	{
		for(Index=0;Index<Trees->NoTaxa;Index++)
			CV->DepVect[Index] = Trees->Taxa[Index]->Dependant;
	}

}

void		FreeTempConVars(TEMPCONVAR* TempCon)
{
	free(TempCon->T1);
	free(TempCon->T2);
	FreeMatrix(TempCon->TMat);

	FreeMatrix(TempCon->X);
	FreeMatrix(TempCon->TranX);
	FreeMatrix(TempCon->NX);
	free(TempCon->Y);

	FreeMatrix(TempCon->RVX);
	FreeMatrix(TempCon->XT);
	FreeMatrix(TempCon->TempV1);
	FreeMatrix(TempCon->TempV2);

	free(TempCon);
}

TEMPCONVAR* AllocTempConVars(OPTIONS *Opt, TREES* Trees)
{
	TEMPCONVAR*	Ret;
	int			Size;
	TAXA		*Taxa;
	int			x,y;

	Ret = (TEMPCONVAR*) SMalloc(sizeof(TEMPCONVAR));

	Ret->T1 = (double*)SMalloc(sizeof(double)*Trees->NoTaxa);
	Ret->T2 = (int*)SMalloc(sizeof(int)*Trees->NoTaxa);
	Ret->TMat = AllocMatrix(Trees->NoTaxa, Trees->NoTaxa);

	if(Opt->AlphaZero == FALSE)
	{
		Ret->X		= AllocMatrix(Trees->NoTaxa, Trees->NoSites + 1);
		Ret->TranX	= AllocMatrix(Trees->NoSites + 1, Trees->NoTaxa);
		Ret->NX		= AllocMatrix(Trees->NoSites + 1, Trees->NoSites + 1);
	}
	else
	{
		Ret->X		= AllocMatrix(Trees->NoTaxa, Trees->NoSites);
		Ret->TranX	= AllocMatrix(Trees->NoSites, Trees->NoTaxa);
		Ret->NX		= AllocMatrix(Trees->NoSites, Trees->NoSites);
	}

	Ret->Y	= (double*)malloc(sizeof(double) * Trees->NoTaxa);

	if(Opt->AlphaZero == FALSE)
		Size = Trees->NoSites+1;
	else
		Size = Trees->NoSites;

	Ret->RVX	= AllocMatrix(Trees->NoTaxa, Size);
	Ret->XT		= AllocMatrix(Size, Trees->NoTaxa);
	Ret->TempV1	= AllocMatrix(Size, Trees->NoTaxa);
	Ret->TempV2	= AllocMatrix(Size, Size);

	if(Opt->AlphaZero == FALSE)
	{
		for(x=0;x<Trees->NoTaxa;x++)
		{
			Taxa = Trees->Taxa[x];
			Ret->RVX->me[x][0] = 1;
			for(y=0;y<Trees->NoSites;y++)
				Ret->RVX->me[x][y+1] = Taxa->ConData[y];
		}
	}
	else
	{
		for(x=0;x<Trees->NoTaxa;x++)
		{
			Taxa = Trees->Taxa[x];
			for(y=0;y<Trees->NoSites;y++)
				Ret->RVX->me[x][y] = Taxa->ConData[y];
		}
	}

	Transpose(Ret->RVX, Ret->XT);

	return Ret;
}

void	InitContinus(OPTIONS *Opt, TREES* Trees)
{
	int		TIndex;

	CheckZeroTaxaBL(Trees);
	SetTreesDistToRoot(Trees);

	if(Opt->Model == M_CONTINUOUS_REG)
		RemoveDependantData(Opt, Trees);

	AddRecNodes(Opt, Trees);

	Trees->TempConVars = AllocTempConVars(Opt, Trees);

	InitEstData(Opt, Trees);

	Opt->InvertV = FALSE;
	if(	Opt->EstDelta == TRUE ||
		Opt->EstKappa == TRUE ||
		Opt->EstLambda== TRUE ||
		Opt->EstOU == TRUE ||
		Opt->UseVarRates == TRUE)
		Opt->InvertV = TRUE;

	if(Opt->Analsis == ANALMCMC)
	{
		for(TIndex=0;TIndex<Trees->NoTrees;TIndex++)
			InitContinusTree(Opt, Trees, TIndex);
	}
}
