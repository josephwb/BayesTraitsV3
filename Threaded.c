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
#include <string.h>
#include <time.h>

#ifdef CLIK_P
	#include <cilk/cilk.h>
	#include <cilk/cilk_api.h>
#endif

#include "TypeDef.h"
#include "Threaded.h"



int		GetThreadNo(void)
{
#ifdef OPENMP_THR
	return omp_get_thread_num();
#endif	
	return 0;
}

int		GetMaxThreads(void)
{
#ifdef OPENMP_THR
	return omp_get_num_procs();
#endif
	
	return 1;
}

void	SetNoOfThreads(int No)
{
#ifdef USE_MKL
	mkl_set_num_threads(No);
#endif

#ifdef OPENMP_THR
	omp_set_num_threads(No);
	return; 
#endif	

#ifdef CLIK_P
	char *TStr;

	TStr = (char*)SMalloc(sizeof(char) * 64);
	sprintf(TStr, "%d", No);

	if (0 != __cilkrts_set_param("nworkers",TStr))
	{
		printf("Failed to set worker count\n");
		exit(1);
	}

	free(TStr);
	 __cilkrts_init();
	return;
#endif

	return ;
}

double	GetSeconds(void)
{
#ifndef OPENMP_THR
	return  (double)time(NULL);
#else
	return omp_get_wtime();
#endif
}
