#ifndef MINFIT
#define MINFIT

#include "Praxis.h"

/*
#define EPSILON 1.0e-8
#define SQREPSILON 1.0e-16
#define N 20
*/
void minfit(int n,double  eps, double tol, double ab[PNSIZE][PNSIZE], double q[PNSIZE]);

#endif
