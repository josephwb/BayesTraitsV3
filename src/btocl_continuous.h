#ifndef BTOCL_CONTINUOUS_H
#define BTOCL_CONTINUOUS_H

#ifdef BTOCL
#include "btocl_runtime.h"
#include "btocl_kernels_bayestraits.h"
#include "btocl_lin.h"

int		btocl_FindInvV(TREES *Trees, TREE* Tree);
void    btocl_VectByKroneckerMult(TREE* Tree);
void	btocl_AllocConVar(CONVAR* ConVar, TREES *Trees);
void	btocl_FreeConVar(CONVAR* ConVar);
#endif

#endif
