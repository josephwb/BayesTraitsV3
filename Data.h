#ifndef DATA_H
#define DATA_H

#include "TypeDef.h"

void		LoadData(char* FileName, TREES *Trees);
void		CheckDataWithModel(char* FileName, TREES *Trees, MODEL Model);
void		PreProcessDataWithModel(TREES *Trees, MODEL Model);

void		PrintData(TREES* Trees, OPTIONS *Opt);
void		PrintTaxaData(OPTIONS *Opt, TREES* Trees);

ANALSIS		GetAnalsis(TREES *Trees);
void		SquashDep(TREES	*Trees);
void		RemoveConMissingData(TREES* Trees);
void		SetTreeAsData(OPTIONS *Opt, TREES *Trees, int TreeNo);
int			EstData(TREES *Trees);

void		FreeData(OPTIONS *Opt);
void		FreeTaxa(TAXA *Taxa, int NoSites);
void		FindSiteSymbols(TREES *Trees, int SiteNo);
void		AddRecNodes(OPTIONS *Opt, TREES *Trees);
char*		SetDescUnknownStates(char S1, char S2);

void		SetDataRegTC(OPTIONS *Opt);
#endif
