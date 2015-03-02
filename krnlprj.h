/*
	krnlprj.h

	$Id: krnlprj.h 30 2005-06-17 15:30:03Z mjs $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

void *pvInitESSEScatterModel(char *pchKrnlFile, char *pchParStr, float
	fPixSize, int nPix, int nSlices, int iCollpaseFac, int bDebugFlag);

void CalcEffScatterSource(void *pvScatterData, float *pfAct, float *pfMap, 
	int nx, int ny, int nz, float *pfEffSrc);

void *pvDeinitScatterPrj(void *pvScatterData);
