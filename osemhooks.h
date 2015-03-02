/*
	osemhooks.h

	This file declares the public hooks API 

	$Id: osemhooks.h 1 2007-05-23 20:58:17Z binhe $
 
 	Copyright 2002-2011 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#ifndef OSEMHOOKS_H
#define OSEMHOOKS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _OsemData
{
	int		iIterationIdx;
	int		iSubsetIdx;
	int		iAngleIdx;
	float *	pfEst;
	float * pfSct;
	float * pfNorm;
	float * pfUpdate;
} sOsemData, * psOsemData;
	
typedef struct _AdvOptions
{
	int		bFastRotate;
	int		iStartIteration;
} sAdvOptions, * psAdvOptions;

int OsemSetAdvOptions(sAdvOptions * pOpt);

int OsemSetPrevIterHook(void (*pPrevIterHookFunc)(sOsemData * psData, void * ctx), void * ctx);
int OsemSetPostIterHook(void (*pPostIterHookFunc)(sOsemData * psData, void * ctx), void * ctx);

int OsemSetPrevSubsetHook(void (*pPrevSubsetHookFunc)(sOsemData * psData, void * ctx), void * ctx);
int OsemSetPostSubsetHook(void (*pPostSubsetHookFunc)(sOsemData * psData, void * ctx), void * ctx);

int OsemSetPrevAngleHook(void (*pPrevAngleHookFunc)(sOsemData * psData, void * ctx), void * ctx);
int OsemSetPostAngleHook(void (*pPostAngleHookFunc)(sOsemData * psData, void * ctx), void * ctx);


/*IrlSetup.c*/
/* Stores arrays of acquisition times and/or detector numbers for each view.
	The detector number is an index into the array of detector data created
	by calls to pvStoreDetectorData. This allows different detectors to have
	different DRFs.
*/
int iSetupExtraPrjViewData(int iNumPrjViews, float *pfAcqTimes, 
			int *piDetectorNums);

/* Stores array of detector data for use by Osem and Genprj. The array
	ppsDetectors is an array of pointers to the data returned by 
	pvStoreDetectorData. This must be called once for each detector used.
	If only one detector is used then IrlOsem will create the detector
	table itself.
*/
int iSetupDetectors(int iNumDetectors, void **ppsDetectors);

/* Stores data for detector used to compute or read DRF */
void *pvStoreDetectorData(
	float fHoleDiam, float fHoleLen, float fBackToDet, /* collimator params */
	float fIntrinsicFWHM,  /* Intrinsic resolution */
	float fSensitivity, /* Sensitivity. Scales output of projector */
	int bFFTconvolve, /* if true DRF convolution uses FFT */
	int bUseGrfInBck, /* if true then GRF is computed using above parameters
								and used in backprojection. This is faster and better
								when DRF has very long tails
							*/
	char *pchDrfTabFname /*if not NULL, then GRF is computed using above 
									parameters. Note that string is not copied
									so memory pointing to it must not be changed.
									It is freed at the end of IrlOsem or IrlGenprj
								*/
);

#ifdef __cplusplus
}
#endif

#endif /* OSEMHOOKS_H */
