/* 
	 IrlSetup.c

	 $Id: IrlSetup.c 124 2011-09-21 15:23:46Z frey $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include <mip/miputil.h>
#include <mip/printmsg.h>
#include <mip/errdefs.h>

#include "irlprivate.h"
#include "osemhooks.h"
#include "drftabio.h"

struct {
	/* Table containing extra data for each projection view that is not stored
		in the PrjViews structure. Data stored includes the detector number
		and the time per projection view. 
		
		iNumPrjViews is the number of projection views for which there is extra
		data. It should either be 0 or equal to the number of views passed
		to IrlOsem.  If iNumPrjViews=0 then, by convention, it is assumed there
		is no extra data per projection view.
		
		The detector number is the index into the ppsDetectors
		array above, and so each must be less than
		sDetectorTable.  or piDetectorNum == NULL
		then this feature is not used and it is assumed
		that all views used the same detector whose parameters and/or drf
		table file name are specified in the Parameters structure.

		The times per projection view are used to
		multiply the projection of the activity distribution and so can
		be in any units. The activity distribution should be in decays per
		the same unit.

		This data is stored in a global variable that is defined in prjviewsetup.c
	*/
	int iNumPrjViews;
	int *piDetectorNums; /* array of detector numbers. */
	float *pfAcqTimes; /* array of times for each projection view */
} sExtraPrjViewDataTable={0,NULL,NULL};

void vSetupParms(IrlParms_t *psIrlParms, Options_t *psOptions, 
		float *pfPrjImage, float *pfAtnMap, float *pfReconImage,
		float *pfScatterEstimate, char *pchScatKrnlFname, Parms_t *psParms)
{
	int ModelAtn, ModelDrf, ModelSrf;

	psParms->NumPixels=psIrlParms->NumPixels;
	psParms->NumBins=psIrlParms->NumPixels;
	psParms->Circular=TRUE;
	psParms->NumRotPixs=psIrlParms->NumPixels;
	psParms->NumAngles=psIrlParms->NumViews;
	psParms->NumSlices=psIrlParms->NumSlices;
	psParms->SliceStart=0;
	psParms->SliceEnd=psIrlParms->NumSlices - 1;
	psParms->BinWidth = psIrlParms->BinWidth;
	psParms->PixelWidth = psIrlParms->BinWidth;
	psParms->SliceThickness = psIrlParms->BinWidth;
	psParms->fAtnScaleFac = psIrlParms->fAtnScaleFac;
	psParms->iNumSrfIterations=psIrlParms->iNumSrfIterations;
	psParms->fScatEstFac = psIrlParms->fScatEstFac;
	psParms->fHoleLen=psIrlParms->fHoleLen;
	psParms->fHoleDiam=psIrlParms->fHoleDiam;
	psParms->fBackToDet=psIrlParms->fBackToDet;
	psParms->fMaxFracErr=0.05; /*Hard Coded  Value!!*/
	psParms->fIntrinsicFWHM=psIrlParms->fIntrinsicFWHM;
	ModelSrf = pchScatKrnlFname == NULL ? 0 : MODEL_SRF;
	ModelDrf = psOptions->bModelDrf ? MODEL_DRF : 0;
	ModelAtn = pfAtnMap == NULL ? 0 : MODEL_ATN;
	
	psParms->iPrjModel= ModelSrf | ModelAtn | ModelDrf;
	psParms->iBckModel= ModelAtn | ModelDrf;
	psParms->iModel= ModelSrf | ModelAtn | ModelDrf;
	psParms->NumAngPerSubset=psIrlParms->NumAngPerSubset;
	psParms->pchNormImageBase=psIrlParms->pchNormImageBase;
	psParms->NumIterations=psIrlParms->NumIterations;
	psParms->StartIteration=0;
	psParms->SaveInterval=1;
	psParms->pchIterSaveString=NULL;
	psParms->iAxialPadLength=psOptions->iAxialPadLength;
	psParms->iAxialAvgLength=psOptions->iAxialAvgLength;
}

void FreeParms(Parms_t *psParms)
/* frees strings and memory for items stored in Parms structure*/
/* note that items that are NULL haven't been used and this is handled
	by IrlFree*/
{
	IrlFree(psParms->pchIterSaveString);
}

View_t *psSetupViews(IrlParms_t *psIrlParms, Options_t *psOptions,
			PrjView_t *psPrjViews)
/* creates and set of views for each angle. This includes things like
	the sine and cosine of the angle and distance from the collimator face
	to the center of rotation (CFCR), Center of rotation (Center)  for each 
	view and the view structure itself.*/
{
	View_t *psViews;
	int iNumAngles;
	double dAng;
	int iAng;

	vPrintMsg(4,"SetupViews Times=%p\n", sExtraPrjViewDataTable.pfAcqTimes);
	iNumAngles = psIrlParms->NumViews;
	psViews = (View_t *)pvIrlMalloc(sizeof(View_t)*iNumAngles,"SetupViews:Views");
	/* compute angles and sine and cosine */
	for(iAng = 0; iAng < iNumAngles; ++iAng){
		psViews[iAng].iAngle = iAng;
		psViews[iAng].Center = (float)psIrlParms->NumPixels/2.0;
		dAng=psViews[iAng].Angle = psPrjViews[iAng].Angle;
		psViews[iAng].Sin = sin(dAng);
		psViews[iAng].Cos = cos(dAng);
		psViews[iAng].CFCR=psPrjViews[iAng].CFCR;
		psViews[iAng].Right = psPrjViews[iAng].Right;
		psViews[iAng].Left = psPrjViews[iAng].Left;
		psViews[iAng].ViewTime=1.0;
		psViews[iAng].Sensitivity=1.0;
		psViews[iAng].psDrfTab = NULL;
#ifdef DEBUG
		fprintf(stderr,"view=%d, angle=%.2f, cfcr=%.2f, left=%d, right=%d\n",
			iAng, psViews[iAng].Angle, psViews[iAng].CFCR, psViews[iAng].Left,
			psViews[iAng].Right);
#endif
	}
	if (sExtraPrjViewDataTable.pfAcqTimes != NULL){
		if (sExtraPrjViewDataTable.iNumPrjViews != iNumAngles)
			vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, 
				"SetupViews", 
			"Number of Prj Views passed to SetupExtrPrjViewData is incorrect");
		for(iAng = 0; iAng < iNumAngles; ++iAng)
			psViews[iAng].ViewTime=sExtraPrjViewDataTable.pfAcqTimes[iAng];
	}
	return psViews;
}


Srf_t *psSrfSetup(Parms_t *psParms, int iSrfCollapseFac, char *pchScatKrnlFname)
/* sets up for srf modeling. Note that this version is for ESSE modeling only
	and for a fixed collapse factor, though scatter modeling can be turned 
	off during the reconstruction using scat_iterations parameter.
*/
{
	Srf_t *psSrf;

	vPrintMsg(4,"SrfSetup\n");
	psSrf=(Srf_t *)pvIrlMalloc(sizeof(Srf_t),"SrfSetup:Srf");
	
	if (iSrfCollapseFac < 1 || iSrfCollapseFac > 4)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, 
			"SrfSetup", 
			"srf_collapse_fac = %f: must be 1 (no collapse) to 4",
			iSrfCollapseFac);
	/* the parm file and parm line are required parameters*/

	/* The ESSE parameters are xy pad factor, z pad factor, number of terms,
		and dimension of the model (must be 3).
		The default is to pad by 2 in both directions and to use 2 terms in
		the expansion of the exponential. */
	psSrf->pvScatterData = pvInitESSEScatterModel( pchScatKrnlFname,
										"2 2 2 3", /*ESSE parms: Hard Coded */
										psParms->PixelWidth,
										psParms->NumRotPixs,
										psParms->NumSlices,
										iSrfCollapseFac, FALSE);
	return(psSrf);
}

Srf_t *psFreeSrf(Parms_t *psParms, Srf_t *psSrf)
{
	if (psSrf != NULL){
		pvDeinitScatterPrj(psSrf->pvScatterData);
		IrlFree(psSrf);
	}
	return NULL;
}

/* relative difference in pixel size in DRF table and reconstruction */
#define PIXSIZE_TOLERANCE 0.01 

/* relative difference allowed in min and max distances*/
#define DIST_TOLERANCE 0.0001

/* compute the relative difference of 2 fp numbers. They are normalized by
	the average of the absolute value and if this is zero then the absolute
	difference is computed */
#define RELDIF(a,b) (2.0*fabs((a)-(b))/(((fabs(a)+fabs(b))==0.0)? 1 : (fabs(a)+fabs(b))))

/* compare floating point numbers, in terms of the relative difference tol
	difference*/
#define CMPRFLOAT(a,b,tol) (RELDIF((a),(b)) <= tol)

static void vCheckDrfTab(DrfTab_t sDrfTab, float fPixWid,
				float fMinDist, float fMaxDist, char *pchFname)
/* checks the DrfTab read from the file to see if the parameters agree
	are suitable. This includes the pixel sizes being the same (required),
	the min and max distances are included in the DrfTab (warning only),
	and the Drf spacing (distance between planes in the Drf Table) and the
	pixel size
*/
{
	float fDrfTabMaxDist;

	/* the relative difference in the pixel size in the Drf tab and the 
		reconstruction must be within PIXSIZE_TOLERANCE of each other */
	if (!CMPRFLOAT(sDrfTab.fPixSize, fPixWid, PIXSIZE_TOLERANCE)) 
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "vCheckDrfTab",
	"The pixel size in Drf Table %s is %.4f, and not %.4f as needed for recon.",
			pchFname, sDrfTab.fPixSize, fPixWid);

	fDrfTabMaxDist = sDrfTab.fStartDistance + 
								(sDrfTab.iNumDistances-1)*sDrfTab.fDrfSpacing;
	/* check that the minimum (fMinDist) and maximum (fMaxDist) for the
		reconstruction are within the the range for the DRF table. To be
		o.k., the mindist must be within the specified tolerance or greater
		than of the Drf table starting distance. Also, the maxdist must be
		within the specified tolerance or less than the drf table max distance*/
	if (! (((CMPRFLOAT(fMinDist,sDrfTab.fStartDistance, DIST_TOLERANCE) ||
			  fMinDist > sDrfTab.fStartDistance) && 
			(CMPRFLOAT(fMaxDist,fDrfTabMaxDist,DIST_TOLERANCE) ||
				fMaxDist < fDrfTabMaxDist))
			|| (fMaxDist > fDrfTabMaxDist)))
		vErrorHandler(ECLASS_WARN, ETYPE_ILLEGAL_VALUE, "vCheckDrfTab",
			"The range of distances in %s is (%.2f to %.2f), \nreconstruction needs (%.2f to %.2f)\nApproximations will be used",
				pchFname, sDrfTab.fStartDistance, fDrfTabMaxDist,
				fMinDist, fMaxDist);
	if (!CMPRFLOAT(sDrfTab.fDrfSpacing,fPixWid, 
			PIXSIZE_TOLERANCE))
		vErrorHandler(ECLASS_WARN, ETYPE_ILLEGAL_VALUE, "vCheckDrfTab",
			"The distance between Drfs in %s is %.4f, not %.4f. NN interp will be used",
			pchFname, sDrfTab.fDrfSpacing, fPixWid);
}

static void vComputeDrfs(DrfTab_t *psDrfTab, 
	float fPixelWidth, float fBinWidth, int iReconMatSize, 
	float fMinDist, float fMaxDist, int iIncKrnlOffsets)
/* computes the DRF tables from the detector parameters for each detector
	stored in the Detectors*/
{
	int iNumDists;
	int bDoIncBlur;
	int iNumDrfOffs = (int)(iReconMatSize*0.5);/* maximum half-size of Drf*/
	float fGap, fHoleDiam, fColThick, fIntrinsicFWHM, fMaxFracErr;

	/* get collimator parameters*/
	fGap = psDrfTab->fGap;
	fColThick = psDrfTab->fColThick;
	fHoleDiam = psDrfTab->fHoleDiam;
	fIntrinsicFWHM = psDrfTab->fIntrinsicFWHM ;

	fMaxFracErr = psDrfTab->fMaxFracErr;

	vPrintMsg(6,"ComputeDrfs: gap=%.2f, collthick=%.2f, intrinsic=%.4f, holediam=%.4f\n",
		fGap, fColThick, fIntrinsicFWHM, fHoleDiam);
	if (iIncKrnlOffsets < 0 || iIncKrnlOffsets > 4)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "vComputeDrfs",
			"Number of offsets in inc kernel must be >=0 and <=4, not %d",
			iIncKrnlOffsets);
	bDoIncBlur = (iIncKrnlOffsets != 0);
	iNumDists = ((fMaxDist - fMinDist)/fPixelWidth + 1 + 0.5); 

	vPrintMsg(6, "minDist = %f, maxDist = %f\n",fMinDist, fMaxDist); 
	vPrintMsg(4, "generating detector response function ...\n");
	vPrintMsg(6, "iNumDists = %d, iNumOffs=%d\n",iNumDists,iNumDrfOffs); 

	psDrfTab->psDrfs = (Drf_t *)pvIrlMalloc(iNumDists*sizeof(Drf_t),
				"ComputeDrfs:psDrfs");
	psDrfTab->iNumDistances = iNumDists;
	psDrfTab->fStartDistance=fMinDist;
	psDrfTab->fDrfSpacing=fPixelWidth;
	psDrfTab->fPixSize=fPixelWidth;
	psDrfTab->bIncKrnl = bDoIncBlur;
	gendrf(PARALLEL,0.0, fColThick,
			fHoleDiam, fGap, fIntrinsicFWHM, fMaxFracErr,
			fBinWidth, fPixelWidth, fMinDist, iNumDists, iNumDrfOffs,
			psDrfTab->psDrfs);
	if (bDoIncBlur){
#ifdef INCBLUR
#ifdef DEBUG
		int i;
		pDrf=psDrfTab->psDrfs; /* points to Drf for iDist=0*/
		psDrfTab->bIncKrnl = TRUE;
		vPrintMsg(9,"first psf sum=%.4g @ %x\n", 
			sum_float(pDrf->Psf, (2*(pDrf->MaxOffsX-1)* (2*(pDrf->MaxOffsY-1)),
				pDrf.Psf);
#endif
	/* the below needs to be fixed if INCBLUR is used since psParms is no
		longer avaiable in this function.
	*/
	vCalcIncKrnls(iNumDists, psParms->NumBins, psParms->NumSlices, 
			iIncKrnlOffsets, iIncKrnlOffsets,
			psDrfTab->psDrfs);
#ifdef DEBUG
	writeimage("firstpsf.im", 
			(2*(pDrf->MaxOffsX-1),
			(2*(pDrf->MaxOffsY-1)),
			pDrf.Psf);
	vPrintMsg(9,"first psf sum after CalcIncKrnls=%.4g @ %x\n", 
		sum_float(pDrf.Psf,(2*(pDrf->MaxOffsX-1)* (2*(pDrf->MaxOffsY-1)),
		pDrf.Psf);
	for(i=0; i<iNumDists; ++i){
		int iOffsX=pDrfTabs[iDet].psDrfs[i].MaxOffsX;
		int iOffsY=pDrfTabs[iDet].psDrfs[i].MaxOffsY;
		float *Psf=pDrfTabs[iDet].psDrfs[i].Psf;

		vPrintMsg(9,
			"dist=%d, iMaxOffs=%d x %d, sum=%.4g b[-1,0,1]=%.4g %.4g %.4g\n",
			i, iOffsX, iOffsY,
			sum_float(Psf, (2*iOffsX-1)*(2*iOffsY-1)),
			Psf[iOffsX-2 + (2*iOffsX-1)*(iOffsY-1)],
			Psf[iOffsX-1 + (2*iOffsX-1)*(iOffsY-1)],
			Psf[iOffsX + (2*iOffsX-1)*(iOffsY-1)]);
	}
#endif
#else
	vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "vComputeDrfs",
		"compiled without incremental blurring, use inc_krnl_offsets=0");
#endif
	}else{
	/* no incremental blurring */
		psDrfTab->bIncKrnl=FALSE;
	}
	vPrintMsg(4, "done drf generation.\n");
	return;  
}


struct {
	/* Table of data for each detector. 
		The array ppsDetectors is an array of pointers to structs
		that contains the information (filename or collimator parameters) needed
		to read or compute the DRF Table for each detector. The detectors are 
		mapped to projection views by sDetectorNumTable.

		This data is stored in a global variable that is defined in prjviewsetup.c
	*/
	int iNumDetectors;
	Detector_t **ppsDetectors;
} sDetectorDataTable= {0,NULL};

int iSetupExtraPrjViewData(
	int iNumPrjViews,
	float *pfAcqTimes,
	int *piDetectorNums
)
{
	int i;
	int iError=0;

	vPrintMsg(4,
		"SetupExtraPrjViewData iNumPrjViews=%d AcqTimes=%c DetectorNums=%c\n",
		iNumPrjViews, pfAcqTimes == NULL ? 'n':'y',
		piDetectorNums == NULL ? 'n':'y');
	/* store the data */
	sExtraPrjViewDataTable.iNumPrjViews=iNumPrjViews;
	/* check for NULLs and only set if the pointers are non-null. This allows
		this function to be called independently to set AcqTimes and DetectorNums
	*/
	if (pfAcqTimes != NULL){
		sExtraPrjViewDataTable.pfAcqTimes=pfAcqTimes;
	}
	if (piDetectorNums != NULL){
		sExtraPrjViewDataTable.piDetectorNums=piDetectorNums;
	}
	/* Check for illegal values */
	
	if (iNumPrjViews < 0 )
		iError=1;

	if (iNumPrjViews > 0 && (pfAcqTimes == NULL && piDetectorNums == NULL))
		iError=1;

	if (iNumPrjViews == 0 && (pfAcqTimes != NULL || piDetectorNums != NULL))
		iError = 1;
	
	if (pfAcqTimes != NULL){
		for(i=0; i<iNumPrjViews;++i){
			if (pfAcqTimes[i] < 0.0){
				iError |= 2;
			}
		}
	}
	if (piDetectorNums != NULL){
		for(i=0; i<iNumPrjViews;++i){
			if (piDetectorNums[i] < 0) {
				iError |= 4;
			}
		}
	}
	if (iError != 0)
		vPrintMsg(6,"SetupExtraPrjViewData, iError=%d\n",iError);
	return iError;
}

void vFreeExtraPrjViewData()
{
	IrlFree(sExtraPrjViewDataTable.pfAcqTimes);
	sExtraPrjViewDataTable.pfAcqTimes = NULL;

	IrlFree(sExtraPrjViewDataTable.piDetectorNums);
	sExtraPrjViewDataTable.piDetectorNums = NULL;

	sExtraPrjViewDataTable.iNumPrjViews=0;
}

void *pvStoreDetectorData(
	float fHoleDiam,
	float fHoleLen,
	float fBackToDet,
	float fIntrinsicFWHM,
	float fSensitivity,
	int bFFTconvolve,
	int bUseGrfInBck,
	char *pchDrfTabFname
)
{
	Detector_t *psDetector;

	psDetector=pvIrlMalloc(sizeof(Detector_t),"StoreDetectorData: sDetector");
	psDetector->fHoleDiam=fHoleDiam;
	psDetector->fHoleLen=fHoleLen;
	psDetector->fBackToDet=fBackToDet;
	psDetector->fIntrinsicFWHM=fIntrinsicFWHM;
	psDetector->fSensitivity=fSensitivity;
	psDetector->pchDrfTabFname=pchDrfTabFname;
	psDetector->bFFTconvolve=bFFTconvolve;
	psDetector->bUseGrfInBck=bUseGrfInBck;
	psDetector->psDrfTab=NULL;
	return (void *)psDetector;
}

int iSetupDetectors(
	int iNumDetectors,
	void **ppsDetectors
)
{
	int iError=0;

	vPrintMsg(4,"iSetupDetectors: NumDets=%d, ppsDetectors=%p\n",
		iNumDetectors,ppsDetectors);
	printf("iSetupDetectors: NumDets=%d, ppsDetectors=%p\n",
		iNumDetectors,ppsDetectors);
	sDetectorDataTable.iNumDetectors=iNumDetectors;
	sDetectorDataTable.ppsDetectors=(Detector_t **)ppsDetectors;

	if (iNumDetectors < 0)
		iError=1;
	if (iNumDetectors > 0 && ppsDetectors == NULL)
		iError=2;
	if (iNumDetectors == 0 && ppsDetectors != NULL)
		iError = 4;
	return iError;
}

void vFreeDetectors()
{
	int iDet;
	Detector_t **ppsDetectors;
		
	/* free data items in each detector */

	ppsDetectors=sDetectorDataTable.ppsDetectors;
	for(iDet=0; iDet < sDetectorDataTable.iNumDetectors; ++iDet){
		IrlFree(ppsDetectors[iDet]->pchDrfTabFname);
		ppsDetectors[iDet]->psDrfTab=psFreeDrf(ppsDetectors[iDet]->psDrfTab);
		IrlFree(ppsDetectors[iDet]);
		ppsDetectors[iDet]=NULL;
	}
	IrlFree(ppsDetectors);
	sDetectorDataTable.ppsDetectors=NULL;
	sDetectorDataTable.iNumDetectors=0;
}

static void vFindMaxMinCFCR(int iNumViews, View_t *psViews, float *pfMinCFCR, 
		float *pfMaxCFCR);
static void vFindMaxMinCFCR(int iNumViews, View_t *psViews, float *pfMinCFCR, 
		float *pfMaxCFCR)
{
	int iAng;
	float fMinCFCR, fMaxCFCR;

	fMinCFCR=0.0;
	fMaxCFCR=0.0;
	if (iNumViews > 0) fMinCFCR=fMaxCFCR=psViews[0].CFCR;
	for(iAng=0; iAng < iNumViews; ++iAng){
		fMinCFCR= fMinCFCR > psViews[iAng].CFCR ?  psViews[iAng].CFCR : fMinCFCR;
		fMaxCFCR= fMaxCFCR < psViews[iAng].CFCR ?  psViews[iAng].CFCR : fMaxCFCR;
	}
	*pfMinCFCR=fMinCFCR;
	*pfMaxCFCR=fMaxCFCR;
	vPrintMsg(6,"maxCFCR=%.1f, minCFCR=%.1f\n",fMaxCFCR, fMinCFCR);
}
		

DrfTab_t *psDrfSetup(
	float fMinCFCR,
	float fMaxCFCR,
	float fPixelWidth,
	float fBinWidth,
	int iReconMatSize,
	float fMaxFracErr,
	Detector_t *psDetectorParms
	);

DrfTab_t *psDrfSetup( float fMinCFCR, float fMaxCFCR, 
	float fPixelWidth, float fBinWidth, int iReconMatSize, float fMaxFracErr,
	Detector_t *psDetectorParms)
/* this reads or computes the Drf tables for each detector (if 
	necessary). All data needed by the projector and backprojector to 
	model the drf is returned as a pointer to a DrfTab_t structure.
	If no drf modeling is requested then this routine returns NULL.
	If psDetectorParms->bFFTconvolve is true, then the convolution is done in the
	Fourier Domain. If psDetector-Parms->bUseGrfInBck is set, then a flag is set
	to use the GRF in the backprojector and the GRF is computed using the 
	collimator parameters and stored in the psBckDrfTab element of the DRF table.
*/
{
	DrfTab_t *pDrfTab;
	char *pchColName=NULL, *pchDrfTitle=NULL; 
	int bDrfFromFile;
	float fMinDist, fMaxDist;
	
	bDrfFromFile = psDetectorParms->pchDrfTabFname != NULL;
	vPrintMsg(7,"DrfFromFile=%s\n",bDrfFromFile ? "true":"false");
	pDrfTab=(DrfTab_t *)pvIrlMalloc(sizeof(DrfTab_t),
			"DrfSetup:DrfTab");
	pDrfTab->psBckDrfTab = NULL;
		
	fMinDist = fMinCFCR - iReconMatSize/2*fPixelWidth + fPixelWidth/2;
	if (fMinDist < 0) fMinDist=fPixelWidth/2;
	fMaxDist = fMaxCFCR + iReconMatSize/2*fPixelWidth - fPixelWidth/2;
	/* set the collimator parameters */
	pDrfTab->fGap=psDetectorParms->fBackToDet;
	pDrfTab->fHoleDiam=psDetectorParms->fHoleDiam;
	pDrfTab->fColThick=psDetectorParms->fHoleLen;
	pDrfTab->fIntrinsicFWHM=psDetectorParms->fIntrinsicFWHM;
	pDrfTab->fMaxFracErr=fMaxFracErr;

	if (bDrfFromFile){
		vPrintMsg(6,"Reading Drf table from file \n");
		/* for this case the drf table files are gotten directly from the
			parameter file drf_tab_file*/
		*pDrfTab = sReadDrfTab(psDetectorParms-> pchDrfTabFname, 
						&pchDrfTitle, &pchColName);
		vCheckDrfTab(*pDrfTab, fPixelWidth, fMinDist, fMaxDist, 
			psDetectorParms->pchDrfTabFname);
		vPrintMsg(6," colname=%s\n   %s\n", pchColName, pchDrfTitle);
		if (psDetectorParms->bUseGrfInBck){
			pDrfTab->psBckDrfTab=(DrfTab_t *)pvIrlMalloc(sizeof(DrfTab_t),
				"DrfSetup:DrfTab");
			(pDrfTab->psBckDrfTab)->fGap=psDetectorParms->fBackToDet;
			(pDrfTab->psBckDrfTab)->fHoleDiam=psDetectorParms->fHoleDiam;
			(pDrfTab->psBckDrfTab)->fColThick=psDetectorParms->fHoleLen;
			(pDrfTab->psBckDrfTab)->fIntrinsicFWHM=psDetectorParms->fIntrinsicFWHM;
			(pDrfTab->psBckDrfTab)->fMaxFracErr=fMaxFracErr;
			vPrintMsg(8,"Computing GRF for Bck\n");
			vComputeDrfs(pDrfTab->psBckDrfTab, fPixelWidth, fBinWidth, 
					iReconMatSize, fMinDist, fMaxDist, 0);
			(pDrfTab->psBckDrfTab)->psBckDrfTab = NULL;
			(pDrfTab->psBckDrfTab)->bFFTconvolve= FALSE;
		}else{
			pDrfTab->psBckDrfTab = NULL; /* use the DRF in both Prj and Bck */
		}
	}else{/* compute the drf for each detector */
		vComputeDrfs(pDrfTab, fPixelWidth, fBinWidth, iReconMatSize, 
			fMinDist, fMaxDist, 0);
	}
	pDrfTab->bFFTconvolve = psDetectorParms->bFFTconvolve;
	if (pDrfTab->bFFTconvolve)
		vPrintMsg(6,"convolve DRF using fft\n");

	/*free the titles returned by vReadDrfParms since they were strdup-ed*/
	if (pchColName != NULL)IrlFree(pchColName);
	if (pchDrfTitle != NULL)IrlFree(pchDrfTitle);
#ifdef DEBUG
	fprintf(stderr,"psDrfTab=%d, psBckDrfTab=%d\n",pDrfTab, pDrfTab->psBckDrfTab);
#endif
	return pDrfTab;
}

void vSetupDrfTables(Parms_t *psParms, Options_t *psOptions,
		View_t *psViews, char *pchDrfTabFname)
/* Computes or reads DRFs for all detectors and then stores pointers
	to the appropriate drf tables in the structure for each view. 
	If multiple detectors (i.e., detectors with different properties such
	as colliamtors) are used, then this is setup by calling vStoreDetectorData
	for each detector and then passing the array of detector data pointers to
	iSetupDetectors. If only a single detector is used then the DRF
	parameters are stored in the parms array or the DrfTabFname is passed in. 
	(Note: If the DrfTabFname points to Null, then a computed grf is used.)
	This routine creates a single detector entry and calls
	iSetupDetectors with it.  DrfTables are either read or computed,
	depending on the options specified, by the routine psDrfSetup.
*/ {
	int iAng, iDet; 
	int iNumDetectors; 
	int *piDetectorNums; 
	int ierr;
	Detector_t **ppsDetectors, *psDetectors; 
	float fMinCFCR, fMaxCFCR;

	if (sDetectorDataTable.iNumDetectors==0){
		/* if no detectors have been setup using iSetupDetectors, then
			there must be only 1 detector. Set it up from the parameters
			array to store pointer to detectors. Note that there is no
			sensitivity parameter in psOptions, so it is assumed to be
			1. */
		ppsDetectors=pvIrlMalloc(sizeof(Detector_t *),
			"vSetupDrfTables:ppsDetectors");
		psDetectors = ppsDetectors[0]= (Detector_t *)pvStoreDetectorData(
				psParms->fHoleDiam,psParms->fHoleLen,psParms->fBackToDet,
				psParms->fIntrinsicFWHM, 1.0,
				psOptions->bFFTConvolve, psOptions->bUseGrfInBck, pchDrfTabFname);
		if ((ierr=iSetupDetectors(1,(void *) ppsDetectors)) != 0)
			vErrorHandler(ECLASS_FATAL, ETYPE_UNKNOWN, "vSetupDrfTables",
				"Error seting up single detector. This shouldn't happen. err=%d",
				ierr);
	}

	vFindMaxMinCFCR(psParms->NumAngles, psViews, &fMinCFCR, &fMaxCFCR);
	iNumDetectors=sDetectorDataTable.iNumDetectors;
	ppsDetectors=sDetectorDataTable.ppsDetectors;
	for(iDet=0; iDet < sDetectorDataTable.iNumDetectors; ++iDet){
		/* Compute the DRF for each detector */
		if (psParms->iPrjModel & MODEL_DRF){
			ppsDetectors[iDet]->psDrfTab = 
				psDrfSetup(fMinCFCR, fMaxCFCR, psParms->PixelWidth,
						psParms->BinWidth, psParms->NumRotPixs,
	  					psParms->fMaxFracErr, ppsDetectors[iDet]);
		}else{
			ppsDetectors[iDet]->psDrfTab = NULL;
		}
	}
	
	piDetectorNums=sExtraPrjViewDataTable.piDetectorNums;
	for(iAng=0; iAng < psParms->NumAngles; ++iAng){
		if (piDetectorNums == NULL){
			psViews[iAng].psDrfTab=ppsDetectors[0]->psDrfTab;
			psViews[iAng].Sensitivity=1.0;
		}else{
			iDet=piDetectorNums[iAng];
			if (iDet < 0 || iDet >= iNumDetectors)
				vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "vSetupDrfTables",
					"Illegal detector specified for view %d (%d)", iAng, iDet);
			psViews[iAng].psDrfTab=ppsDetectors[iDet]->psDrfTab;
			psViews[iAng].Sensitivity=ppsDetectors[iDet]->fSensitivity;
		}
	}
}

