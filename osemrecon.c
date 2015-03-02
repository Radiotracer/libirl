/*
	osemrecon.c

	$Id: osemrecon.c 132 2014-07-04 21:51:37Z frey $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#ifndef WIN32
#include <unistd.h>
#endif
#ifdef WIN32
#define unlink _unlink
#endif

#include <mip/miputil.h>
#include <mip/errdefs.h>
#include <mip/printmsg.h>

#define NORM_EXTENSION ".dat"

#ifdef DEBUG 
#include <mip/imgio.h>
#endif

#include "irlprivate.h"
#include "osemhooks.h"
static void (*g_pPrevIterHookFunc)(sOsemData * psData, void * ctx) = NULL;
static void (*g_pPostIterHookFunc)(sOsemData * psData, void * ctx) = NULL;
static void (*g_pPrevAngleHookFunc)(sOsemData * psData, void * ctx) = NULL;
static void (*g_pPostAngleHookFunc)(sOsemData * psData, void * ctx) = NULL;
static void (*g_pPrevSubsetHookFunc)(sOsemData * psData, void * ctx) = NULL;
static void (*g_pPostSubsetHookFunc)(sOsemData * psData, void * ctx) = NULL;
static void * g_pCtx;
static sOsemData g_sOsemHookData;
static sAdvOptions g_sAdvOpts = {1,0}; //fastrotate , start iteration #

int OsemSetAdvOptions(sAdvOptions * pOpt)
{
	memcpy(&g_sAdvOpts, pOpt, sizeof(sAdvOptions));
	return 1;
}

int OsemSetPrevIterHook(void (*pFunc)(sOsemData * psData, void * ctx), void * ctx)
{
	g_pPrevIterHookFunc = pFunc;
	g_pCtx  = ctx;
	return 1;
}

int OsemSetPostIterHook(void (*pFunc)(sOsemData * psData, void * ctx), void * ctx)
{
	g_pPostIterHookFunc = pFunc;
	g_pCtx  = ctx;
	return 1;
}

int OsemSetPrevAngleHook(void (*pFunc)(sOsemData * psData, void * ctx), void * ctx)
{
	g_pPrevAngleHookFunc = pFunc;
	g_pCtx  = ctx;
	return 1;
}

int OsemSetPostAngleHook(void (*pFunc)(sOsemData * psData, void * ctx), void * ctx)
{
	g_pPostAngleHookFunc = pFunc;
	g_pCtx  = ctx;
	return 1;
}

int OsemSetPrevSubsetHook(void (*pFunc)(sOsemData * psData, void * ctx), void * ctx)
{
	g_pPrevSubsetHookFunc = pFunc;
	g_pCtx  = ctx;
	return 1;
}

int OsemSetPostSubsetHook(void (*pFunc)(sOsemData * psData, void * ctx), void * ctx)
{
	g_pPostSubsetHookFunc = pFunc;
	g_pCtx  = ctx;
	return 1;
}

static void SetMeasuredPrjBinsToOne(float *pfOutPrj,
						float *pfMeasPrj,
						Parms_t *psParms,
						View_t *View)
	{
	/* modified 3/25/98 as per Dan Kadrmas to treat measured bins with
	negative value as nonmeasured*/
#define mod_ptr(iSlice, iBin) ((iSlice)*nRotPixs+(iBin))
	/* pixel width in units of bin width */ 
	float PWid = psParms->PixelWidth / psParms->BinWidth;
	float PixelLeft, PixelRight;
	int	BinLeft,	BinRight;
	int	iSlice, iPixel;
	int 	nRotPixs = psParms->NumRotPixs;
	int 	nBins = psParms->NumBins;
	int 	nSlices = psParms->NumSlices;

	set_float(pfOutPrj, nSlices*nRotPixs, 0.0);

	PixelLeft = View->Left;

	vPrintMsg(8,"SetMeasPrjBinsToOne: npix=%d,nslice=%d,nbin=%d,start=%.3f\n",
			nRotPixs, nSlices, nBins,sum_float(pfOutPrj,nSlices*nRotPixs));
	vPrintMsg(8,"    view->Left = %.2f\n",View->Left);
	for (iPixel=0; iPixel<nRotPixs; iPixel++) {
		PixelRight = PixelLeft + PWid;
		BinLeft  = floor((double)PixelLeft);
		BinRight = floor((double)PixelRight);
		if (BinLeft > View->Right) break;

/* Modified 4/10/98 by Dan Kadrmas
		Changed limits to include all (ModPrj)Bins that are at least
		partly measured by UsrPrj.  This helps to reduce truncation artifacts
		with convergent beam collimators.  Caution:  there may be a tendency
		to cause artifacts when using Bins that are barely in the FOV (ie, if
		only 1% of the bin was measured, for example). */
		if ((BinRight>0) && (BinLeft<nBins)){
			for (iSlice=0; iSlice<nSlices; iSlice++) {
				if (pfMeasPrj[mod_ptr(iSlice,iPixel)] >= 0.0)
					pfOutPrj[mod_ptr(iSlice,iPixel)] = 1.0;
			 }
		}
		PixelLeft = PixelRight;
	}
	vPrintMsg(8,"     end=%.3f\n",sum_float(pfOutPrj,nSlices*nRotPixs));
#undef mod_ptr
}

static void SetNonMeasuredBinsToZero(float *pfPrj,
							float *pfMeasPrj,
							Parms_t *psParms,
							View_t *View)
{
/* modified 3/25 as per Dan Kadrmas to treat measured bins with
	negative value as nonmeasured*/
#define mod_ptr(iSlice, iBin) ((iSlice)*nRotPixs+(iBin))
	/* pixel width in units of bin width */ 
	float PWid = psParms->PixelWidth / psParms->BinWidth;

	float PixelLeft, PixelRight;
	int	BinLeft,	BinRight;
	int	iSlice, iPixel;
	int 	nRotPixs = psParms->NumRotPixs;
	int 	nBins = psParms->NumBins;
	int 	nSlices = psParms->NumSlices;

	PixelLeft = View->Left;

	for (iPixel=0; iPixel<nRotPixs; iPixel++) {
		PixelRight = PixelLeft + PWid;
		BinLeft  = floor((double)PixelLeft);
		BinRight = floor((double)PixelRight);

/* Modified 4/10/98 by Dan Kadrmas
		Changed limits to include all (ModPrj)Bins that are at least
		partly measured by UsrPrj.  This helps to reduce truncation artifacts
		with convergent beam collimators.  Caution:  there may be a tendency
		to cause artifacts when using Bins that are barely in the FOV (ie, if
		only 1% of the bin was measured, for example). */
		if ((BinRight<=0) || (BinLeft>=nBins)){
			for (iSlice=0; iSlice<nSlices; iSlice++) {
				pfPrj[mod_ptr(iSlice,iPixel)] = 0.0;
			 }
		}else{
			for(iSlice=0; iSlice<nSlices; iSlice++){
				if (pfMeasPrj[mod_ptr(iSlice,iPixel)] < 0.0)
					pfPrj[mod_ptr(iSlice,iPixel)] = 0.0;
			}
		}
		PixelLeft = PixelRight;
	}
#undef mod_ptr
}

static void writedat(char *fname, int iNbytes, void *pBuf)
{
	FILE *fp;
	int sz;
	
	fp=fopen(fname,"wb");
	if (fp==NULL){
		vErrorHandler(ECLASS_FATAL, ETYPE_IO, "writedat",
			"IO error opening to %s for writing", fname);
	}
	sz=fwrite(pBuf, 1, iNbytes, fp);
	if (sz != iNbytes || ferror(fp))
		vErrorHandler(ECLASS_FATAL, ETYPE_IO, "writedat",
			"IO error while writing to %s", fname);
	fclose(fp);
}

static int readdat(char *fname, int iStart, int iNbytes, void *pBuf)
{
	FILE *fp;
	int sz;

	fp=fopen(fname,"rb");
	if (fp==NULL){
		vErrorHandler(ECLASS_FATAL, ETYPE_IO, "readdat",
			"IO error opening to %s for reading", fname);
	}
	if (iNbytes > 0){
		if (fseek(fp, iStart, SEEK_SET) != 0)
			vErrorHandler(ECLASS_FATAL, ETYPE_IO, "readdat",
				"IO error seeking to position %d in %s ", iStart, fname);
	}
	sz=fread(pBuf, 1, iNbytes, fp);
	fclose(fp);
	return sz;
}

static NormData_t *psComputeNormImages(
	char *pchBase, /*Base name for norm images*/
	Parms_t *psParms, /* general parameters*/
	View_t *psViews, /*per-view data for reconstruction*/
	Srf_t *psSrf, /* data for scatter modeling*/
	int iNumSubsets, /* number of subsets*/
	int **ppiSubsetAngles, /* array of pointers to array containing angle 
										increments numbers for each subset */
	float *pfCalcPrj,/* projection image to be flled with ones and backproject*/
	float *pfPrjImage,/*Measured projections. Need these to help figure out
								which pixels were measured*/
	float *pfAtnMap, /* attenuation map if performing atten comp*/
	float *pfNormImage)/* Image used to generate each norm image*/
{
	int iSubset;
	int iSetMember;
	int iAngle;
	int iPix;
	int iNumSlices = psParms->NumSlices;
	int iNumPix = psParms->NumPixels;
	int iNumRotPix = psParms->NumRotPixs;
	int iNumAngles= psParms->NumAngles;
	int iNumAngPerSet = psParms->NumAngPerSubset;
	int iReconImageSize = iNumPix*iNumPix*iNumSlices;
	int iPrjSize = iNumRotPix * iNumSlices;
	int bNormInMemory=FALSE;
	char *pchNormBase=NULL;
	char *pchNormName=NULL;
	float *pfNormImages=NULL;/* used to store norm images if they are stored in
									memory*/
	NormData_t *psNormData;
	
	vPrintMsg(8,
		"prjsize=%d, ReconImageSize=%d, NumRotPix=%d, NumSlices=%d,NumPix=%d\n",
		iPrjSize,iReconImageSize, iNumRotPix, iNumSlices, iNumPix);
	vPrintMsg(4,"ComputeNormImages\n");
	if (pchBase == NULL){
		vPrintMsg(8,"Norm In Memory\n");
		bNormInMemory=TRUE;
		pfNormImages = vector(iReconImageSize*iNumSubsets,
								"ComputeNormImages:NormImages (in memory)");
		set_float(pfNormImages, iReconImageSize*iNumSubsets,0.0);
	}else{
		/* make base file name. File name is format:
			base.norm.subset.extension, where:
				base is passed in as pchBase 
				subset is the subset number
				extension is defined by NORM_EXTENSION
		*/
		/*to be safe, make output string extra 25 characters long: 
			5 for .norm, 
			5 for .extension, 6 for subset number, one for / and 1 for null*/
		pchNormBase = pvIrlMalloc(+strlen(pchBase)+strlen(NORM_EXTENSION)+15,
			"ComputeNormImages:NormBase");
		pchNormName = pvIrlMalloc(strlen(pchBase)+strlen(NORM_EXTENSION)+25,
			"ComputeNormImages:NormName");
		sprintf(pchNormBase,"%s.norm",pchBase);
		vPrintMsg(6,"norm image basename=%s\n",pchNormBase);
	}
	for(iSubset=0; iSubset < iNumSubsets; ++iSubset){
		set_float(pfNormImage, iReconImageSize, 0.0);
		for(iSetMember = 0; iSetMember < iNumAngPerSet; ++iSetMember){
			iAngle = ppiSubsetAngles[iSubset][iSetMember];
			if (iAngle >= iNumAngles)/* this shouldn't happen */
				vErrorHandler(ECLASS_FATAL,ETYPE_ILLEGAL_VALUE,"ComputeNormImages",
		"unexpected angle %d while computing norm matrix, subset %d, member %d:",
					iAngle, iSubset, iSetMember);
			SetMeasuredPrjBinsToOne(pfCalcPrj, pfPrjImage + iPrjSize*iAngle,
				psParms, &psViews[iAngle]);
			/* Taking the time per view into account requires backprojecting
				the time per view instead of all ones.
			*/
			scale_float(pfCalcPrj, iPrjSize, psViews[iAngle].ViewTime);
			BackProject(iNumSlices, iNumPix, iNumRotPix, psParms->iBckModel,
							psParms->PixelWidth, iAngle, 
							psViews[iAngle].Angle, 
							psViews[iAngle].CFCR, psViews[iAngle].psDrfTab, psSrf,
							pfAtnMap, pfCalcPrj, pfNormImage, 1.0, 1.0, 
							psParms->iAxialPadLength, psParms->iAxialAvgLength,
							FALSE, g_sAdvOpts.bFastRotate); 
			vPrintMsg(8,
				"subset=%d, ang=%d, sum(prj)=%.3f,sum(calc)=%.3f, sum(norm)=%.3f\n",
					iSubset, iAngle, 
					sum_float(pfPrjImage +iPrjSize*iAngle,iPrjSize),
					sum_float(pfCalcPrj,iPrjSize),
					sum_float(pfNormImage,iReconImageSize));
			vPrintMsg(7,"  Done angle %d\n",iAngle);
		}
		/* we want to store the image whose pixels are the inverse of the
			norm image. treat all pixel values less than 1e-6 as zeros.
			Note that in osemmw Dan added code that checks for small
			pixel values near the edge of the image. I haven't done that
			here. In the future one might check the first pixel in each row to
			see if it is small and if so, set to zero*/
		for(iPix=0; iPix<iReconImageSize; ++iPix){
			pfNormImage[iPix] = pfNormImage[iPix] < 1e-6 ? 
											0.0 : 1.0/pfNormImage[iPix];
		}
		if (bNormInMemory){
			sum_float(pfNormImages + iReconImageSize*iSubset, iReconImageSize);
			memcpy((void *)(pfNormImages + iReconImageSize*iSubset),
				pfNormImage,sizeof(float)*iReconImageSize);
		}else{
			sprintf(pchNormName,"%s.%d%s",pchNormBase, iSubset,NORM_EXTENSION);
			writedat(pchNormName, iNumPix*iNumPix*iNumSlices*sizeof(float),
				pfNormImage);
		}
		vPrintMsg(6,"Computed norm image for subset %d\n",iSubset);
	}
	psNormData = pvIrlMalloc(sizeof(NormData_t),"ComputeNormImage:NormData");
	psNormData->pfNormImages = pfNormImages;
	psNormData->pchNormImageBase = pchNormBase;
	psNormData->pchNormImageFileName = pchNormName;
	/* save these values from psParms so we only need psNormData to
		get the images*/
	psNormData->bSaveNormImages = FALSE;
	psNormData->bNormInMemory = bNormInMemory;
	psNormData->iNumPix = iNumPix;
	psNormData->iNumSlices = iNumSlices;
	return psNormData;
}

static void vGetNormImage(NormData_t *psNormData, int iSubset, 
					float *pfNormImage)
{
	int iNumSlices = psNormData->iNumSlices;
	int iNumPix = psNormData->iNumPix;
	int iReconImageSize;
	int iDatLen;

	iReconImageSize = iNumPix*iNumPix*iNumSlices;
	if (psNormData->bNormInMemory){
		memcpy(pfNormImage, 
					psNormData->pfNormImages + iSubset*iReconImageSize, 
					iReconImageSize*sizeof(float));
	}else{
		sprintf(psNormData->pchNormImageFileName, "%s.%d%s",
					psNormData->pchNormImageBase, iSubset,NORM_EXTENSION);
		iDatLen=readdat(psNormData->pchNormImageFileName, 0, 
						iReconImageSize*sizeof(float), pfNormImage);
		if (iDatLen != iReconImageSize*sizeof(float))
			vErrorHandler(ECLASS_FATAL,ETYPE_ILLEGAL_VALUE,"GetNormImage",
				"Error reading norm matrix for subset %d from %s",
					iSubset, psNormData->pchNormImageFileName);
	}
}

static void vDoneWithNormImages(int iNumSubsets, NormData_t *psNormData)
/* free memory and delete files, as needed, use to store norm
	images */
{
	int iSubset;

	if (psNormData->bNormInMemory){
		IrlFree(psNormData->pfNormImages);
	}else{
		if (!psNormData->bSaveNormImages){
			vPrintMsg(4,"deleting norm images: base=%s\n",
				psNormData->pchNormImageBase);
			for(iSubset = 0; iSubset < iNumSubsets; ++iSubset){
				sprintf(psNormData->pchNormImageFileName, "%s.%d%s",
							psNormData->pchNormImageBase, iSubset,NORM_EXTENSION);
				vPrintMsg(6,"unlinking :%s:\n",psNormData->pchNormImageFileName);
				unlink(psNormData->pchNormImageFileName);
			}
		}
		IrlFree(psNormData->pchNormImageBase);
		IrlFree(psNormData->pchNormImageFileName);
	}
	IrlFree(psNormData);
}

static int iComputeSubsets(
	Parms_t *psParms, /* uses NumAngles and NumAngPerSet, computes iNumSubsets*/
	int ***pppiSubsetAngles)/* stores pointer to array of pointers to array
										of integes. Each array of integers contains the
										angle indices for a subset. For example,
										(*pppiSubsetAngles)[iSubset][i] will be the
										angle index for the i-th member of subset 
										iSubset*/
/* computes and returns the number of subsets and the angles that belong
	to each subset. */
{
	int iSubset;
	int iMember;
	int iAngle;
	int iNumSubsets;
	int iNumAngles = psParms->NumAngles;
	int iNumAngPerSet = psParms->NumAngPerSubset;
	int **ppiSubsetAngles;
		
	/* NumAngPerSet should have already been checked in
		vGetReconParms, but check it here anyway*/
	if (iNumAngles % iNumAngPerSet)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "iComputeSubsets",
			"%s\n NumAngles=%d, NumAngPerSet=%d",
		"No. of angles not an integer multiple of the no. of angles per subset",
			psParms->NumAngles, psParms->NumAngPerSubset);
	iNumSubsets = iNumAngles / iNumAngPerSet;
	if (iNumSubsets != 1 && iNumSubsets % 2)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "iComputeSubsets",
			"Number of subsets (%d) for num_ang_per_set=%d is not even\n",
			iNumSubsets, psParms->NumAngPerSubset);
	/* allocate memory for the array of pointers*/
	ppiSubsetAngles = (int **)pvIrlMalloc(iNumSubsets * sizeof(int *), 
											"ComputeSubsets:ppiSubsetAngles");
	/* and allocate memory for the angle index array for each subset*/
	for(iSubset = 0; iSubset < iNumSubsets; ++iSubset)
		ppiSubsetAngles[iSubset] = ivector(iNumAngPerSet, 
											"ComputeSubsets:SubsetAngles");
	/* the subset scheme is originally from Dave Lalush. The first angle
		in the first subset is 0. The first angle in the next subset is
		iNumSubsets/2. The next subset is 1, the next subset starts with  
		iNumSubsets/2+1, and so on. The second angle in any subset is
		the first angle plus iNumSubsets. */
	for(iSubset=0; iSubset < iNumSubsets; ++iSubset){
		/* the following is that start angle index for the subset. The
			of the subset number by 2 gives 0 for the first 2 subsets,
			1 for the next 2 subsets. The mod operation is 0 for even
			subsets and 1 for odd subsets. This results in the
			subset scheme mentioned above */
		iAngle = iSubset / 2 + (iSubset % 2)*iNumSubsets/2;
		vPrintMsg(8,"subset %d: ",iSubset);
		for(iMember=0; iMember < iNumAngPerSet; ++iMember){
			ppiSubsetAngles[iSubset][iMember] = iAngle;
			vPrintMsg(8,"%d ", iAngle);
			iAngle += iNumSubsets;
		}
		vPrintMsg(8,"\n");
	}
	vPrintMsg(8,"\n");
	*pppiSubsetAngles = ppiSubsetAngles;
	return iNumSubsets;
}

static void vFreeSubsets(int iNumSubsets, int **ppiSubsetAngles)
{
	int iSubset;

	for(iSubset = 0; iSubset < iNumSubsets; ++iSubset)
		IrlFree(ppiSubsetAngles[iSubset]);
	IrlFree(ppiSubsetAngles);
}

static void FindNextChangeIter(char *pchIterSaveString,
										int *piChangeIter, int *piInterval)
{
	char *pchStart, *pchEnd;
	char ch;

	/*Skip initial whitespace*/
	pchStart = pchIterSaveString + strspn(pchIterSaveString," \t,");
	vPrintMsg(6,"SaveString=%s:\n",pchIterSaveString);
	/*now look for digits*/
	pchEnd=pchStart+strspn(pchStart,"0123456789");
	if (pchEnd == pchStart){
		vErrorHandler(ECLASS_FATAL, ETYPE_SYNTAX, "FindNextChangeIter",
			"missing Change Iteration");
	}
	/* we have a string of digits, now convert it to integer */
	ch = *pchEnd;
	*pchEnd = '\0';
	*piChangeIter=atoi(pchStart);
	*pchEnd = ch;
	pchStart = pchEnd + strspn(pchEnd, " \t,");
	if (*pchStart != '/'){
		if (*pchStart != '\0' && !isdigit(*pchStart)){
			vErrorHandler(ECLASS_FATAL, ETYPE_SYNTAX, "FindNextChangeIter",
				"unexpected error parsing iteration save string: %c",
				*pchStart);
		}
		/* the user did not specify a save interval, assume we save
			only at the next change iteration */
		*piInterval = 0;
	}else{
		/*the user has specified an increment*/
		pchStart++;
		/* skip whitespace */
		pchStart = pchStart + strspn(pchStart," \t,");
		/*now look for digits*/
		pchEnd=pchStart+strspn(pchStart,"0123456789");
		/* we have a string of digits, now convert it to integer */
		ch = *pchEnd;
		*pchEnd = '\0';
		*piInterval=atoi(pchStart);
		*pchEnd = ch;
		/*skip trailing whitespace*/
		pchStart = pchEnd + strspn(pchEnd, " \t,");
	}
	strcpy (pchIterSaveString,pchStart);
}static int bCheckIfSaveIteration(int iIter, Parms_t *psParms)
/* returns true if the image for iteration iIter is to be saved, false if not */
/* the image is saved if it is the last iteration, or if
	the number of iterations since the iteration interval was
	last changed is an integer multiple of the save interval.
	A single save interval can be specified in the SaveInterval
	parameter, or a string of iterations/intervals can be specified
	in the save_iterations parameter. The save iterations is a list
	of iterations to save and save intervals between the iterations.
	For example, 1 2 5 10/5
	Would save iterations 1 2 5 10, and every 5th one after that */
{
	static int iLastIterChange=0;
	static int iCrntInterval=0;
	static int iNextIterChange=0;
	static int iNextInterval=0;
	static int bFirst=TRUE;
	char *pchIterSaveString=psParms->pchIterSaveString;

	if (bFirst){
		bFirst=FALSE;
		if (pchIterSaveString == NULL || *pchIterSaveString == '\0'){
			iNextIterChange = iIter-1;
			iNextInterval=psParms->SaveInterval;
		}
	}
	while (iIter >= iNextIterChange){
		iLastIterChange=iNextIterChange;
		iCrntInterval=iNextInterval;
		if (pchIterSaveString == NULL || *pchIterSaveString == '\0'){
			iNextIterChange = 99999;
		}else{
			FindNextChangeIter(pchIterSaveString,
				&iNextIterChange, &iNextInterval);
		}
		vPrintMsg(6, "Current Save Interval=%d. On iter. %d will change to %d\n",
							iCrntInterval, iNextIterChange, iNextInterval);
	}
	if (iIter == psParms->NumIterations) return (TRUE);
	if (iIter == iLastIterChange) return(TRUE);
	if (iCrntInterval == 0) return (FALSE);
	return((iIter - iLastIterChange) % iCrntInterval == 0);
}

void vOsemRecon(
	Parms_t *psParms, /* input parameters read and set up in routines above*/
	View_t *psViews, /* one view structure for each angle containing things like
							 the COR, ROR, angle, etc.*/
	Srf_t *psSrf, /* parameters and data needed for scatter comp. Null if
							no scatter compensation is needed*/
	float *pfPrjImage, /*modified projection data shifted into array with
								size iNumRotPix and center at iNumRotPix/2
							*/
	float *pfScatterEstimate, /* scatter estimate for "subtraction" scatter
							compensation in "reorderd" format.
							This is NULL if not scatter estimate is to
							be added*/
	float *pfAtnMap, /* Attenuation map in units of per pixel in reordered
							 format. NULL if atten comp is not to be done*/
	float *pfReconImage,
	void (*pIterationCallbackFunc)(int IterationNumber,float *pfCurrentEst))
{
	int iPix;
	float *pfTempImage=NULL;/* used to store normalization matrix and update
										image */
	float *pfPrjScatEst=NULL;
	float *pfCalcPrj;
	NormData_t *psNormImageData;
	int iNumSubsets;
	int iIteration;
	int iSubset;
	int iSetMember;
	int iPrjModel;
	int iNumAngPerSet = psParms->NumAngPerSubset;
	int iAngle;
	int bModelScatThisIter;
	int bModelScatNextIter;
	int bPrjModelSrf = psParms->iPrjModel & MODEL_SRF;
	int **ppiSubsetAngles;/* array of array of pointers to angles in each subset.
							ppiSubsetAngles[iSubset][i] will be the ith angle index to
							use for subset iSubset.*/
	int iNumSlices = psParms->NumSlices;
	int iNumPix = psParms->NumPixels;
	int iNumRotPix = psParms->NumRotPixs;
	int iNumAngles= psParms->NumAngles;
	int iReconImageSize = iNumPix*iNumPix*iNumSlices;
	int iPrjSize = iNumRotPix * iNumSlices;
#ifdef DEBUG
	char pchName[200];
	IMAGE *pPrjImage;
	IMAGE *pRatioPrjImage;
#endif
	float * pfNormImage = vector(iReconImageSize, "vOsemRecon:NormImage");

	PrintTimes("Start OSEM");

	pfCalcPrj = vector(iNumRotPix*iNumSlices, "vOsemRecon:CalcPrj");
	if (bPrjModelSrf && psParms->iNumSrfIterations > 0){
		/* for intermittent RBS we need a place to store the scatter estimate
			for use in the next iteration*/
		pfPrjScatEst = vector(iPrjSize * iNumAngles, "OsemRecon:PrjScatEst");
		set_float(pfPrjScatEst,iPrjSize * iNumAngles, 0.0);
	}
	pfTempImage = vector(iReconImageSize, "vOsemRecon:TempImage");
	iNumSubsets = iComputeSubsets(psParms, &ppiSubsetAngles);
	psNormImageData = psComputeNormImages(psParms->pchNormImageBase, 
								psParms, psViews, psSrf,
								iNumSubsets, ppiSubsetAngles, 
								pfCalcPrj, pfPrjImage, pfAtnMap, pfTempImage);
	PrintTimes("Done Computing Norm Matrix");
	/* if we are modeling the SRF in the projector, set bModelScatNextIter 
		to the proper value for the Start Iteration */
	if (bPrjModelSrf){
		vPrintMsg(7,"modeling scatter: iNumSrfIterations=%d\n",
			psParms->iNumSrfIterations);
		bModelScatNextIter = psParms->StartIteration + 1 <= 
										psParms->iNumSrfIterations;
	}

	g_sOsemHookData.pfEst  = pfReconImage;
	g_sOsemHookData.pfSct  = pfPrjScatEst;
	g_sOsemHookData.pfNorm = pfNormImage;
	g_sOsemHookData.pfUpdate  = pfTempImage;
	psParms->StartIteration = g_sAdvOpts.iStartIteration;
	for(iIteration=psParms->StartIteration+1; 
			iIteration<= psParms->NumIterations;
			++iIteration){
		/* effects to model this iteration. Update this if scatter isn't modeled
			this iteration */
		
		// Call Prev-Iteration callback function if defined
		g_sOsemHookData.iIterationIdx = iIteration;
		if (g_pPrevIterHookFunc != NULL)
			(*g_pPrevIterHookFunc)(&g_sOsemHookData, g_pCtx);
		
		vPrintMsg(8,"sum of recon image at start of iteration=%.2f\n",
			sum_float(pfReconImage, iReconImageSize));
#ifdef DEBUG
#endif
		iPrjModel = psParms->iPrjModel;
		if (bPrjModelSrf){
			/* if we are modeling the SRF, check whether we need it on or off
				for this iteration */
			bModelScatThisIter = bModelScatNextIter;

			if (iIteration == psParms->NumIterations){
				/* if this is the last iteration, there is no need to figure 
					out if we are modeling it the next iteration, so keep the current
					value. */
				bModelScatNextIter = bModelScatThisIter;
			}else{
				/* this isn't the last iteration, so figure out if we are
					modeling it on the next iteration*/
				bModelScatNextIter = iIteration + 1 <= psParms->iNumSrfIterations;
			}
			vPrintMsg(6,"Iteration %d: this iteration: %s model scatter\n",
				iIteration, bModelScatThisIter ? "do" : "don't");
			vPrintMsg(6,"              next iteration: %s model scatter\n",
				bModelScatNextIter ? "do" : "don't");
			/* update iPrjModel flag for this iteration to scatter on
				or off as needed*/
			iPrjModel = bModelScatThisIter ? psParms->iPrjModel :
														psParms->iPrjModel & ~MODEL_SRF;
		}
#ifdef DEBUG
		sprintf(pchName,"prj.%d.im",iIteration);
		vPrintMsg(6,"prj file name=%s\n",pchName);
		pPrjImage=
			imgio_openimage(pchName,'n',&iNumRotPix, &iNumSlices, &iNumAngles);
		sprintf(pchName,"ratioprj.%d.im",iIteration);
		pRatioPrjImage=
			imgio_openimage(pchName,'n',&iNumRotPix, &iNumSlices, &iNumAngles);
		vPrintMsg(6,"ratioprj file name=%s\n",pchName);
#endif
		for(iSubset=0; iSubset < iNumSubsets; ++iSubset){
			// Call Pre-Subset callback function if defined
			g_sOsemHookData.iSubsetIdx = iSubset;
			if (g_pPrevSubsetHookFunc != NULL)
				(*g_pPrevSubsetHookFunc)(&g_sOsemHookData, g_pCtx);
			
			/* zero the update image for the subset*/
			set_float(pfTempImage, iReconImageSize, 0.0);
			for(iSetMember = 0; iSetMember < iNumAngPerSet; ++iSetMember){
				iAngle = ppiSubsetAngles[iSubset][iSetMember];
				
				// Call Pre-Angle callback function if defined
				g_sOsemHookData.iAngleIdx = iAngle;
				if (g_pPrevAngleHookFunc != NULL)
					(*g_pPrevAngleHookFunc)(&g_sOsemHookData, g_pCtx);
				
				if (iAngle >= iNumAngles)/* this shouldn't happen */
					vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "OsemRecon",
						"unexpected angle %d in iteration %d, subset %d, member %d:",
						iAngle, iIteration, iSubset, iSetMember);
				/* project into pfCalcPrj */
				vPrintMsg(8,"model srf=%d\n",iPrjModel&MODEL_SRF);
				Project(iNumSlices, iNumPix, iNumRotPix, iPrjModel,
							psParms->PixelWidth,
							iAngle, 
							psViews[iAngle].Angle,
							psViews[iAngle].CFCR,
							psViews[iAngle].psDrfTab, psSrf, 
							pfAtnMap, pfReconImage, pfCalcPrj, 
							1.0, 1.0, 
							psParms->iAxialPadLength, psParms->iAxialAvgLength,
							FALSE, FALSE, g_sAdvOpts.bFastRotate);
				/* scale the calculated projection by the time per view.
					Note that we do not do this before the backprojection
					because if done both in the calculation of the normalization
					matrix and the backprojection of the ratios the factor
					will cancel out. Also note that the scatter estimate
					is not scaled by the time because it is assumed that
					the input estimate has already been scaled 
				*/
				scale_float(pfCalcPrj, iPrjSize, 
					psViews[iAngle].ViewTime*psViews[iAngle].Sensitivity);
#ifdef DEBUG
				vPrintMsg(8,"prj %d imgsum=%.6g, prjsum=%.6g\n", iAngle,
								sum_float(pfReconImage, iReconImageSize),
								sum_float(pfCalcPrj, iPrjSize));
#endif
				if (bPrjModelSrf){
					if (!bModelScatThisIter){
						/*if we did not model scatter this iteration, 
							we add the previous scatter estimate to the 
							projection data */
#ifdef DEBUG
						vPrintMsg(8,"sum before adding scatter estimate=%.2f\n",
							sum_float(pfCalcPrj, iPrjSize));
#endif
						for(iPix=0; iPix < iPrjSize; ++iPix){
							pfCalcPrj[iPix] += pfPrjScatEst[iPrjSize*iAngle + iPix];
						}
#ifdef DEBUG
						vPrintMsg(8,"sum after adding scatter estimate=%.2f\n",
							sum_float(pfCalcPrj, iPrjSize));
#endif
					}else if (!bModelScatNextIter){
						/* Here we have modeled scatter this iteration, but are
							not modeling it next iteration. So we need to perform
							the projection again with scatter modeling turned off.
							Instead of projecting into CalcPrj, project into
							PrjScatEst and the subtract this from CalcPrj to get
							the scatter estimate*/
						Project(iNumSlices, iNumPix, iNumRotPix, 
							psParms->iPrjModel & ~MODEL_SRF,
							psParms->PixelWidth,
							iAngle, 
							psViews[iAngle].Angle,
							psViews[iAngle].CFCR,
							psViews[iAngle].psDrfTab, psSrf, 
							pfAtnMap, pfReconImage, 
							pfPrjScatEst+iPrjSize*iAngle, 
							1.0, 1.0, 
							psParms->iAxialPadLength, psParms->iAxialAvgLength,
							TRUE, /* no need to rotate attenuation map*/
							TRUE, /* or the ReconImage)*/
							g_sAdvOpts.bFastRotate); /* Fast Rotate */
						/* Scale by time per view */
						scale_float(pfPrjScatEst+iPrjSize*iAngle, iPrjSize, 
								psViews[iAngle].ViewTime);
						/* and compute the scatter estimate by subtracting
						 	the no scatter projection from the projection with
							scatter */
#ifdef DEBUG
						vPrintMsg(8,"sum with scatter=%.2f\n", 
							sum_float(pfCalcPrj, iPrjSize));
						vPrintMsg(8,"sum w/o scatter=%.2f\n", 
							sum_float(pfPrjScatEst+iPrjSize*iAngle, iPrjSize));
#endif
						for(iPix=0; iPix < iPrjSize; ++iPix){
							pfPrjScatEst[iPrjSize*iAngle + iPix] = 
								pfCalcPrj[iPix] - pfPrjScatEst[iPrjSize*iAngle + iPix];
						} 
					}
				}
				if (pfScatterEstimate != NULL){
					/* there is a scatter component to add to the projection data*/
					for(iPix=0; iPix < iPrjSize; ++iPix){
						pfCalcPrj[iPix] += 
							pfScatterEstimate[iPrjSize*iAngle + iPix];
					}
				}
#ifdef DEBUG
				imgio_writeslices(pPrjImage, iAngle,iAngle,pfCalcPrj);
#endif
				/* take ratio of measured and calculated projections */
				/* set non-measured bins to zero so we ignore them when we take
					the ratio*/
				SetNonMeasuredBinsToZero(pfCalcPrj, pfPrjImage, 
					psParms, psViews + iAngle);
				for(iPix=0; iPix < iPrjSize; ++iPix){
					/* only take ratio for positive values of CalcPrj. That way
						we will set ratio for non-measured bins to zero and they
						won't affect update*/
					if(pfCalcPrj[iPix] > 1e-20 ){
						pfCalcPrj[iPix] = pfPrjImage[iPrjSize*iAngle + iPix] /
												pfCalcPrj[iPix];
					}else{
						pfCalcPrj[iPix] = 0.0;
					}
				}
				/* pfCalcPrj now contains the error projections (ratio of measured
					and estimated projections. We backproject this into TempImage.*/
				/* Taking into account the time per view requires scaling the
					backprojected projection by the acquisition time per view
				*/
				scale_float(pfCalcPrj, iPrjSize, psViews[iAngle].ViewTime);
				BackProject(iNumSlices, iNumPix, iNumRotPix, psParms->iBckModel,
								psParms->PixelWidth, iAngle, 
								psViews[iAngle].Angle, 
								psViews[iAngle].CFCR, psViews[iAngle].psDrfTab, psSrf,
								pfAtnMap, pfCalcPrj, pfTempImage, 1.0, 1.0, 
								psParms->iAxialPadLength, psParms->iAxialAvgLength,
								TRUE, g_sAdvOpts.bFastRotate);
								/* since we just performed projection above,
										 the rotated atteunation map is for the correct
										 angle so we tell BackProject that the
										 attenuation map is already rotated*/

				// Call Post-Angle callback function if defined
				if (g_pPostAngleHookFunc != NULL)
					(*g_pPostAngleHookFunc)(&g_sOsemHookData, g_pCtx);
			}/* iSetMember Loop*/
			
			/* now we need to get the normalization image for this angle. We
				actually get the element-by-element	inverse of the normalization
				image so we don't need to check	for zeros if we used division*/
			vGetNormImage(psNormImageData, iSubset, pfNormImage);
#ifdef DEBUG
			sprintf(pchName,"norm.%d.im",iSubset);
			writeimage(pchName, iNumPix, iNumSlices, iNumPix, pfNormImage);
#endif
			/* do the pixel by pixel multiplication of the current estimate
				(in pfReconImage) by the update image (in pfTempImage). 
			    And do the normalization by multiplying (since the elements in
				pfNormImage are the inverse of the elements in the norm image*/
			for(iPix=0; iPix < iReconImageSize; ++iPix){
				pfReconImage[iPix] *= pfTempImage[iPix] * pfNormImage[iPix];
			}
#ifdef WDEBUG
			sprintf(pchName,"partupdate.%d.im",iSubset);
			writeimage(pchName, iNumPix, iNumSlices, iNumPix, pfReconImage);
#endif

			// Call Post-Subset callback function if defined
			if (g_pPostSubsetHookFunc != NULL)
				(*g_pPostSubsetHookFunc)(&g_sOsemHookData, g_pCtx);
			
			vPrintMsg(6,"Done subset %d\n",iSubset);
		} /* iSubset Loop */
#ifdef DEBUG
		imgio_closeimage(pPrjImage);
		imgio_closeimage(pRatioPrjImage);
#endif
		
		// Call Post-Iteration callback function if defined
		if (g_pPostIterHookFunc != NULL)
			(*g_pPostIterHookFunc)(&g_sOsemHookData, g_pCtx);
		
		/* done the iteration. Check to see if it needs to be saved */
		if (bCheckIfSaveIteration(iIteration, psParms)){
			reorder(iNumPix, iNumSlices, iNumPix, pfReconImage, pfTempImage);
			if (pIterationCallbackFunc != NULL)
				(*pIterationCallbackFunc)(iIteration, pfTempImage);
		}
		
		vPrintMsg(4,"----Done iteration %d\n", iIteration);
		PrintTimes("  iteration times:");
		vPrintMsg(4,"\n");
	} /* iIteration Loop */
	/* clean up*/
	/* free up memory and delete disk files used to store normalization matrices,
		as needed*/
	vPrintMsg(8,"sum of pfRecon after all iterations=%.4g (before reorder) ",
		sum_float(pfReconImage, iNumSlices*iNumPix*iNumPix));
	vDoneWithNormImages(iNumSubsets, psNormImageData);
	vFreeSubsets(iNumSubsets, ppiSubsetAngles);
	IrlFree(pfCalcPrj);
	IrlFree(pfTempImage);
	IrlFree(pfPrjScatEst);
	IrlFree(pfNormImage);
	vPrintMsg(8,", %.4g (after freeing memory)\n",
		sum_float(pfReconImage, iNumSlices*iNumPix*iNumPix));
}
