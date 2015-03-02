/*
	@file irlosem.c
	
	@code $Id: irlosem.c 121 2011-09-14 02:57:26Z frey $ @endcode
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>

#include <mip/errdefs.h>
#include <mip/printmsg.h>
#include <mip/realft.h>
#include <mip/miputil.h>

#include "irlprivate.h"


int IrlOsem(
	IrlParms_t *psIrlParms, /* reconstruction parameters, see below */
	Options_t *psOptions, /*structure containing reconstruction options */
	PrjView_t *psPrjViews, /*array of View_t structures, one for each projection
							 angle (see below) in the same order as the actual
							 images */
	char *pchDrfTabFname, /*filename containing the drf table. If NULL, and
									drf modeling is requested, then the grf is
									computed on the fly using the parameters contained
									in the IrlParms structure (fHoleDiam, fHoleLen,
									fBackToDet, and fIntrinsicFwhm). Note that the
									DRF table file must have the same pixel size as
									the reconstructed image.
								*/
									
	char *pchScatKrnlFname, /*filename containing the scatter kernel file.
									  if NULL, then scatter is not modeled.
									  The pixel. The scatter kernel should be for the
									  appropriate isotope and must have the same pixel
									  size as the reconstructed image after factoring
									  in the SrfCollapseFac. For example, if the
									  pixel size in the  reconstructed image is 0.6 cm,
									  and a collapse factor of 2 is used, then the
									  pixel size in the scatter kernel must be 1.2 cm
									*/
	void (*pIterationCallbackFunc)(int IterationNumber, float *pfCurrentEst),
	       /*callback function that is called at the end of every iteration.
			   The callback function has two arguments: the iteration number
				and a pointer to the Current Estimate image. If this is Null,
				then no callback will be made */
	float *pfScatterEstimate, /* Estimate of scatter that will be added to
										  the reconstructed image estimate after
										  the computed projection. Must be the same
										  size and pixel order as the
										  projection image. If NULL, then no scatter
										  estimate is used.
										  Notes: 1. The estimate is scaled by fScatEstFac
										            before being added.
													2. This is indendent of scatter modeling
													   and can be used in addition or in
														lieu of it.
									  */
	float *pfAtnMap,				/* pointer to attenuation map image. This
											must be the same size as the reconstructed
											image (NumPixels*NumPixels*NumSlices).
											The pixel ordering is (fastest to slowest)
											x, y, z (slice).
										*/
	float *pfPrjImage,		 	/* pointer to projection image. This must have
											the size NumPixels*NumPixels*NumSlices. The
											pixel ordering is (fastest to slowest):
											transaxial bin, slice, view
										*/
	float *pfReconImage,			/* pointer to reconstructed image initial
											estimate (if bReconIsInitEst is true) and
										 	memory where the reconstructed image will
											be stored.
										*/
	char *pchLogFname,			/* file name where warning, information, and
	                              fatal error messages are written. If NULL,
											this will not be created */
	char *pchMessageFname		/* file name where status messages are printed.
											if Null, these are not printed */

										  
)
{
	int err_num=0;
	jmp_buf jmp_env;
	float *pfTemp=NULL;

	Parms_t sParms;
	Srf_t *psSrf=NULL;
	FILE *fpLogFile=NULL, *fpMessageFile=NULL;
	View_t *psViews=NULL;
	int iNumPrjPixels;
	int iNumReconPixels;


	/* set up to trap fatal errors. Note that the longjmp that uses this
		setjmp environment is _only_ called from the ErrorHandler routine.
	*/
	if ((err_num=setjmp(jmp_env))){
			/* tells ErrorHandler to jump to above on fatal error */
			fprintf(stderr,"FATAL ERROR Trapped: err_num=%d\nErrMsg=%s\n",
			err_num, pchErrorMsg());
		vPrintMsg(1,"aborted\n");
		goto abort;
	}
	vTrapFatalErrs(&jmp_env); 

	if (pchLogFname == NULL){
		vSetErrLog(NULL);
	}else if (strcmp(pchLogFname,"stderr")==0){
		vSetErrLog(stderr);
	}else{
		fpLogFile=fopen(pchLogFname,"w");
		if (fpLogFile == NULL){
			vSetErrLog(NULL);
			vErrorHandler(ECLASS_FATAL, ETYPE_IO, 
				"IrlOsem","Error opening log file %s", pchLogFname);
			goto abort;
		}
		vSetErrLog(fpLogFile);
	}

	vSetMsgLevel(psOptions->iMsgLevel);
	if (pchMessageFname == NULL){
		vSetAllMsgFilePtr(NULL);
	}else if (strcmp(pchMessageFname,"stderr")==0){
		vSetAllMsgFilePtr(stderr);
	}else{
		fpMessageFile=fopen(pchMessageFname,"w");
		if (fpMessageFile == NULL){
			vErrorHandler(ECLASS_FATAL, ETYPE_IO, 
				"IrlOsem","Error opening messages file %s", pchMessageFname);
			goto abort;
		}
		vSetAllMsgFilePtr(fpMessageFile);
	}

	/* Set up the reconstruciton parameters array, sParms */
	vSetupParms(psIrlParms, psOptions, 
		pfPrjImage, pfAtnMap, pfReconImage, pfScatterEstimate, 
			pchScatKrnlFname, &sParms);
	psViews=psSetupViews(psIrlParms, psOptions, psPrjViews);

	
	vPrintMsg(6,"model=%d, prjmodel=%d, bckmodel=%d\n",
		sParms.iModel, sParms.iPrjModel, sParms.iBckModel);
	if (sParms.iModel & MODEL_SRF){
	  psSrf = psSrfSetup(&sParms, psIrlParms->SrfCollapseFac, pchScatKrnlFname);
	}else{
		psSrf=NULL;
	}

	if (sParms.iModel & MODEL_DRF)
		vPrintMsg(6,"SetupDrfTables\n");
		vSetupDrfTables(&sParms, psOptions, psViews, pchDrfTabFname);
	{
		int iAng;
		for(iAng=0; iAng<sParms.NumAngles; ++iAng)
			vPrintMsg(8,"  %d=%p\n",iAng,psViews[iAng].psDrfTab);
	}
	
	iNumReconPixels=sParms.NumPixels*sParms.NumPixels*sParms.NumSlices;
	iNumPrjPixels=sParms.NumPixels*sParms.NumSlices*sParms.NumAngles;

    /* reorder the attenuation map to the internally-used format */
	if (pfAtnMap != NULL){
		pfTemp=vector(iNumReconPixels,"IrlOsem:Temp");
		reorder(sParms.NumSlices, sParms.NumPixels, sParms.NumPixels,
				pfAtnMap, pfTemp);
		memcpy(pfAtnMap, pfTemp, 
		  sParms.NumPixels*sParms.NumPixels*sParms.NumSlices*sizeof(float));
		IrlFree(pfTemp);
	}

	if (! psOptions->bReconIsInitEst){
		/* compute or get initial estimate */
		pfReconImage=pfComputeInitialEstimate(&sParms, 
								psOptions->bUseContourSupport,
								psOptions->fAtnMapThresh, psViews,
								pfPrjImage, pfAtnMap, pfReconImage);
	}else{
		/* reorder the initial estimate to the internally-used format */
		pfTemp=vector(sParms.NumPixels*sParms.NumPixels*sParms.NumSlices,
								"IrlOsem:Temp");
		reorder(sParms.NumSlices, sParms.NumPixels, sParms.NumPixels,
					pfReconImage, pfTemp);
		memcpy(pfReconImage, pfTemp, 
			sParms.NumPixels*sParms.NumPixels*sParms.NumSlices*sizeof(float));
		IrlFree(pfTemp);
	}

	/*vInitPrj sets up temporary variables used by the projector and
	  backprojector*/
	InitPrj(sParms.NumSlices, sParms.NumPixels, sParms.NumRotPixs, 
		sParms.iModel);
	
	/*irl_osem computes subsets,normalization matrices. It projects+backprojects
	  It saves the reconstruction from each iteration and all the projections 
	  from the current iteration*/
	PrintTimes("Done Initialization");


	/* make sure there are no negative pixels in the projections */
	vZeroNeg(pfPrjImage, iNumPrjPixels);

	/* scale the attenuation map and zero any negative values in it */
	if (pfAtnMap != NULL){
		if (psIrlParms->fAtnScaleFac > 0 && psIrlParms->fAtnScaleFac != 1.0){
			scale_float(pfAtnMap, iNumReconPixels, psIrlParms->fAtnScaleFac);
		}
		vZeroNeg(pfAtnMap, iNumReconPixels);
	}

	/* scale the scatter estimate, if provided, and zero the negatives in it */
	if (pfScatterEstimate != NULL){
		if (psIrlParms->fScatEstFac > 0 && psIrlParms->fScatEstFac != 1.0){
			scale_float(pfScatterEstimate, iNumPrjPixels, 
				psIrlParms->fScatEstFac);
		}
		vZeroNeg(pfScatterEstimate, iNumPrjPixels);
	}
	/* make sure that that there are no negatives in the initial estimate
		(which is stored in the reconstruction image */
	vZeroNeg(pfReconImage, iNumReconPixels);
	/* return the initial estimate to the callback function */
	pfTemp=vector(sParms.NumPixels*sParms.NumPixels*sParms.NumSlices,
							"IrlOsem:Temp");
	reorder(sParms.NumPixels, sParms.NumSlices, sParms.NumPixels,
				pfReconImage, pfTemp);
	(*pIterationCallbackFunc)(0, pfTemp);
	IrlFree(pfTemp);

	/* do the reconstruction. Input images are the prj image, atn map and
		the initial estimate (passed via pfReconImage)*/

	vOsemRecon(&sParms, psViews, psSrf, pfPrjImage, pfScatterEstimate, 
				pfAtnMap, pfReconImage, pIterationCallbackFunc);
	vPrintMsg(8,"Total counts in pfReconImage= %.4g (after vOsemRecon)\n",
		sum_float(pfReconImage,
		sParms.NumSlices*sParms.NumPixels*sParms.NumPixels));
	/* reorder the output image */
	pfTemp=vector(sParms.NumPixels*sParms.NumPixels*sParms.NumSlices,
							"IrlOsem:Temp");
	reorder(sParms.NumPixels, sParms.NumSlices, sParms.NumPixels,
				pfReconImage, pfTemp);
	memcpy(pfReconImage, pfTemp, 
		sParms.NumPixels*sParms.NumPixels*sParms.NumSlices*sizeof(float));
	vPrintMsg(8,"Total counts in pfReconImage= %.4g (after reorder)\n",
		sum_float(pfReconImage,
		sParms.NumSlices*sParms.NumPixels*sParms.NumPixels));
	IrlFree(pfTemp);
	vPrintMsg(4,"done\n");
abort:
	vFreeExtraPrjViewData();
	vFreeDetectors();

	IrlFree(psViews);
	psViews=NULL;
	vPrintMsg(6,"freeing Srf Data\n");
	psSrf=psFreeSrf(&sParms, psSrf);
	vPrintMsg(6,"freeing Parameters\n");
	FreeParms(&sParms);
	vPrintMsg(6,"freeing Projectors\n");
	FreePrj();

	vPrintMsg(6,"freeing fft memory\n");
	realft_free();
//	vPrintMsg(6,"freeing everything else\n");
//	vFreeAllMem();
	vPrintMsg(7,"done freeing memory\n");
	if (fpLogFile != NULL){
		fclose(fpLogFile);
		vSetErrLog(NULL);
	}
	if (fpMessageFile != NULL){
		fclose(fpMessageFile);
		vSetAllMsgFilePtr(NULL);
	}
	return(err_num);
}

char *pchIrlErrorString(void)
{
	return pchErrorMsg();
}
