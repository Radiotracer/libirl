/*
	@file irlgenprj.c
	
	@code $Id: irlgenprj.c 0 2005-05-14 21:29:05Z duyong $ @endcode
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

#ifndef FASTROTATE
#define FASTROTATE 1
#endif

int IrlGenprj(
	IrlParms_t *psIrlParms, /* projection/reconstruction parameters, see below */
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
	void (*pPrjViewCallbackFunc)(int PrjViewNumber, float *pfCurrentPrjView),
						       /*callback function that is called at the end of each 
									projection view. The callback function has two arguments: 
									the view number and a pointer to the Current projection 
									view image. If this is Null, then no callback will be made 
								 */
	float *pfScatterEstimate, /* Estimate of scatter that will be added to
										  the the computed projection. Must be the same
										  size and pixel order as the
										  projection image. If NULL, then no scatter
										  estimate is used.
										  Notes: 1. The estimate is scaled by fScatEstFac
										            before being added.
													2. This is indendent of scatter modeling
													   and can be used in addition or in
														lieu of it.
													3. Currently only supports same bin size as
														pixel size. The interpolation from different
														bin size is expected to be done by user
														before calling irlGenprj
									  */
	float *pfActImage,			/* pointer to image of activity distribution
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
											transaxial bin, slice, view. Currently only 
											supports using of bin-size equals pixel-size.
											For different bin and pixel size, the interpolation
											is expected to be done by user after calling irlGenprj 
										*/
	float fPrimaryFac,			// Primary factor, use 0 when generate scatter estimation only
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
	int iNumActPixels;
   int i, iAngle;
   float *pfCalcPrj=NULL;

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
	}else if (strcmp(pchLogFname,"stderr")){
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
	}else if (strcmp(pchMessageFname,"stderr")){
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

	/* Set up the projection parameters array, sParms */
	vSetupParms(psIrlParms, psOptions, 
		pfPrjImage, pfAtnMap, pfActImage, pfScatterEstimate, 
			pchScatKrnlFname, &sParms);
	psViews=psSetupViews(psIrlParms, psOptions, psPrjViews);

	
//	fprintf(stderr,"model=%d, prjmodel=%d, bckmodel=%d\n",
//		sParms.iModel, sParms.iPrjModel, sParms.iBckModel);
	if (sParms.iModel & MODEL_SRF){
	  psSrf = psSrfSetup(&sParms, psIrlParms->SrfCollapseFac, pchScatKrnlFname);
	}else{
		psSrf=NULL;
	}
	if (sParms.iModel & MODEL_DRF)
		vSetupDrfTables(&sParms, psOptions, psViews, pchDrfTabFname);

	/* reorder the activity image to the internally-used format */
	pfTemp=vector(sParms.NumPixels*sParms.NumPixels*sParms.NumSlices,
								"IrlOsem:Temp");
	reorder(sParms.NumSlices, sParms.NumPixels, sParms.NumPixels,
					pfActImage, pfTemp);
	memcpy(pfActImage, pfTemp, 
			sParms.NumPixels*sParms.NumPixels*sParms.NumSlices*sizeof(float));
	IrlFree(pfTemp);

	/*vInitPrj sets up temporary variables used by the projector */
	InitPrj(sParms.NumSlices, sParms.NumPixels, sParms.NumRotPixs, 
		sParms.iModel);
	
	/*irl_genprj projects and saves the projections for each view */
	PrintTimes("Done Initialization");

	iNumActPixels=sParms.NumPixels*sParms.NumPixels*sParms.NumSlices;
	iNumPrjPixels=sParms.NumPixels*sParms.NumSlices*sParms.NumAngles;

	/* make sure there are no negative pixels in the active images */
	vZeroNeg(pfActImage, iNumActPixels);

	/* scale the attenuation map and zero any negative values in it */
	if (pfAtnMap != NULL){
		if (psIrlParms->fAtnScaleFac > 0 && psIrlParms->fAtnScaleFac != 1.0){
			scale_float(pfAtnMap, iNumActPixels, psIrlParms->fAtnScaleFac);
		}
		vZeroNeg(pfAtnMap, iNumActPixels);
    /* reorder the attenuation map to the internally-used format */
    pfTemp=vector(iNumActPixels,"IrlOsem:Temp");
    reorder(sParms.NumSlices, sParms.NumPixels, sParms.NumPixels,
          pfAtnMap, pfTemp);
    memcpy(pfAtnMap, pfTemp, 
      sParms.NumPixels*sParms.NumPixels*sParms.NumSlices*sizeof(float));
    IrlFree(pfTemp);
	}
	/* scale the scatter estimate, if provided, and zero the negatives in it */
	if (pfScatterEstimate != NULL){
		if (psIrlParms->fScatEstFac > 0 && psIrlParms->fScatEstFac != 1.0){
			scale_float(pfScatterEstimate, iNumPrjPixels, 
				psIrlParms->fScatEstFac);
		}
		vZeroNeg(pfScatterEstimate, iNumPrjPixels);
	}

	/* set all pixels in projections to zero */
   set_float(pfPrjImage, iNumPrjPixels, 0.0);

	/* return the initial estimate to the callback function */
	pfTemp=vector(sParms.NumPixels*sParms.NumPixels*sParms.NumSlices,
							"IrlOsem:Temp");
	reorder(sParms.NumPixels, sParms.NumSlices, sParms.NumPixels,
				pfActImage, pfTemp);
	IrlFree(pfTemp);

	/* do the projection */
   /* Here the projection image dimension should be same as activity images */
   /* interpolation from pixel-size projection to bin-size projection should 
      done in application level by users */
   pfCalcPrj = vector(sParms.NumRotPixs * sParms.NumSlices,"IrlGenprj:CalcPrj");

   for(iAngle = 0; iAngle < sParms.NumAngles; iAngle++){
      set_float(pfCalcPrj, sParms.NumRotPixs * sParms.NumSlices, 0.0);
      vPrintMsg(6,"Projecting angle %d\n",iAngle);
      Project(sParms.NumSlices, sParms.NumPixels, sParms.NumRotPixs, 
              sParms.iPrjModel,
              sParms.PixelWidth, 
              iAngle,
              psViews[iAngle].Angle, 
              psViews[iAngle].CFCR,
              psViews[iAngle].psDrfTab, psSrf,
              pfAtnMap, pfActImage, 
              pfCalcPrj,
              fPrimaryFac,   /* the primary factor is handled in scatter est. */
              1.0,   /* scaling of attenuation map is done when reading in image */ 
              sParms.iAxialPadLength, sParms.iAxialAvgLength,
              FALSE, FALSE,  /* need to rotate act and atn for each angle */
              FASTROTATE );
       vPrintMsg(6,"   Done angle %d\n", iAngle);
       if (pfScatterEstimate != NULL) {
          /* there is a scatter component to add to the projection data*/
          for(i=0; i < sParms.NumRotPixs * sParms.NumSlices; ++i){
              pfCalcPrj[i] +=
                     pfScatterEstimate[sParms.NumRotPixs * sParms.NumSlices*iAngle + i];
               }
       }
       /* save calculated projection view into output projection marix */
		 /* scale by time per view and Sensitivity and copy to output image*/
		 for(i=0; i < sParms.NumRotPixs * sParms.NumSlices; ++i){
			 pfPrjImage[sParms.NumRotPixs * sParms.NumSlices*iAngle + i] = 
			 	pfCalcPrj[i]*psViews[iAngle].ViewTime*psViews[iAngle].Sensitivity; 
       }
       if(pPrjViewCallbackFunc != NULL)
           (*pPrjViewCallbackFunc)(iAngle, 
			  			pfPrjImage+sParms.NumRotPixs*sParms.NumSlices*iAngle); 
       
   }

	scale_float(pfPrjImage, iNumPrjPixels, 1.0/(float)sParms.NumAngles);
	vPrintMsg(8,"Total counts in pfPrjImage= %.4g (after Genprj)\n",
		sum_float(pfPrjImage,
		sParms.NumSlices*sParms.NumPixels*sParms.NumAngles));
	/* reorder the activity image and attenuation map */
	pfTemp=vector(sParms.NumPixels*sParms.NumPixels*sParms.NumSlices,
							"IrlOsem:Temp");
	reorder(sParms.NumPixels, sParms.NumSlices, sParms.NumPixels,
				pfActImage, pfTemp);
	memcpy(pfActImage, pfTemp, 
		sParms.NumPixels*sParms.NumPixels*sParms.NumSlices*sizeof(float));
   if(pfAtnMap != NULL){
      reorder(sParms.NumPixels, sParms.NumSlices, sParms.NumPixels,
              pfAtnMap, pfTemp);
      memcpy(pfAtnMap, pfTemp,
             sParms.NumPixels*sParms.NumPixels*sParms.NumSlices*sizeof(float));
   }
	IrlFree(pfTemp);
	vPrintMsg(4,"done\n");
abort:
	IrlFree(pfCalcPrj);
	IrlFree(psViews);
	psViews=NULL;
	vFreeExtraPrjViewData();
	vFreeDetectors();
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
