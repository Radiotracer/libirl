/*
	prjbck.c

	$Id: prjbck.c 114 2011-01-31 16:51:08Z frey $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <math.h>

#include <mip/miputil.h>
#include <mip/printmsg.h>

#include "krnlprj.h"
#include "irlprivate.h"

#define EPS (1.0e-5)

#define _3dptr(iRow, iSlice, iCol) ((iRow)*iNumSlices*iNumRotPixs+\
													(iSlice)*iNumRotPixs+(iCol))
#define _2dptr(iSlice, iCol) ((iSlice)*iNumRotPixs+(iCol))
#define _23dptr(iRow, i) ((iRow)*iNumSlices*iNumRotPixs+(i))

/* Init mallocs space for matrices being used in projection or backprojection
	It uses user-defined options like attenuation or scatter to decide whether 
	to allocate memory for some of the matrices. All the matrices are in the
	structure sTempImages*/
static struct {
	float *RotReconImage;
	float *UnRotReconImage;
	float *RotAtnImage;
	float *AtnImage;
	float *ConvImage;
	float *ScatterEstimate;
	float *SumAtnTable;
} sTempImages = {NULL,NULL,NULL,NULL,NULL,NULL};

extern float *pfgConvMatrix;
extern int igConvMatrixSize;

void InitPrj(int iNumSlices, int iNumPixels, int iNumRotPixs, int iModel)
{	
	sTempImages.RotReconImage = vector(iNumRotPixs*iNumSlices*iNumRotPixs,
		"RotRecon:Init");
	sTempImages.UnRotReconImage=vector(iNumRotPixs*iNumSlices*iNumRotPixs,
		"UnRotRecon:Init");
	
	if(iModel & MODEL_ATN || iModel & MODEL_SRF){
		sTempImages.RotAtnImage = vector(iNumRotPixs*iNumSlices*iNumRotPixs,
			"RotAtn:Init");
		sTempImages.SumAtnTable = 
			vector(iNumSlices*iNumRotPixs,"SumAtnTable:Init");
		sTempImages.AtnImage = 
			vector(iNumSlices*iNumRotPixs,"SumAtnTable:Init");
	}
	
	if(iModel & MODEL_DRF)
		sTempImages.ConvImage = vector(iNumSlices*iNumRotPixs,"ConvImage:Init");
	
	if(iModel & MODEL_SRF)
		sTempImages.ScatterEstimate=vector(iNumSlices*iNumRotPixs*iNumRotPixs,
									"ScatterEstimate:Init");

	vInitFastRotate(iNumPixels);
}

void Project(
	int iNumSlices,
	int iNumPixels,
	int iNumRotPixs,
	int iPrjModel,
	float fPixelWidth,
	int iAngle,
	float fAngle,
	float fCFCR,
	DrfTab_t *psDrfTab,
	Srf_t *psSrf, 
	float *atnmatrix, 
	float *reconmatrix, 
	float *calcprj, 
	float fPrimaryFac, 
	float fAtnFac, 
	int iAxialPadLength, // Number of slices to pad at top and bottom to reduce
								// Edge effects (passed to ConvolveDrfYpad)
	int iAxialAvgLength, //Number of slices to average over when performing
								//axial padding
	int bAtnMapRotated, 
	int bActImgRotated,
	int bFastRotate)
	/* note that after exit the rotated image and rotated attenuation maps
		are still stored in sTempImages.pfRotReconImage and 
		sTempImages.pfRotAtnImage, respectively. This way if the same angle
		is to be projected again (as when we are computing the scatter estimate
		when we are not modeling scatter next iteration) we won't need to
		rotate the ReconImage. Also, when we are backprojecting the same
		angle we just projected we won't have to rotate the attenuation map
		(if attenuation or scatter are modeled)
	*/
{
	int  iRow, i; 
	int iDist;
	int bModelAtn = iPrjModel & MODEL_ATN;
	int bModelSrf = iPrjModel & MODEL_SRF;
	int bModelDrf = iPrjModel & MODEL_DRF;
	float *pfProjPixels; 
	float *pfRotPixels;
	float fStartDist;
	float fDist;
	float *pfSumAtnTable = sTempImages.SumAtnTable;
	float *pfRotAtnImage = sTempImages.RotAtnImage;
	float *pfRotReconImage = sTempImages.RotReconImage;
	float *pfScatterEstimate = sTempImages.ScatterEstimate;
	float *pfAtnImage = sTempImages.AtnImage;
	Drf_t *drf;
	float fScaleFac;/* factor to scale the final projections by. This is
							used to implement the scaling of the primaries by
							fPrimaryFac. It is set to 1 if scatter is modeled as
							the SRF modeling takes care of the primary factor
							scaling itself. If no scatter modeling is being
							done, then it is set to fPrimaryFac*/
	
	if (psDrfTab != NULL) 
		drf=psDrfTab->psDrfs;
	else
		drf=NULL;
	fScaleFac = bModelSrf ? 1.0 : fPrimaryFac;

	if (!bActImgRotated){
		if (bFastRotate){
			vShearsRotate(0.5*M_PI+fAngle, reconmatrix, 
					sTempImages.RotReconImage, 0, 
					iNumSlices, iNumPixels);
		}else{
			rot3dpar(0.5*M_PI+fAngle, reconmatrix, 
					sTempImages.RotReconImage, 0, 
					iNumSlices, iNumPixels, iNumRotPixs);
		}
	}else{
	}
#ifdef DEBUG
	vPrintMsg(8,"sum after rotating=%.6g\n",
		sum_float(sTempImages.RotReconImage, iNumRotPixs*iNumRotPixs*iNumSlices));
#endif

	if(bModelAtn || bModelSrf){
		set_float(pfSumAtnTable,iNumSlices*iNumRotPixs,0.0);
		if (bAtnMapRotated){
			vPrintMsg(8, 
			"Project: skipping atten. map rotation for angle %d\n", iAngle);
		}else{
			vPrintMsg(8, 
				"Project: rotating atten. map for angle %d\n", iAngle);
			if (bFastRotate){
				vShearsRotate(0.5*M_PI+fAngle, atnmatrix, 
						pfRotAtnImage, 0, 
						iNumSlices, iNumPixels);
			}else{
				rot3dpar(0.5*M_PI+fAngle, atnmatrix, pfRotAtnImage, 0,
					iNumSlices, iNumPixels, iNumRotPixs);
			}
		}
	}
	if (bModelSrf){
		/* for effective scatter source modeling, we compute the 
				scatter source, then add it to the true source */
		CalcEffScatterSource(psSrf->pvScatterData, pfRotReconImage,
				pfRotAtnImage, iNumRotPixs, iNumRotPixs, iNumSlices,
				pfScatterEstimate);
		for(i=0; i<iNumSlices*iNumRotPixs*iNumRotPixs; ++i)
			pfScatterEstimate[i] = 
						pfRotReconImage[i]*fPrimaryFac +
						pfScatterEstimate[i];
			/* Scatter estimate now contains data to be projected*/
			pfRotPixels = pfScatterEstimate;
	}else{
		/* without scatter modeling, RotReconImage contains data 
			to be projected*/
		pfRotPixels = sTempImages.RotReconImage;
	}
	
	/* fStartDist is the distance from the center of rotation to the center
			of the first row of pixels. This is used to compute the distance
		for subsequent rows of pixels and hence the index into the DRF table
	*/
	if (bModelDrf){
		fStartDist =  fCFCR - fPixelWidth*iNumRotPixs/2 + 
									0.5*fPixelWidth;
		vPrintMsg(8,"    Start distance=%.2f, start index=%d\n",
			fStartDist, 
			(int)floor((fStartDist-psDrfTab->fStartDistance)/
								psDrfTab->fDrfSpacing + 0.5));
	}
	set_float(calcprj,iNumRotPixs*iNumSlices,0.0);
	for (iRow=0; iRow < iNumRotPixs; iRow++) {
		if (bModelDrf){
			/* fDist is distance from collimator to center of current row */	
			fDist=fStartDist + iRow*fPixelWidth; 
			/* iDist is the index into the Drf Table */
			iDist=(int)floor((fDist-psDrfTab->fStartDistance)/
							psDrfTab->fDrfSpacing + 0.5);
			iDist = MAX(0,iDist);
			iDist = MIN(iDist, psDrfTab->iNumDistances-1);
		}

		pfProjPixels = &(pfRotPixels[iRow*iNumSlices*iNumRotPixs]);
			/*pfProjPixels will contain a pointer to the pixels to be projected. It
			  starts out by pointing to the attenuation pixels, but this will 
			  point to the scatter+primary estimate and the convolved estimate as
			  different levels of modelling are included*/
			
		if(bModelAtn || bModelSrf)   /*sum till middle of pixel*/
			for (i=0; i < iNumSlices* iNumRotPixs; i++)
					pfSumAtnTable[i] += (0.5 * pfRotAtnImage[_23dptr(iRow, i)]);

		/* do the attenuation before convolving with the DRF. This improves
			accuracy for DRFs with long tails. Do into pfAtnImage. 
		*/
		if(bModelAtn){
			for (i=0; i < iNumSlices * iNumRotPixs; i++) 
				pfAtnImage[i] = pfProjPixels[i] *
									exp(-pfSumAtnTable[i]*fAtnFac);
			pfProjPixels=pfAtnImage;

		}
		/*if DRF is defined then use convolve. Note that the convolution
		 * will be done in spatial or frequency domain based on the value
		 * of psDrfTab->bFFTconvolve. the result will be put 
			into ConvImage and the projection mage pointer will be updated
			to point to this.
		 */
		if(bModelDrf){
			if (fDist > 0.0){/* for fdistance < don't do any convolution*/
				//fprintf(stderr,"Prj: ");
				ConvolveDrfYpad(pfProjPixels,iNumRotPixs, iNumSlices, 
					&drf[iDist],
					psDrfTab->bFFTconvolve,
					iAxialPadLength, iAxialAvgLength,
               sTempImages.ConvImage);
					pfProjPixels = sTempImages.ConvImage;
			}
		}/* note that if fDist is less than zero, pfProjPixels stays the same*/
		/* Finally, add the pixels for this plane to the calculated projections */
		for (i=0; i < iNumSlices * iNumRotPixs; i++)
			calcprj[i] += pfProjPixels[i]*fScaleFac;
				
		/*add the remaining half pixel's attenuation */
		if(bModelAtn || bModelSrf)
			for (i=0; i < iNumSlices * iNumRotPixs; i++)
				pfSumAtnTable[i] += 0.5 * pfRotAtnImage[_23dptr(iRow,i)];
	}
}

void BackProject(
	int iNumSlices,
	int iNumPixels,
	int iNumRotPixs,
	int iBckModel,
	float fPixelWidth,
	int iAngle,
	float fAngle,
	float fCFCR,
	DrfTab_t *psDrfTab, 
	Srf_t *psSrf, 
	float *atnmatrix,float* ratioprj, 
	float *imagematrix,
	float fPrimaryFac,
 	float fAtnFac,
	int iAxialPadLength, // Number of slices to pad at top and bottom to reduce
								// Edge effects (passed to ConvolveDrfYpad)
	int iAxialAvgLength, //Number of slices to average over when performing
								//axial padding
	int bAtnMapRotated,
	int bFastRotate)
{
	int iRow, iCol, iSlice, i, nRowSkip, nColSkip;
	float *pfBackProjPixels;
	float fStartDist;
	float fDist;
	int iDist;
	Drf_t *drf;
	float *pfSumAtnTable  = sTempImages.SumAtnTable;
	float *pfRotAtnImage = sTempImages.RotAtnImage;
	float *rotreconmatrix = sTempImages.RotReconImage;
	float *unrotreconmatrix = sTempImages.UnRotReconImage;
	float *pfAtnImage = sTempImages.AtnImage;
	/*float *imagematrix;*/ 
	int bModelAtn =iBckModel & MODEL_ATN;
	int bModelSrf = iBckModel & MODEL_SRF;
	int bModelDrf = iBckModel & MODEL_DRF;
	float fScaleFac;/* factor to scale the final projections by. This is
			 				used to implement the scaling of the primaries by
							fPrimaryFac. It is set to 1 if scatter is modeled as
							the SRF modeling takes care of the primary factor
							scaling itself. If no scatter modeling is being
							done, then it is set to fPrimaryFac*/

	/*pfBackProjPixels is the pointer that points to first the projection 
		(ratio) image and then updated as more levels of modeling are included 
		in the backprojector. 
		If attenuation or scatter are needed, rotate the attenuation matrix.
		If scatter is needed then use add_scat to add scatter to the projection 
		data and update the image pointer pfBackProjPixels to Scatter estimate,
		If drf needed, then for a particular row of the backprojection matrix 
		convolve the projection data with the psf of the drf for that distance
		Finally, put the projections into each row of the backprojection matrix*/

	if (psDrfTab != NULL && psDrfTab->psBckDrfTab != NULL)
		/* Use the Grf table for backprojection */
		psDrfTab = psDrfTab->psBckDrfTab;
	
	if (psDrfTab != NULL) 
		drf=psDrfTab->psDrfs;
	else
		drf=NULL;
	//fprintf(stderr,"pDrfTab=%d, drfs=%d\n",psDrfTab, drf);
	fScaleFac = bModelSrf ? 1.0 : fPrimaryFac;
	if(bModelAtn || bModelSrf){
		if (bAtnMapRotated){
			vPrintMsg(8, 
				"BackProject: skipping atten. map rotation for angle %d\n",
				iAngle);
		}else{
			vPrintMsg(8, 
				"BackProject: rotating atten. map for angle %d\n",
				iAngle);
			if (bFastRotate){
				vShearsRotate(0.5*M_PI+fAngle, atnmatrix, 
						pfRotAtnImage, 0, 
						iNumSlices, iNumPixels);
			}else{
				rot3dpar(0.5*M_PI+fAngle,atnmatrix, pfRotAtnImage, 0,
						iNumSlices, iNumPixels, iNumRotPixs);
			}
		}
		set_float(pfSumAtnTable,iNumSlices*iNumRotPixs,0.0);
	}

	/* fStartDist is the distance from the center of rotation to the center
			of the first row of pixels. This is used to compute the distance
		for subsequent rows of pixels and hence the index into the DRF table
	*/
	if (bModelDrf)
		fStartDist =  fCFCR - fPixelWidth*iNumRotPixs/2 + 
								0.5*fPixelWidth;
	for(iRow=0; iRow < iNumRotPixs; iRow++){
		/* pfBackProjPixels always points to what is to be added to the 
			backprojection image for the current plane
		*/
		pfBackProjPixels = ratioprj; 

		if (bModelDrf){
			/* fDist is distance from collimator to center of current row */	
			fDist=fStartDist + iRow*fPixelWidth; 
			/* iDist is the index into the Drf table */
			iDist=(int)floor((fDist - psDrfTab->fStartDistance)/
									psDrfTab->fDrfSpacing + 0.5);
			iDist = MAX(0,iDist);
			iDist = MIN(iDist, psDrfTab->iNumDistances-1);
		}

		if(bModelAtn || bModelSrf)
			for(i=0; i< iNumSlices*iNumRotPixs;i++)
				pfSumAtnTable[i] += (0.5 * pfRotAtnImage[_23dptr(iRow,i)]);


		/* first attenuation the pixels to be backprojected */
		if(bModelAtn){
			for (i=0; i<iNumSlices*iNumRotPixs; i++){ 
				pfAtnImage[i] = pfBackProjPixels[i]*exp(-pfSumAtnTable[i]*fAtnFac);
			}
			/* the values to be backprojected are now in AtnImage*/
			pfBackProjPixels = pfAtnImage;
		} 

		if(bModelDrf){
			if (fDist > 0.0){/* only do convolution for distances > 0*/
				//fprintf(stderr,"Bck: ");
				ConvolveDrfYpad(pfBackProjPixels,iNumRotPixs, iNumSlices, 
					&drf[iDist],
					psDrfTab->bFFTconvolve,
					iAxialPadLength, iAxialAvgLength,
               sTempImages.ConvImage);
				pfBackProjPixels = sTempImages.ConvImage;
			}
		}

		for (i=0; i< iNumSlices*iNumRotPixs; i++)
			rotreconmatrix[_23dptr(iRow, i)] = pfBackProjPixels[i]*fScaleFac;
		if(bModelAtn || bModelSrf)
			for (i=0; i<iNumSlices*iNumRotPixs; i++)
				pfSumAtnTable[i] += (0.5 * pfRotAtnImage[_23dptr(iRow, i)]);
	}

	if (bModelSrf){
		/* for effective scatter source modeling, we compute the 
				scatter source, then add it to the true source */
		vPrintMsg(8,"calculating effective source\n");
		CalcEffScatterSource(psSrf->pvScatterData, sTempImages.RotReconImage, 
			pfRotAtnImage,
			iNumRotPixs, iNumRotPixs, iNumSlices, sTempImages.ScatterEstimate );
		for(i=0; i<iNumSlices*iNumRotPixs*iNumRotPixs; ++i)
			sTempImages.RotReconImage[i] = 
										fPrimaryFac*sTempImages.RotReconImage[i]+
										sTempImages.ScatterEstimate[i];
	}

	if (bFastRotate){
		vShearsRotate(-0.5*M_PI-fAngle, rotreconmatrix, 
				unrotreconmatrix, 1, 
				iNumSlices, iNumPixels);
	}else{
		rot3dpar(-fAngle-0.5*M_PI, rotreconmatrix, 
				unrotreconmatrix,1, 
				iNumSlices, iNumPixels, iNumRotPixs);
	}

	nRowSkip = nColSkip= (iNumRotPixs - iNumPixels)/2;
	for(iRow=0;iRow<iNumPixels;iRow++)
		for(iSlice=0; iSlice<iNumSlices; iSlice++)
			for(iCol=0;iCol<iNumPixels;iCol++)
				imagematrix[iRow*iNumSlices*iNumPixels+iSlice*iNumPixels + iCol] += 
					unrotreconmatrix[_3dptr(iRow+nRowSkip, iSlice, iCol+nColSkip)];
	for(i=0; i < iNumPixels*iNumSlices*iNumPixels; i++)
		if(imagematrix[i] < EPS)
			imagematrix[i] = 0.0;
}
//Temporary matrices used in conv2d and conv2d_ypad
float *pfgConvMatrix = NULL; 
int igConvMatrixSize = 0;

void FreePrj()
{
	/* note that if these are NULL, as will be the case when an effect
		isn't modeled, IrlFree will do nothing*/
	IrlFree(sTempImages.AtnImage);
	sTempImages.AtnImage=NULL;
	IrlFree(sTempImages.RotReconImage);
	sTempImages.RotReconImage=NULL;
	IrlFree(sTempImages.UnRotReconImage);
	sTempImages.UnRotReconImage=NULL;
	IrlFree(sTempImages.RotAtnImage);
	sTempImages.RotAtnImage=NULL;
	IrlFree(sTempImages.SumAtnTable);
	sTempImages.SumAtnTable=NULL;
	IrlFree(sTempImages.ConvImage);
	sTempImages.ConvImage=NULL;
	IrlFree(sTempImages.ScatterEstimate);
	sTempImages.ScatterEstimate=NULL;
	IrlFree(pfgConvMatrix);
	pfgConvMatrix=NULL;
	vFreeFastRotate();
}
