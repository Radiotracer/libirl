/*
	initial_est.c
	
	$Id: initial_est.c 79 2007-02-27 04:28:57Z mjs $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <stdlib.h>

#include <mip/miputil.h>
#include <math.h>
#include <mip/errdefs.h>
#include <mip/printmsg.h>

#ifdef WIN32
#define hypot _hypot
#endif

#include "irlprivate.h"

/* purpose is to compute initial estimate and put it in the portion of the
   image that can contain the object as limted by the camera orbit*/

static void vZeroBehind(int iNdim, float fSin, float fCos, 
		float fRadius, float *pfPix);

float *pfComputeInitialEstimate(
		Parms_t *psParms, 
		int bUseContourSupport,
		float fAtnMapThresh,
		View_t *Views, 
		float *pfPrjImage, 
		float *pfAtnMap,
		float *pfInitImage)
	/* computes the initial estimate by putting the average into either
		the entire image or the region defined by the orbit (everything behind
		plane of face of camera is zeroed for each orbit. Also, the
		region where the attenuation values are below a user-specified threshold 
		can be zeroed. A pointer to the initial image is returned. Note that
		the output image is in "reordered" format (i.e w/y direction varying 
		slowest)
		
	*/
{
	float fPixelRadRot;
	float fAvg;
	int iAng, iSlice, iCol, iRow, iPix, nPixels, nSlices,nAngs;
	int UseContour;
	int bUseMapSupport=FALSE;
	float fMapThresh;
	float *pfWriteMatrix=NULL;
	float *pf2dImage=NULL;
	float fNumReconPix;

	vPrintMsg(5,"ComputeInitialEstimate\n");
	nPixels = psParms->NumPixels;
	nSlices=psParms->NumSlices;
	nAngs = psParms->NumAngles;
	if (pfInitImage == NULL)
		pfInitImage = vector(nPixels*nSlices*nPixels,
						"ComputeInitialEstimate:InitImage");
	fNumReconPix = (float)(nPixels*nPixels*nSlices);
	/* for a circular reconstruction, only a cylindrical region is used*/
	if (psParms->Circular)
		fNumReconPix *= M_PI/4.0;
		
	fAvg = sum_float(pfPrjImage, psParms->NumRotPixs*nSlices*nAngs)/
				(fNumReconPix*nAngs);
	vPrintMsg(6,"average pixel value in reconstructed initial image=%.3g\n",
			fAvg);
	UseContour = bUseContourSupport;
	vPrintMsg(6, "use contour=%d, atnmapthresh=%d (%.3f)\n", 
		UseContour, pfAtnMap != NULL, fAtnMapThresh);
	if (pfAtnMap != NULL){
		fMapThresh=fAtnMapThresh;
		fMapThresh *= psParms->PixelWidth; /* convert to per pixel */
		if (fMapThresh < 0.0){
			vErrorHandler(ECLASS_WARN, ETYPE_ILLEGAL_VALUE, 
			 	"ComputeInitialEstimate", 
				"atnmap_support_thresh must be > 0, not %.3f. Ignoring",
				fMapThresh);
		}else{
			bUseMapSupport = fMapThresh > 0.0;
		}
	}

	if (!UseContour){
		/* use average value in entire image for each isotope */
		set_float(pfInitImage, nPixels*nSlices*nPixels, fAvg);
	}else{
		/* put average value only in pixels that are inside the contour
			defined by the detector face */
		/* first build a 2d image that excludes region behind the detector*/
		pf2dImage = vector(nPixels*nPixels,"ComputeInitialEstimate: 2dImage");
		/* set it all to average value*/
		set_float(pf2dImage, nPixels*nPixels,fAvg);

		for (iAng=0; iAng<nAngs; iAng++){
			/* loop over angles and set everything behind detector to zero */
			fPixelRadRot = Views[iAng].CFCR / psParms->PixelWidth;
			vZeroBehind(nPixels,Views[iAng].Sin, Views[iAng].Cos, fPixelRadRot,
				pf2dImage);
		}
		/* now copy the pixels into the reconstructed image pixels */
		for (iSlice=0; iSlice<nSlices; iSlice++){
			for(iRow=0; iRow < nPixels; ++iRow){
				for (iCol=0; iCol<nPixels; iCol++){
				  /* note that the InitImage is in "reordered format" where
				  	  pixels are stored with x (iCol), then z (iSlice), they y 
					  (iRow) varying from fastest to slowest*/
				  pfInitImage[iRow*nSlices*nPixels + iSlice*nPixels + iCol] = 
				  				pf2dImage[iCol+iRow*nPixels];
				}
			}
		}
	}
	if (bUseMapSupport){
		/* with UseMapSupport we set all pixels equal to zero for which
			the AtnMap falls below fMapThresh */
		/* note that the Atn map must already be in "reordered" format */
		for(iPix=0; iPix < nSlices*nPixels*nPixels; ++iPix){
		  if (pfAtnMap[iPix] < fMapThresh){
			  pfInitImage[iPix] = 0.0;	
			}
		}
	}

	if (psParms->Circular){
		/* for circular reconstructions, zero all pixels outside the circle*/
		/* do this by brute force: check each pixel to see if it has
			a distance of less than fRadius from the center of the recon image*/
		/* we make sure that the entire pixel is in the image */
		float fRadius = (float)nPixels*0.5;
		float fX, fY;
		float fCent;

		fCent=fRadius;
		for(iRow=0; iRow < nPixels; ++iRow){
			fY = (float)iRow - fCent;
			if (fY > 0) fY += 1;
			for (iCol=0; iCol<nPixels; iCol++){
				fX = (float)iCol - fCent;
				if (fX > 0) fX += 1;
				if (hypot(fX, fY) > fRadius ){
					for (iSlice=0; iSlice<nSlices; iSlice++){
					  pfInitImage[iRow*nSlices*nPixels+iSlice*nPixels+iCol] = 0.0;
					}
				}
			}
		}
	}

	if (pfWriteMatrix != NULL)free_vector(pfWriteMatrix);
	if (pf2dImage != NULL)free_vector(pf2dImage);
	return(pfInitImage);
}

static void vZeroBehind(int iNdim, float fSin, float fCos, 
		float fRadius, float *pfPix)
{
	/* this zeroes all pixels behind the plane located at the angle defined
		by fSin and fCos and fRadius. It works in a brute force way by taking
		the dot product of the position vector with the normal to the plane
		(Sin, Cos) and checking if this is greater than fRadius */
	int ix, iy;
	float fX, fY;
	float fR;

	for(iy=0; iy<iNdim; ++iy){
		fY = iy-(float)iNdim/2.0;
		for(ix=0; ix<iNdim; ++ix){
			fX = ix-(float)iNdim/2.0;
			fR = floor(fX*fCos + fY*fSin + 0.5);
			if (fR > fRadius)
				pfPix[ix+iy*iNdim] = 0.0;
		}
	}
}
