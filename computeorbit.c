/*
	computeorbit.c

	$Id: computeorbit.c 125 2011-09-26 15:02:12Z frey $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mip/miputil.h>

#include "irlprivate.h"

#define ERR_MALLOC 1
#define ERR_ORBIT_TYPE 2

int ComputeOrbit(float *pfPixels, //2D image used to compute orbit
	IrlParms_t *psParms,
	float fMargin, //Gap between edge of patient and camera face
						// that is used in computing the CFCR for each view
						//For circular orbits, this is treated as the CFCR
	float *pfMargins,//If non-null, fMargin is ignored and this must be a
						// vector of length NumViews, each of which contains
						// a margin around the patient that is added to the distance
						 // of closest approach  
	float fThresh, 	//Threshold specifying used in computing edge of object
						//If positive, then it is the percent of maximum
						//If negative, then it is an absolute threshold

	int OrbitType, //0: circular orbit. fMargin is treated as the CFCR for
						//   each view is placed in the CFCR element for each view
						//1: noncircular orbit with camera at each view moved 
						//		independently so it is fMargin from the edge of the
						//		patient as defined by the pixel values in fImage
						//2: Noncircular orbit with camera in "L-mode"
						//   configuration. It is assumed that the camera moves
						//   in and out as a unit with the vertex of the camera
						//   bisected by a line through the center of rotation.
						//   It is also assumed that camera view  i and i+NumViews/2
						//   are from the same camera configuration.
	PrjView_t *psViews//Array of view structures. The CFCR is computed
						 //for each view and put in the CFCR element for the
						 // view. These must already have the projection angle
						 // placed in the view structure
	)
{
	int iNdim;
	int iNumRotPix;
	float fCos, fSin;
	float fAngle;
	float fPixSize;
	int iNang;
	int iPix;
	float *pfRotPix;
	float fRadius;
	float fSum;
	float fPixThresh=0.0;
	float fMax;
	int iAng;
	int iRow;

	iNang = psParms->NumViews;
	iNdim = psParms->NumPixels;
	fPixSize = psParms->BinWidth;
	
	if (OrbitType == 0){
		//Circular orbit
		for(iAng=0; iAng < iNang; ++iAng){
			psViews[iAng].CFCR=fMargin;
		}
		return 0;
	}else if (OrbitType > 2){
		//only orbit types 0, 1, and 2 are supported
		return ERR_ORBIT_TYPE;
	}

	/* Threshold the image, getting rid of small pixel values */
	if (fThresh <= 0.0)
		fPixThresh=-fThresh;
	else{
		fMax = pfPixels[0];
		for(iPix=1; iPix < iNdim*iNdim ; ++iPix)
			fMax = pfPixels[iPix] > fMax ? pfPixels[iPix] : fMax;
		fPixThresh=fMax * fThresh;
	}
			
	for(iPix=0; iPix < iNdim*iNdim ; ++iPix)
		if (pfPixels[iPix] <= fPixThresh) pfPixels[iPix] = 0.0;

	//Shears Rotate only supports input and output matrix the same size
	//iNumRotPix=iNdim;
	iNumRotPix = (int)ceil(sqrt(2.0)*(double)iNdim);
	if ((iNumRotPix % 2) != (iNdim % 2))iNumRotPix++;
	pfRotPix = (float *)malloc(iNumRotPix*iNumRotPix*sizeof(float));
	if (pfRotPix == NULL) return (ERR_MALLOC);

	for(iAng=0; iAng < iNang; ++iAng){
		fAngle=psViews[iAng].Angle;
		fCos=(float)cos(fAngle);
		fSin=(float)sin(fAngle);
		//Rotate the image so the rows are parallel to the detector.
		//Row 0 will be the furthest row from the center of rotation in the
		// direction of the detector.
		//Using rot3dpar and a larger rotation matrix handles the case where
		// the object is off-center
		//vShearsRotate(0.5*M_PI + fAngle, pfPixels, pfRotPix, 0, 1, iNdim);
		rot3dpar((float)(0.5*M_PI + fAngle), pfPixels, pfRotPix, 0, 1, iNdim,
			iNumRotPix);

		/* check for nonzero columns, which are now stored with the index
			varying fastest and iRow=0 closest to the detector. Do this
			by summing the pixels in a Row. If they are nonzero stop*/
		// since all pixels less than the threshold are set to zero, all we
		// need to do to check if there are pixels above the threshold in a row
		// is to make sure that the sum is > 0. Step through the rows starting
		// at the detector.
		fSum = 0.0;
		for(iRow = 0; iRow < iNumRotPix/2 && fSum <= 0.0;	 ++iRow){
			fSum = sum_float(pfRotPix+iRow*iNumRotPix, iNumRotPix);
		}
		fRadius = ((float)iNumRotPix/2.0 - (float)iRow)*fPixSize;
		if (pfMargins != NULL)
			fRadius += pfMargins[iAng];
		else
			fRadius += fMargin;
		psViews[iAng].CFCR=fRadius;
	}
	if (OrbitType == 2){
		// Noncircular orbit in L-mode with vertex of camera bisected by line 
		// connecting center of rotation and vertex point. It is assumed that
		// Projection i and i+NumViews/2 are from the same detector position.
		// In this case the two detectors don't move independently. So, we use
		// the radius for one detector and use it to compute the radius for the
		// other detector. If this computed radius is greater than the one
		// allowing the closest approach computed above, 
		for(iAng=0; iAng < iNang/2; ++iAng){
			if (psViews[iAng].CFCR > psViews[iAng+iNang/2].CFCR){
				psViews[iAng+iNang/2].CFCR = psViews[iAng].CFCR;
			}else{
				psViews[iAng].CFCR = psViews[iAng+iNang/2].CFCR;
			}
		}
	}
	free(pfRotPix);
	return 0;
}
