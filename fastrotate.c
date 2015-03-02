/*
	fastrotate.c

	$Id: fastrotate.c 94 2010-01-22 01:33:19Z frey $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <math.h>

#include <mip/miputil.h>
#include <mip/printmsg.h>
#include <mip/errdefs.h>
#include "irlprivate.h"

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif
#ifdef WIN32
#define ceilf ceil
#define floorf floor
#endif
#define PTR(X,Y,Slice) ((Y)*NumSlices*NumPixels+(Slice)*NumPixels+(X))

/* File globals initialized by vInitFastRotate(iNumPixels) */

static int *piStart=NULL, *piEnd=NULL;

/***************************************************************************

	3-Pass Method of Shears fast rotation method:

	AW Paeth, "A fast algorithm for general raster rotation", Graphics
	Interface, pp. 77-81, 1986.

	A Tanaka, M Kameyama, S Kazama, and O Watanabe, "A rotation method
	for raster image using skew transformation", in Proc. IEEE Conf. Comp.
	Vision Patt. Recog., pp.272-277, 1986.

	EVR DiBella, AB Barclay, RL Eisner, and RW Schafer, "A comparison of
	rotation-based methods for iterative reconstruction algorithms", IEEE
	Trans. Nucl. Sci., 43(6), pp. 3370-3376, 1996.

	Rotation formula:

	Rotation is achieved through 3 steps:
		1) Shift by -r tan(theta/2) in x
		2) Shift by r sin(theta) in y
		3) Shift by -r tan(theta/2) in x again
	where r = perpendicular distance from COR to row(column) being shifted

	The rotation is only used for angles from -45 deg. to +45 deg.  For other
	angles, e.g. 60 deg., the image is first rotated by 90 eg. (simply done
	in place by swapping (i,j) indices for (NumPixels-j-1,i), then rotated to the
	final angle, e.g. -30 deg. in our example.

	Note that this version rotates pixels that are in X, Z, Y order. To speed
	things up it copies them slice-by-slice into an array in X,Y order, and
	rotates that slice


*************************************************************************/

void vShearsRotateCopy (float fAngle, float *pfOldMatrix, float *pfRotMatrix,
		int PB_FLAG, int NumSlices, int NumPixels)
{
	int i, j, jj, jk, k, iTmp, jTmp;
	int mDim, iAdd;
	int iStart, iEnd, jStart, jEnd;
	int b90, b180;
	float *pfShiftedX, *pfShiftedY;
	float fCntr, tCntr, fDelta;
	float yDist, xDist, f1, f2;
	static int *piStart=NULL, *piEnd=NULL;
#ifdef DEBUG0
	char pchSaveImageName[256];
#endif

	/*	PrintTimes("Start vRotate"); */
	
	if (piStart == NULL || piEnd == NULL){
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, 
			"vShearsRotateCopy", "piStart || piEnd have not been initialized");
	}
	
	set_float(pfRotMatrix,NumPixels*NumPixels*NumSlices,0.0);

#ifdef DEBUG0
      sprintf(pchSaveImageName,"rot0.%f.im", fAngle);
      writeimage(pchSaveImageName,NumPixels,NumPixels,1,&(pfOldMatrix[NumPixels*NumPixels*33]));
#endif

/*** Using the 3-pass method of shears, only want to rotate within
	the range of -45 deg to +45 deg ***/
	b90 = b180 = 0;
	while (fAngle > 0.75*M_PI) {
	 fAngle -= M_PI; /* Must flip-flop i,j to NumPixels-i-1,NumPixels-j-1 */
	 b180++;
	 if (b180 == 2) b180=0;
	}
	while (fAngle <= -M_PI/4.0) {
	 fAngle += M_PI;
	 b180++;
	 if (b180 == 2) b180=0;
	}
	if (fAngle > M_PI/4.0) {
	 /* Rotate image by 90 deg. 1st (flip-flop i,j to NumPixels-j-1,i), then rotate */
	 b90 = 1;
	 fAngle -= M_PI/2.0; /* Must rotate image by 90 deg. for this to work */
	}

/*** Find out how big we need to make the intermediate matrices ***/
	fCntr = ((float)NumPixels-1.0)*0.5;
	iAdd = 0;
	for (j=NumPixels/4; j<NumPixels/2; j++) {
	 yDist = (float)j - fCntr;
	 fDelta = fabs(yDist * tan((double)(fAngle*0.5)));
	 fDelta -= fCntr - sqrt(fCntr*fCntr - yDist*yDist);
	 iTmp = ceilf(fDelta);
	 if (iTmp > iAdd) iAdd = iTmp;
	}
	mDim = NumPixels + 2*iAdd;
	tCntr = ((float)mDim-1.0)*0.5;

/*** Allocate space for temporary matrices ***/
	pfShiftedX = (float *)malloc(mDim*mDim*sizeof(float));
	set_float(pfShiftedX,mDim*mDim,0.0);
	pfShiftedY = (float *)malloc(mDim*mDim*sizeof(float));
	set_float(pfShiftedY,mDim*mDim,0.0);

/*** Work on 1 slice at a time to reduce intermediate memory requirements ***/
	for (k=0; k<NumSlices; k++) {

/*** Perform first x shift ***/
	 if (b180 == 0) {
	  if (b90 == 0) {

		for (j=0; j<NumPixels; j++) {
		 yDist = (float)j - fCntr;
		 fDelta = -1.0 * yDist * tan((double)(fAngle*0.5));
		 iTmp = floor((double)fDelta);
		 f2 = fDelta - (float)iTmp;
		 f1 = 1.0 - f2;
		 iStart = piStart[j] - iTmp;
		 iEnd = piEnd[j] - iTmp - 1;

		 iTmp -= iAdd;
		 jk = j + iAdd;
		 for (i=iStart+iAdd; i<iEnd+iAdd; i++) {
		  pfShiftedX[i+mDim*jk] = f1 * pfOldMatrix[PTR(i+iTmp,j,k)] +
			f2 * pfOldMatrix[PTR(i+iTmp+1,j,k)];
		 }
	   }

	  } else { /* b90 = 1; flip-flop by 90 deg. */

		for (j=0; j<NumPixels; j++) {
		 jj = NumPixels - j - 1;
		 yDist = (float)j - fCntr;
		 fDelta = -1.0 * yDist * tan((double)(fAngle*0.5));
		 iTmp = floorf(fDelta);
		 f2 = fDelta - (float)iTmp;
		 f1 = 1.0 - f2;
		 iStart = piStart[j] - iTmp;
		 iEnd = piEnd[j] - iTmp - 1;

		 iTmp -= iAdd;
		 jk = j + iAdd;
		 for (i=iStart+iAdd; i<iEnd+iAdd; i++) {
		  pfShiftedX[i+mDim*jk] = f1 * pfOldMatrix[PTR(jj,i+iTmp,k)] +
			f2 * pfOldMatrix[PTR(jj,i+iTmp+1,k)];
		 }
		}

	  } /* end if b90 */

	 } else { /* b180 == 1 */

	  if (b90 == 0) {

		for (j=0; j<NumPixels; j++) {
		 jj = NumPixels - j - 1;
		 yDist = (float)j - fCntr;
		 fDelta = -1.0 * yDist * tan((double)(fAngle*0.5));
		 iTmp = floorf(fDelta);
		 f2 = fDelta - (float)iTmp;
		 f1 = 1.0 - f2;
		 iStart = piStart[j] - iTmp;
		 iEnd = piEnd[j] - iTmp - 1;

		 iTmp -= iAdd;
		 jk = j + iAdd;
		 iTmp += 1;
		 for (i=iStart+iAdd; i<iEnd+iAdd; i++) {
		  pfShiftedX[i+mDim*jk] = f1 * pfOldMatrix[PTR(NumPixels-i-iTmp,jj,k)] +
			f2 * pfOldMatrix[PTR(NumPixels-i-iTmp-1, jj, k)];
		 }
		}

	  } else { /* b90 = 1; flip-flop by 90 deg. */

		for (j=0; j<NumPixels; j++) {
		 yDist = (float)j - fCntr;
		 fDelta = -1.0 * yDist * tan((double)(fAngle*0.5));
		 iTmp = floorf(fDelta);
		 f2 = fDelta - (float)iTmp;
		 f1 = 1.0 - f2;
		 iStart = piStart[j] - iTmp;
		 iEnd = piEnd[j] - iTmp - 1;

		 iTmp -= iAdd;
		 jk = j + iAdd;
		 iTmp += 1;
		 for (i=iStart+iAdd; i<iEnd+iAdd; i++) {
		  pfShiftedX[i+mDim*jk] = f1 * pfOldMatrix[PTR(j,NumPixels-i-iTmp,k)]+
			f2 * pfOldMatrix[PTR(j,NumPixels-i-iTmp-1,k)];
		 }
		}

	  } /* end if b90 */

	 } /* end if b180 */

#ifdef DEBUG0
	if (k == 33) {
      sprintf(pchSaveImageName,"rot1.%f.im", fAngle);
      writeimage(pchSaveImageName,mDim,mDim,1,pfShiftedX);
	}
#endif

/*** Perform y shift ***/
	 for (i=0; i<mDim; i++) {

	  xDist = (float)i - tCntr;
	  fDelta = xDist * sin((double)(fAngle));
	  jTmp = floorf(fDelta);
	  f2 = fDelta - (float)jTmp;
	  f1 = 1.0 - f2;

	  jStart = 0 - jTmp;
	  if (jStart < 0) jStart = 0;
	  jEnd = mDim - jTmp - 1;
	  if (jEnd > mDim) jEnd = mDim;

	  for (j=jStart; j<jEnd; j++) {
		pfShiftedY[i+mDim*j] = f1 * pfShiftedX[i+mDim*(j+jTmp)] +
			f2 * pfShiftedX[i+mDim*(j+jTmp+1)];
	  }

	 }

#ifdef DEBUG0
	if (k == 33) {
      sprintf(pchSaveImageName,"rot2.%f.im", fAngle);
      writeimage(pchSaveImageName,mDim,mDim,1,pfShiftedY);
	}
#endif

/*** Perform second x shift ***/
	 for (j=0; j<NumPixels; j++) {

	  yDist = (float)j - fCntr;
	  fDelta = -1.0 * yDist * tan((double)(fAngle*0.5));
	  iTmp = floorf(fDelta);
	  f2 = fDelta - (float)iTmp;
	  f1 = 1.0 - f2;

	  /* Avoid edge effects, make all backprojected images (e.g. normalization
		matrices and update matrices) 1 pixel smaller on each side; this avoids
		edge-truncation errors encountered during the linear interpolations */
	  iStart = piStart[j]+1;
	  iEnd = piEnd[j]-1;

	  iTmp += iAdd;
	  jj = j + iAdd;
	  for (i=iStart; i<iEnd; i++) {
		pfRotMatrix[PTR(i,j,k)] = f1 * pfShiftedY[i+iTmp+mDim*jj] +
			f2 * pfShiftedY[i+iTmp+1+mDim*jj];
	  }

	 }

/*** Done with this slice, continue loop over all slices ***/
	}

#ifdef DEBUG0
      sprintf(pchSaveImageName,"rot3.%f.im", fAngle);
      writeimage(pchSaveImageName,NumPixels,NumPixels,1,&(pfRotMatrix[NumPixels*NumPixels*33]));
#endif

	free(pfShiftedX);
	free(pfShiftedY);

/*	PrintTimes("End vRotate"); */

	return;
}

/***************************************************************************

	3-Pass Method of Shears fast rotation method:

	AW Paeth, "A fast algorithm for general raster rotation", Graphics
	Interface, pp. 77-81, 1986.

	A Tanaka, M Kameyama, S Kazama, and O Watanabe, "A rotation method
	for raster image using skew transformation", in Proc. IEEE Conf. Comp.
	Vision Patt. Recog., pp.272-277, 1986.

	EVR DiBella, AB Barclay, RL Eisner, and RW Schafer, "A comparison of
	rotation-based methods for iterative reconstruction algorithms", IEEE
	Trans. Nucl. Sci., 43(6), pp. 3370-3376, 1996.

	Rotation formula:

	Rotation is achieved through 3 steps:
		1) Shift by -r tan(theta/2) in x
		2) Shift by r sin(theta) in y
		3) Shift by -r tan(theta/2) in x again
	where r = perpendicular distance from COR to row(column) being shifted

	The rotation is only used for angles from -45 deg. to +45 deg.  For other
	angles, e.g. 60 deg., the image is first rotated by 90 eg. (simply done
	in place by swapping (i,j) indices for (NumPixels-j-1,i), then rotated to the
	final angle, e.g. -30 deg. in our example.

*************************************************************************/

void vShearsRotate (float fAngle, float *pfOldMatrix, float *pfRotMatrix,
		int PB_FLAG, int NumSlices, int NumPixels)
{
	int i, j, jj, jk, k, iTmp, jTmp;
	int mDim, iAdd;
	int iStart, iEnd, jStart, jEnd;
	int b90, b180;
	float *pfShiftedX, *pfShiftedY;
	float fCntr, tCntr, fDelta;
	float yDist, xDist, f1, f2;
#ifdef DEBUG0
	char pchSaveImageName[256];
#endif
	
	/*	PrintTimes("Start vRotate"); */
	
	if (piStart == NULL || piEnd == NULL){
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, 
			"vShearsRotate", "piStart || piEnd have not been initialized");
	}
	
	set_float(pfRotMatrix,NumPixels*NumPixels*NumSlices,0.0);

#ifdef DEBUG0
      sprintf(pchSaveImageName,"rot0.%f.im", fAngle);
      writeimage(pchSaveImageName,NumPixels,NumPixels,1,&(pfOldMatrix[NumPixels*NumPixels*33]));
#endif

/*** Using the 3-pass method of shears, only want to rotate within
	the range of -45 deg to +45 deg ***/
	b90 = b180 = 0;
	while (fAngle > 0.75*M_PI) {
	 fAngle -= M_PI; /* Must flip-flop i,j to NumPixels-i-1,NumPixels-j-1 */
	 b180++;
	 if (b180 == 2) b180=0;
	}
	while (fAngle <= -M_PI/4.0) {
	 fAngle += M_PI;
	 b180++;
	 if (b180 == 2) b180=0;
	}
	if (fAngle > M_PI/4.0) {
	 /* Rotate image by 90 deg. 1st (flip-flop i,j to NumPixels-j-1,i), then rotate */
	 b90 = 1;
	 fAngle -= M_PI/2.0; /* Must rotate image by 90 deg. for this to work */
	}

/*** Find out how big we need to make the intermediate matrices ***/
	fCntr = ((float)NumPixels-1.0)*0.5;
	iAdd = 0;
	for (j=NumPixels/4; j<NumPixels/2; j++) {
	 yDist = (float)j - fCntr;
	 fDelta = fabs(yDist * tan((double)(fAngle*0.5)));
	 fDelta -= fCntr - sqrt(fCntr*fCntr - yDist*yDist);
	 iTmp = ceil(fDelta);
	 if (iTmp > iAdd) iAdd = iTmp;
	}
	mDim = NumPixels + 2*iAdd;
	tCntr = ((float)mDim-1.0)*0.5;

/*** Allocate space for temporary matrices ***/
	pfShiftedX = (float *)malloc(mDim*mDim*sizeof(float));
	set_float(pfShiftedX,mDim*mDim,0.0);
	pfShiftedY = (float *)malloc(mDim*mDim*sizeof(float));
	set_float(pfShiftedY,mDim*mDim,0.0);

/*** Work on 1 slice at a time to reduce intermediate memory requirements ***/
	for (k=0; k<NumSlices; k++) {

/*** Perform first x shift ***/
	 if (b180 == 0) {
	  if (b90 == 0) {

		for (j=0; j<NumPixels; j++) {
		 yDist = (float)j - fCntr;
		 fDelta = -1.0 * yDist * tan((double)(fAngle*0.5));
		 iTmp = floorf(fDelta);
		 f2 = fDelta - (float)iTmp;
		 f1 = 1.0 - f2;
		 iStart = piStart[j] - iTmp;
		 iEnd = piEnd[j] - iTmp - 1;

		 iTmp -= iAdd;
		 jk = j + iAdd;
		 for (i=iStart+iAdd; i<iEnd+iAdd; i++) {
		  //pfShiftedX[i+mDim*jk] = f1 * pfOldMatrix[i+iTmp+NumPixels*(j+NumPixels*k)] +
		  pfShiftedX[i+mDim*jk] = f1 * pfOldMatrix[PTR(i+iTmp,j,k)] +
			f2 * pfOldMatrix[PTR(i+iTmp+1,j,k)];
		 }
	   }

	  } else { /* b90 = 1; flip-flop by 90 deg. */

		for (j=0; j<NumPixels; j++) {
		 jj = NumPixels - j - 1;
		 yDist = (float)j - fCntr;
		 fDelta = -1.0 * yDist * tan((double)(fAngle*0.5));
		 iTmp = floorf(fDelta);
		 f2 = fDelta - (float)iTmp;
		 f1 = 1.0 - f2;
		 iStart = piStart[j] - iTmp;
		 iEnd = piEnd[j] - iTmp - 1;

		 iTmp -= iAdd;
		 jk = j + iAdd;
		 for (i=iStart+iAdd; i<iEnd+iAdd; i++) {
		  pfShiftedX[i+mDim*jk] = f1 * pfOldMatrix[PTR(jj,i+iTmp,k)] +
			f2 * pfOldMatrix[PTR(jj,i+iTmp+1,k)];
		 }
		}

	  } /* end if b90 */

	 } else { /* b180 == 1 */

	  if (b90 == 0) {

		for (j=0; j<NumPixels; j++) {
		 jj = NumPixels - j - 1;
		 yDist = (float)j - fCntr;
		 fDelta = -1.0 * yDist * tan((double)(fAngle*0.5));
		 iTmp = floorf(fDelta);
		 f2 = fDelta - (float)iTmp;
		 f1 = 1.0 - f2;
		 iStart = piStart[j] - iTmp;
		 iEnd = piEnd[j] - iTmp - 1;

		 iTmp -= iAdd;
		 jk = j + iAdd;
		 iTmp += 1;
		 for (i=iStart+iAdd; i<iEnd+iAdd; i++) {
		  pfShiftedX[i+mDim*jk] = f1 * pfOldMatrix[PTR(NumPixels-i-iTmp,jj,k)] +
			f2 * pfOldMatrix[PTR(NumPixels-i-iTmp-1, jj, k)];
		 }
		}

	  } else { /* b90 = 1; flip-flop by 90 deg. */

		for (j=0; j<NumPixels; j++) {
		 yDist = (float)j - fCntr;
		 fDelta = -1.0 * yDist * tan((double)(fAngle*0.5));
		 iTmp = floorf(fDelta);
		 f2 = fDelta - (float)iTmp;
		 f1 = 1.0 - f2;
		 iStart = piStart[j] - iTmp;
		 iEnd = piEnd[j] - iTmp - 1;

		 iTmp -= iAdd;
		 jk = j + iAdd;
		 iTmp += 1;
		 for (i=iStart+iAdd; i<iEnd+iAdd; i++) {
		  pfShiftedX[i+mDim*jk] = f1 * pfOldMatrix[PTR(j,NumPixels-i-iTmp,k)]+
			f2 * pfOldMatrix[PTR(j,NumPixels-i-iTmp-1,k)];
		 }
		}

	  } /* end if b90 */

	 } /* end if b180 */

#ifdef DEBUG0
	if (k == 33) {
      sprintf(pchSaveImageName,"rot1.%f.im", fAngle);
      writeimage(pchSaveImageName,mDim,mDim,1,pfShiftedX);
	}
#endif

/*** Perform y shift ***/
	 for (i=0; i<mDim; i++) {

	  xDist = (float)i - tCntr;
	  fDelta = xDist * sin((double)(fAngle));
	  jTmp = floorf(fDelta);
	  f2 = fDelta - (float)jTmp;
	  f1 = 1.0 - f2;

	  jStart = 0 - jTmp;
	  if (jStart < 0) jStart = 0;
	  jEnd = mDim - jTmp - 1;
	  if (jEnd > mDim) jEnd = mDim;

	  for (j=jStart; j<jEnd; j++) {
		pfShiftedY[i+mDim*j] = f1 * pfShiftedX[i+mDim*(j+jTmp)] +
			f2 * pfShiftedX[i+mDim*(j+jTmp+1)];
	  }

	 }

#ifdef DEBUG0
	if (k == 33) {
      sprintf(pchSaveImageName,"rot2.%f.im", fAngle);
      writeimage(pchSaveImageName,mDim,mDim,1,pfShiftedY);
	}
#endif

/*** Perform second x shift ***/
	 for (j=0; j<NumPixels; j++) {

	  yDist = (float)j - fCntr;
	  fDelta = -1.0 * yDist * tan((double)(fAngle*0.5));
	  iTmp = floorf(fDelta);
	  f2 = fDelta - (float)iTmp;
	  f1 = 1.0 - f2;

	  /* Avoid edge effects, make all backprojected images (e.g. normalization
		matrices and update matrices) 1 pixel smaller on each side; this avoids
		edge-truncation errors encountered during the linear interpolations */
	  iStart = piStart[j]+1;
	  iEnd = piEnd[j]-1;

	  iTmp += iAdd;
	  jj = j + iAdd;
	  for (i=iStart; i<iEnd; i++) {
		pfRotMatrix[PTR(i,j,k)] = f1 * pfShiftedY[i+iTmp+mDim*jj] +
			f2 * pfShiftedY[i+iTmp+1+mDim*jj];
	  }

	 }

/*** Done with this slice, continue loop over all slices ***/
	}

#ifdef DEBUG0
      sprintf(pchSaveImageName,"rot3.%f.im", fAngle);
      writeimage(pchSaveImageName,NumPixels,NumPixels,1,&(pfRotMatrix[NumPixels*NumPixels*33]));
#endif

	free(pfShiftedX);
	free(pfShiftedY);

/*	PrintTimes("End vRotate"); */

	return;
}

#undef PTR

/*
	Allocate the ivectors for the static globals piStart and piEnd.
	This is only done once per reconstruction. This function is called from
	InitPrj() in prjbck.c
*/

void vInitFastRotate(int NumPixels)
{
	int j, iLen;
	float fTmp, fCntr;

	vPrintMsg(2," Setting up in-plane circular reconstruction space\n");
	piStart = ivector(NumPixels, "fastrotate: piStart");
	piEnd = ivector(NumPixels, "fastrotate: piEnd");
	fCntr = (NumPixels-1)*0.5;
	for (j=0; j<NumPixels; j++) {
		fTmp = ((float)j - fCntr);
		iLen= 2*floor(sqrt(fCntr*fCntr-fTmp*fTmp))+2; /* 2 <= iLen <= nPixels */
		piStart[j] = (NumPixels - iLen)/2;
		piEnd[j] = piStart[j] + iLen;
	}
	vPrintMsg(2," Done setup\n");
}

void vFreeFastRotate()
{
	if (piStart != NULL){
		free_vector(piStart);
		piStart=NULL;
	}
	if (piEnd != NULL){
		free_vector(piEnd);
		piStart=NULL;
	}
}
