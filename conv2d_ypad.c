/*
	conv2d_ypad.c 

	$Id: conv2d_ypad.c 83 2007-05-25 18:20:44Z binhe $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <mip/miputil.h>

#ifdef DEBUG
#include <mip/imgio.h>
#endif

#include "irlprivate.h"

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

/*
	convolves the array a with b. The length of a is given by aSizeY * aSizeX.
	the array b has an odd*odd number of elements and the center is
	at position b[bCtrY][bCtrX]. Thus the total length of the blurring function
	is (2*bCtrY+1) * (2*bCtrX+1). The length of the output array is 
	aSizeY * aSizeX.

	Before convolving, the array is copied into an array that is padded on both
	ends by iYpad rows. For each end, the padded portion is filled with the
	average of the iYavg rows at the end of the image
*/

extern float *pfgConvMatrix;
extern int igConvMatrixSize;

void conv2d_ypad(float *a, int aSizeX, int aSizeY, 
	float *b, int bCtrX, int bCtrY,
	int iPadY, int iAvgY,
	float *out)
{
#define a_ptr(aY, aX) ((aY) * aSizeX + (aX))
#define b_ptr(bY, bX) ((bY) * (2*bCtrX+1) + (bX))
#define t_ptr(tY, tX) ((tY) * (aSizeX+2*bCtrX) + (tX))

	int	aY, aX, bY, bX;
	int iNewXdim, iNewYdim;
	float *tmp;

	/* for other size b's, place the array a in tmp with bCtrX zeros on 
		 left and right to avoid need to test for array bounds*/
	iPadY = MAX(iPadY,0);
	iPadY = MIN(iPadY,bCtrY);
	iAvgY = MAX(iAvgY,0);
	iAvgY = MIN(iAvgY, aSizeY);
	iNewXdim = aSizeX + 2*bCtrX;
	iNewYdim = aSizeY + 2*iPadY;

	if (pfgConvMatrix  == NULL || 
		(igConvMatrixSize < iNewXdim*iNewYdim)){
		pfgConvMatrix = (float *)
			pvIrlRealloc(pfgConvMatrix, iNewXdim*iNewYdim*sizeof(float), 
				"conv2d:convmatrix");
		igConvMatrixSize = iNewXdim*iNewYdim;
	}
	tmp = pfgConvMatrix;
	set_float(pfgConvMatrix, igConvMatrixSize, 0.0);
	pfPad2d(a, bCtrX, iPadY, aSizeX, aSizeY, pfgConvMatrix);
#ifdef DEBUG
	writeimage("padded.im",iNewXdim, iNewYdim, 1, pfgConvMatrix);
#endif
	if (iAvgY > 0 && iPadY > 0){
		//Compute average at top of image and place it in first row
		pfAvgRows(tmp+iNewXdim*iPadY, iNewXdim, iAvgY, tmp);
		//Copy this to the remaining rows in the pad region
		pfCopyRows(tmp, iNewXdim, iPadY-1, tmp+iNewXdim);
		//Compute average at bottom of image and place it at the bottom of
		// the padded image.
		pfAvgRows(tmp + iNewXdim*(iPadY+aSizeY-iAvgY), iNewXdim, iAvgY,
			tmp + iNewXdim*(iNewYdim-1));
		pfCopyRows(tmp + iNewXdim*(iNewYdim-1), iNewXdim, iPadY-1,
			tmp + iNewXdim*(iPadY +aSizeY));
	}
#ifdef DEBUG
	writeimage("padded_1.im",iNewXdim, iNewYdim, 1, pfgConvMatrix);
#endif
	set_float(out, aSizeX*aSizeY, 0.0);
	
	//do the convolution
	for (bY=(-bCtrY); bY<=bCtrY; bY++){
		for (aY=0; aY<aSizeY; aY++){
			if (aY-bY+iPadY>=0 && aY-bY<=aSizeY+iPadY-1){
				for (bX=(-bCtrX); bX<=bCtrX; bX++){
					for (aX=0; aX<aSizeX; aX++){
						out[a_ptr(aY, aX)] += b[b_ptr(bY+bCtrY, bX+bCtrX)] * 
							tmp[t_ptr(iPadY+aY-bY, aX-bX+bCtrX)];
					}
				}
			}
		}
	}

#undef a_ptr
#undef b_ptr
#undef t_ptr
}
