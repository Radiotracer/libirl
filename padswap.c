/*
	padswap.c

	$Id: padswap.c 58 2006-01-30 21:06:01Z mjs $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mip/miputil.h>

#include "irlprivate.h"

#ifndef MIN
#define MIN(a,b) (((a) > (b)) ? (b) : (a))
#endif                                                                          

void PadSwap(float *pfIn, /* pointer to inXdim  by inYdim input array */
				int inXdim, int inYdim, /* size of input array */
				float *pfOut, /* pointer to outXdim by outYdim output array */
				int outXdim, int outYdim /* size of output array*/)
/* swaps quadrants and pads w/zeros as necessary to put 2d arrays in 
	format required for FFT, assuming that the center (x=y=0) of
	pfIn is (xdim/2, ydim/2). The center pixel is maped to array position (0,0) 
	in the output array. The following is for a 1D array w/in and out sizes
	the same
	A B C D E=>  C D E A B
	If the output array is larger than the input, then it is padded with
	zeros. Consider the 1D case where the input size is 5 and the output
	is 7. Then we have
	A B C D E => C D E 0 0 A B 
	If the output is smaller than the input, then the central part is
	extracted. For example, if the input is 7 elements and the output is
	5, we have:
	A B C D E F G => D E F B C

	For even sizes, the quadrants are simply swapped:
	A B C D => C D A B (4 to 4)
	A B C D => C D 0 0 A B (4 to 6)
	A B C D E F => D E B C (6 to 4)
*/
{
	/* indexes in x and y direction to input and output arrays */
	int yOut, xOut, xIn, yIn; 
	int yInCtr, yOutCtr, xInCtr, xOutCtr; /* these point to the centers of the
														  in and out arrays in each direction
													  */
	int i;
	int xOutOdd, yOutOdd; /* these are 1 if the length of the output
										array in the corresponding direction is odd 
										as the 2nd part of the array has to be 
										start at iOutCtr+1*/
	int xOutShift, yOutShift;	/* amount of shift needed in corresponding 
											direction so 2nd part of Psf ends up at 
											end of array.  */
	int xInShift, yInShift;		/* amount of shift in input needed so 2nd part
											of Psf taks from the end of the array */


	xOutOdd = outXdim % 2;
	yOutOdd = outYdim % 2;

	yInCtr = inYdim/2;
	yOutCtr = outYdim/2;
	xInCtr = inXdim/2;
	xOutCtr = outXdim/2;
	/* output shifting is only needed when the output array is bigger than the
		input array */
	xOutShift = (outXdim > inXdim) ? (outXdim - inXdim)/2 + 1 : 0;
	yOutShift = (outYdim > inYdim) ? (outYdim - inYdim)/2 + 1 : 0;
	/* input shifting is needed when the output array is smaller than the input*/
	xInShift = (outXdim < inXdim) ? (inXdim - outXdim)/2 : 0;
	yInShift = (outYdim < inYdim) ? (inYdim - outYdim)/2 : 0;
	/*
	fprintf(stderr,"OutShift=(%d,%d)\n",xOutShift,yOutShift);
	*/

	/*zero output array*/
	for(i=0; i<outXdim*outYdim; i++) pfOut[i]=0.0;

	/*loop in each quadrant of the input array*/
	/*1st and 2nd quadrants, including x-axis*/
	for(yIn=yInCtr, yOut=0; 
			yOut < (yOutCtr+yOutOdd) && yIn < inYdim; 
			++yIn, ++yOut){
		/* quadrant 1, including y axis */
		for(xIn=xInCtr, xOut=0; xOut < (xOutCtr+xOutOdd) && xIn < inXdim; 
				++xIn, ++xOut){
			pfOut[xOut + outXdim*yOut] = pfIn[xIn + inXdim * yIn];
		}
		/* quadrant 2,  */
		for(xIn = xInShift, xOut=xOutCtr+xOutOdd + xOutShift; 
				xIn < xInCtr && xOut<outXdim; 
				++xIn, ++xOut){
			pfOut[xOut + outXdim*yOut] = pfIn[xIn + inXdim * yIn];
		}

	}
	/*3rd and 4th quadrants, excluding x axis */
	for(yIn = yInShift, yOut=yOutCtr+yOutOdd + yOutShift; 
			yIn < yInCtr && yOut<outYdim; 
			++yIn, ++yOut){
		/* quadrant 3, including y axis */
		for(xIn=xInCtr, xOut=0; 
				xOut < (xOutCtr+xOutOdd) && xIn < inXdim; 
				++xIn, ++xOut){
			pfOut[xOut + outXdim*yOut] = pfIn[xIn + inXdim * yIn];
		}
		/* quadrant 4,  */
		for(xIn = xInShift, xOut=xOutCtr+xOutOdd+xOutShift; 
				xIn < xInCtr && xOut<outXdim; 
				++xIn, ++xOut){
			pfOut[xOut + outXdim*yOut] = pfIn[xIn + inXdim * yIn];
		}
	}
}

void vCopy2dArray(float *pfInPix, int nx, int ny,
	float *pfOutPix, int outnx, int outny)
/* copies the image with size nx, ny into an image with size outnx, outny.
	Note that this can be used for padding or unpadding as unaccessed
	parts of the output array are set to zero. No swapping is done.
*/
{
	int ix, iy;
	int endx, endy;

	endx = MIN(nx,outnx);
	endy = MIN(ny,outny);

	for(iy=0; iy < endy; ++iy){
		for(ix=0; ix < endx; ++ix){
			pfOutPix[ix + outnx * iy ] = 
				pfInPix[ix + nx * iy];
		}
		/* zero the rest of the row */
		for(ix = endx; ix < outnx; ++ix)
			pfOutPix[ix + outnx * iy] = 0.0;
	}
	/* zero the rest of the column*/
		for(iy = endy*outnx; iy < outny*outnx; ++iy)
			pfOutPix[iy] = 0.0;
}

#ifdef STANDALONE

#include <imgio.h>

main(int argc, char **argv)
{
	int i, swapflag, newxdim, newydim, xdim, ydim, zdim; 
	char *inname, *outname;
	float *invec, *outvec;

	if (--argc != 5){
		fprintf(stderr, 
			"usage: padswap swapflag im.in newxdim newydim im.out \n");
		fprintf(stderr, "\tswapflag!=0 -> swap quadrants while padding");
		exit(1);
	}

	swapflag = atoi(*++argv);
	inname = *++argv;
	newxdim = atoi(*++argv);
	newydim = atoi(*++argv);
	outname = *++argv;

	invec = readimage2d(inname, &xdim, &ydim);
	outvec = vector(newxdim*newydim,"outvec");
	if (swapflag)
		PadSwap(invec,xdim,ydim,outvec,newxdim,newydim);
	else
		vCopy2dArray(invec,xdim,ydim,outvec,newxdim,newydim);
	writeimage(outname,newxdim,newydim,1,outvec);
	free_vector(invec);
	free_vector(outvec);
}
#endif
