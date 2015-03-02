/*
	conv2dPsfFFT.c

	$Id: conv2dPsfFFT.c 58 2006-01-30 21:06:01Z mjs $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <mip/miputil.h>
#include <mip/errdefs.h>
#include <mip/printmsg.h>
#include <mip/realft.h>
#include "irlprivate.h"

static void multrealft2d(int nx, int ny, float *pfIn1, float *pfIn2, float *pfOut)
/* multiplies 2 complex FTs of real functions that are stored in a special format for in-place
   real-format FTs. For an 2D FT of an input real array, the following elements are real:
		(0,0), (0,Nyquist), (Nyquist,0) and (Nyquist,Nyquist). 
	This is written so that pfOut can be the same as pfIn1 or pfIn2*/
{
	float fReal, fImag, f1, f2;
	int i;

	/* there are 4 real values that are stuffed in special places
	that require special treatment. The first two are in the first
	and second elements of the array.
	*/
	pfOut[0] = pfIn1[0]*pfIn2[0];
	pfOut[1] = pfIn1[1]*pfIn2[1];

	/*The other two are stored in real and imaginary parts of the bin
	  corresponding to the Nyquist frequency in the y-direction for 
	  x frequency=0. If we compute and store these, then we can
	  do regular complex multiplation for the entire of the array except
	  the first two elements, then overwrite with these two saved special 
	  values.
	*/
	f1 = pfIn1[(ny/2)*nx]*pfIn2[(ny/2)*nx];
	f2 = pfIn1[(ny/2)*nx+1]*pfIn2[(ny/2)*nx+1];

	for(i=2; i< nx*ny; i+=2){
		fReal = pfIn1[i]*pfIn2[i] - pfIn1[i+1]*pfIn2[i+1];
		fImag = pfIn1[i]*pfIn2[i+1] + pfIn1[i+1]*pfIn2[i];
		pfOut[i] = fReal;
		pfOut[i+1] = fImag;
	}

	/* now overwrite the two special values f1 and f2*/
	pfOut[(ny/2)*nx] = f1;
	pfOut[(ny/2)*nx+1] = f2;
}

void conv2dPsfFFT(float *pfIn, int iInXdim, int iInYdim, 
					float *pfPsfFT, int iPsfXdim, int iPsfYdim,
					float *pfOut)
/* 
	Convolves the 2D spatial-domain image in pfIn with the Fourier domain
		PSF in pfPsfFT. The pfPsfFT is assumed to be fourier tranformed
		using the format returned by realft3d applied to a 2d image with
		size iPsfXdim by iPsfYdim.

		The input image pfIn has a size that is iInXdim by iInYdim and is padded 
		or extracted to the same size as the input Psf. Note that the extraction
		is from the upper left corner and this is probably not what is wanted. 
		So, the caller should make sure that the size of the DRF that is being
		used is appropriatey padded so it is at least as large as the input
		image.

		The output image is returned in pfOut which has the same size as the 
		input.
*/
{
	float *pfPadded;

	/* allocate memory for the padded image*/

	/*
	fprintf(stderr,"conv2dPsf %d x %d\n",iPsfXdim, iPsfYdim);
	*/
	pfPadded = vector(iPsfXdim*iPsfYdim, "conv2dFFT:pfPadded");
	/* copy the image to the correct size array. This pads with zeros
		if needed and extracts if needed. The image is placed in the
		upper left hand corner
	*/
	vCopy2dArray(pfIn, iInXdim, iInYdim, pfPadded, iPsfXdim, iPsfYdim);
	//writeimage("padded.im", iPsfXdim, iPsfYdim, 1, pfPadded);
	/* take the real-to-complex FFT */
	realft3d(pfPadded, iPsfXdim, iPsfYdim, 1);

	//fprintf(stderr,"pfFT[0]=%f,%f\n",pfFT[0], pfFT[1]);
	//fprintf(stderr,"pfPsfFT[0]=%f,%f\n",pfPsfFT[0], pfPsfFT[1]);

	/* now multiply the FT of the image with the FT of the psf 
		The multiplication is stored in the input array 
	*/
	multrealft2d( iPsfXdim, iPsfYdim, pfPadded, pfPsfFT, pfPadded);
	/*
	fprintf(stderr,"filtered pfFT[0]=%f,%f\n",pfFT[0], pfFT[1]);
	*/

	/* now take the inverse FT of the product and store it in pfPadded */
	realinvft3d(pfPadded, iPsfXdim, iPsfYdim, 1);
	/*
	fprintf(stderr,"sum Filtered Image=%f\n",
		sum_float(pfPadded, iPsfXdim*iPsfYdim));
	*/
	
	/*the image is now convolved, but its size is too big. Extract the pixels
	  we need from the upper left hand corner */
	vCopy2dArray(pfPadded, iPsfXdim, iPsfYdim, pfOut, iInXdim, iInYdim);
	/*
	fprintf(stderr,"sum unpadded Filtered Image=%f\n",
		sum_float(pfOut, iInXdim*iInYdim));
	*/
	/* free all vectors */
	free_vector (pfPadded);
}

