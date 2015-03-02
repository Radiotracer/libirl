/*
	convolve.c

	$Id: convolve.c 58 2006-01-30 21:06:01Z mjs $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <mip/errdefs.h>
#include <mip/printmsg.h>
#include <mip/realft.h>
#include <mip/miputil.h>

#include "irlprivate.h"

#define DRFPADFAC 2

void ConvolveDrf(float *pfIn, /* input 2d slice */
					int inX, int inY, /* size of input slice */
					Drf_t *psDrf, /* Drf to use in convolution */
					int bFFTConvolve,
					float *pfOut /* image to store convolved image in */
				)
/* Do convolution of input array pfIn (inX,inY) with point spread 
		stored in the psDrf structure. If bFFTConvolve is true, the
		the convolution is done in Fourier Domain. If the transform is
		to be done in the Fourier domain, the psDrf->bFFTDone is checked.
		If it is true, this means the psDrf->Psf is already in the Fourier
		domain. If not, the PSF is FT-ed and replaces the spacial PSF.
		Note that the size is also updated and bFFTDone is set.
*/

{
	int	xdim, ydim, psf_xdim, psf_ydim;
	int iPsfX, iPsfY;
	float *PsfTmp, *pfPsf;

	/*
	fprintf(stderr,
		"ConvolveDRF: fftconvolve=%d, ffted=%d, psfsize=(%d,%d), sumpsf=%.3f\n",
		bFFTConvolve, psDrf->bFFTdone, 2*psDrf->MaxOffsX-1, 2*psDrf->MaxOffsY-1, 
		sum_float(psDrf->Psf,(2*psDrf->MaxOffsX-1)*(2*psDrf->MaxOffsY-1)));
	*/

	if (psDrf == NULL)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "DrfConvolve",
				"can't convolve using NULL DRF");
	iPsfX = psDrf->MaxOffsX;
	iPsfY = psDrf->MaxOffsY;
	pfPsf = psDrf->Psf;
	/*
	fprintf(stderr,"psf=%d, (%dx%d)",pfPsf, iPsfX, iPsfY);
	*/
	/* check bFFTConvolve, do convolution in freq domain if set */
	if(!bFFTConvolve){
		conv2d(pfIn, inX, inY, pfPsf, iPsfX-1, iPsfY-1, pfOut);
	}else{
		if(!psDrf->bFFTdone){
			/* the dimension of PSF is 2*iPsfX-1 by 2*iPsfY-1 when it is in
				spatial domain. Conv2dFFT expects the actual size of the psf 
				and a pointer to the upper left corner. So, correct for that here*/
			psf_xdim = 2 * iPsfX - 1;
			psf_ydim = 2 * iPsfY - 1;
			xdim = DRFPADFAC * inX; 
			ydim = DRFPADFAC * inY; 
			/*
			fprintf(stderr,"pad drf psf=(%dx%d) new=(%dx%d), old=(%dx%d)\n", 
							iPsfX, iPsfY, xdim, ydim, psf_xdim, psf_ydim);
			*/
			PsfTmp = vector(xdim*ydim,"PsfTmp in convolve");
			PadSwap (pfPsf,psf_xdim, psf_ydim, PsfTmp,xdim,ydim);
			free_vector(pfPsf);
			realft3d ( PsfTmp, xdim, ydim, 1);
			psDrf->bFFTdone = 1;
			psDrf->MaxOffsX = iPsfX = xdim; 
			psDrf->MaxOffsY = iPsfY = ydim; 
			pfPsf = psDrf->Psf = PsfTmp;
		}
		if( DRFPADFAC*inX != iPsfX || DRFPADFAC*inY!=iPsfY )
			vErrorHandler(ECLASS_FATAL,ETYPE_IO,"convolve",
						"input dim(%d,%d) != psf dim(%d,%d)",inX,inY,iPsfX,iPsfY);
		/* the dimension of PSF is iPsfX by iPsfY if in freq domain since
				it is zero-padded and quadrant-swapped*/
		xdim = iPsfX; 
		ydim = iPsfY; 
		conv2dPsfFFT(pfIn,inX,inY,pfPsf,xdim,ydim, pfOut);
	}
}

#ifdef STANDALONE

#include <mip/imgio.h>

int main(int argc, char **argv)
{
	int xdim1, ydim1, xdim2, ydim2; 
	int bFFTconv=0, in2_freq=0;
	char *inname1, *inname2, *outname, *outpsf;
	float *in1, *in2, *outvec;
	Drf_t Drf;

	if (--argc != 6){
	fprintf(stderr,"usage: convolve in1.im in2.im im.out flag1 flag2 outpsf\n");
	fprintf(stderr,"\tconvolution using multiplication in freq domain\n");
	fprintf(stderr,"\tflag1=1, convolution in freq domain, 0 spatial\n");
	fprintf(stderr,"\tflag2=1, in2 in freq domain already, 0 spatial\n");
	exit(1);
	}

	inname1 = *++argv;
	inname2 = *++argv;
	outname = *++argv;
	bFFTconv = atoi(*++argv);
	in2_freq = atoi(*++argv);
	outpsf = *++argv;

	in1 = readimage2d(inname1, &xdim1, &ydim1);
	in2 = readimage2d(inname2, &xdim2, &ydim2);
	Drf.Psf = in2;
	if (! in2_freq){
		Drf.MaxOffsX = xdim2/2 + 1;
		Drf.MaxOffsY = ydim2/2 + 1;
	}else{
		Drf.MaxOffsX = xdim2;
		Drf.MaxOffsY = ydim2;
		printf("input psf is %dx%d\n",xdim2, ydim2);
	}
	Drf.bFFTdone=in2_freq;
	outvec = vector(xdim1*ydim1,"outvec");
	ConvolveDrf(in1,xdim1,ydim1,&Drf,bFFTconv, outvec);
	writeimage(outname,xdim1,ydim1,1,outvec);
	writeimage(outpsf,Drf.MaxOffsX, Drf.MaxOffsY,1,Drf.Psf);
	return(0);
}
#endif
