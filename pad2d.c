/*
	pad2d.c

	$Id: pad2d.c 58 2006-01-30 21:06:01Z mjs $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <mip/miputil.h>
#include "irlprivate.h"

/*	
	Revision 1.2  2004/07/27 02:40:26  frey
	Fixed bug in pfPad2d
	Change pfAvgRows so that it doesn't modify output row if iNrows is 0

	Revision 1.1  2004/07/27 02:07:00  frey
	Initial revision
*/

float *pfPad2d(float *pfInPix, 
					int iPadX, int iPadY,
					int iXdim, int iYdim, 
					float *pfOutPix)
/* Puts the 2d input array into an array array that is 
	(Xdim + 2*iPadX)x(Ydim + 2*iPadY). The columns 0..iPadX-1 and 
	(iPadX+Xdim)...(Xdim+2*iPadX-1) are filled with zeros. Similarly,
	the rows 0...iPadY-1 and (iPadY+Ydim)...(Ydim+2*iPad-1) are filled with
	zeros. If pfOutPix is not null, then the appropriate size array will be
	allocated. If pfOutPix is not null, then it must point to an appropriately
	sized output matrix.
*/
{
	int iX, iY, iNewXdim, iNewYdim;
#define InPtr(iX, iY) ((iY)*(iXdim) + (iX))
#define OutPtr(iX, iY) ((iY)*(iNewXdim) + (iX))
	iPadX = MAX(iPadX, 0);
	iPadY = MAX(iPadY, 0);
	iNewXdim = iXdim+2*iPadX;
	iNewYdim = iYdim+2*iPadY;
	if (pfOutPix == NULL)
		pfOutPix = vector(iNewXdim*iNewYdim,"Pad2d: OutPix");
	//Zero the padded rows
	for (iY=0 ; iY < iPadY; iY++){
		for (iX=0; iX<iNewXdim; iX++) {
			pfOutPix[OutPtr(iX, iY)] = 0.0;
			pfOutPix[OutPtr(iX, iY+iYdim+iPadY)] = 0.0;
		}
	}
	//Zero the padded columns
	for (iY=0 ; iY < iYdim; iY++){
		for(iX=0; iX < iPadX; ++iX){
			pfOutPix[OutPtr(iX, iY+iPadY)] = 0.0;
			pfOutPix[OutPtr(iX+iXdim+iPadX, iY+iPadY)] = 0.0;
		}
	}
	//Copy the input pixels
	for (iY=0 ; iY < iYdim; iY++){
		for (iX=0; iX<iXdim; iX++) {
			pfOutPix[OutPtr(iX+iPadX, iY+iPadY)] = pfInPix[InPtr(iX, iY)];
		}
	}
	return pfOutPix;
#undef InPtr
#undef OutPtr
}

float *pfUnPad2d(float *pfInPix, 
						int iPadX, int iPadY, 
						int iOutXdim, int iOutYdim, 
						float *pfOutPix)
{
	int iX, iY, iInXdim, iInYdim;
#define InPtr(iX, iY) ((iY)*(iInXdim) + (iX))
#define OutPtr(iX, iY) ((iY)*(iOutXdim) + (iX))
	iPadX = MAX(iPadX, 0);
	iPadY = MAX(iPadY, 0);
	iInXdim = iOutXdim+2*iPadX;
	iInYdim = iOutYdim+2*iPadY;
	if (pfOutPix == NULL)
		pfOutPix = vector(iOutXdim*iOutYdim,"UnPad2d: OutPix");
	for (iY=0 ; iY < iOutYdim; iY++){
		for (iX=0; iX<iOutXdim; iX++) {
			pfOutPix[OutPtr(iX, iY)] = pfInPix[InPtr(iX+iPadX, iY+iPadY)];
		}
	}
	return pfOutPix;
#undef InPtr
#undef OutPtr
}

float *pfAvgRows(float *pfInPix, int iXdim, int iNRows, float *pfOutPix)
//Average iNRows with iXdim elements in pfInPix and save average in vector
//pfOutPix
// If iNRows <= 0, then pfOutPix will not be touched
{
	int iX, iY;
#define InPtr(iX, iY) ((iY)*(iXdim) + (iX))
	iNRows = MAX(iNRows,0);
	if (pfOutPix == NULL){
		pfOutPix = vector(iXdim,"AvgRows: OutPix");
		set_float(pfOutPix, iXdim, 0.0);
	}

	//Nothing to average, set the output to 0
	if (iNRows == 0){
		return pfOutPix;
	}
		
	//Sum all the rows into the output vector
	for (iY=0 ; iY < iNRows; iY++){
		for (iX=0; iX<iXdim; iX++) {
			pfOutPix[iX] += pfInPix[InPtr(iX, iY)];
		}
	}

	if (iNRows > 1)
		for (iX=0; iX<iXdim; iX++) 
			pfOutPix[iX] /= (float)iNRows;
	return pfOutPix;
#undef InPtr
}

float *pfCopyRows(float *pfInPix, int iXdim, int iNRows, float *pfOutPix)
//Copy iXdim element input vector pfInPix to iNRows of rows of output vector
//pfOutPix
{
	int iX, iY;
#define OutPtr(iX, iY) ((iY)*(iXdim) + (iX))
	if (iNRows <= 0)return pfOutPix;
	if (pfOutPix == NULL)
		pfOutPix = vector(iXdim*iNRows,"AvgRows: OutPix");

	for (iY=0 ; iY < iNRows; iY++){
		for (iX=0; iX<iXdim; iX++) {
			pfOutPix[OutPtr(iX,iY)] += pfInPix[iX];
		}
	}
	return pfOutPix;
#undef OutPtr
}

#ifdef STANDALONE

#include <stdio.h>
#include <stdlib.h>
#include <mip/imgio.h>

main(int argc, char **argv)
{
	char *pchIn, *pchOut;
	int iXdim, iYdim;
	int iNewXdim, iNewYdim;
	float *pfIn, *pfPad, *pfUnPad;
	int iPadX, iPadY, iNavg;

	if (--argc != 5){
		fprintf(stderr,"usage: pad2d padx pady navg in.im out.im\n");
		exit(1);
	}
	iPadX = atoi(*++argv);
	iPadY = atoi(*++argv);
	iNavg = atoi(*++argv);
	pchIn = *++argv;
	pchOut = *++argv;
	if (iPadX < 0 || iPadY < 0){
		fprintf(stderr,"padx and pady must be nonnegative");
		exit(1);
	}

	pfIn = readimage2d(pchIn, &iXdim, &iYdim);
	pfPad = pfPad2d(pfIn, iPadX, iPadY, iXdim, iYdim, NULL);
	iNewXdim = iXdim + 2*iPadX;
	iNewYdim = iYdim + 2*iPadY;
	writeimage("pad.im", iNewXdim, iNewYdim, 1, pfPad);

	//Compute average at top of image and place it in first row of output image
	pfAvgRows(pfPad+iNewXdim*iPadY, iNewXdim, iNavg, pfPad);
	pfCopyRows(pfPad, iNewXdim, iPadY-1, pfPad + iNewXdim);

	//Compute average at bottom of image and place it in row at bottom of image
	pfAvgRows(pfPad+iNewXdim*(iPadY+iYdim-iNavg), iNewXdim, iNavg, 
					pfPad+iNewXdim*(iNewYdim-1));
	pfCopyRows(pfPad+ iNewXdim*(iNewYdim-1), iNewXdim, iPadY-1, 
					pfPad + iNewXdim*(iPadY+iYdim));

	writeimage(pchOut, iNewXdim, iNewYdim, 1, pfPad);

	pfUnPad = pfUnPad2d(pfPad, iPadX, iPadY, iXdim, iYdim, pfIn);
	writeimage("unpad.im", iXdim, iYdim, 1, pfUnPad);
}

#endif //STANDALONE
