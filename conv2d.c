/*
	conv2d.c

	$Id: conv2d.c 83 2007-05-25 18:20:44Z binhe $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <mip/miputil.h>
#include "irlprivate.h"

/*
	convolves the array a with b. The length of a is given by aSizeY * aSizeX.
	the array b has an odd*odd number of elements and the center is
	at position b[bCtrY][bCtrX]. Thus the total length of the blurring function
	is (2*bCtrY+1) * (2*bCtrX+1). The length of the output array is 
	aSizeY * aSizeX.
*/

extern float *pfgConvMatrix; 
extern int igConvMatrixSize;
void conv2d(float *a, int aSizeX, int aSizeY, 
				float *b, int bCtrX, int bCtrY,
				float *out)
{
#define a_ptr(aY, aX) ((aY) * aSizeX + (aX))
#define b_ptr(bY, bX) ((bY) * (2*bCtrX+1) + (bX))
#define t_ptr(tY, tX) ((tY) * (aSizeX+2*bCtrX) + (tX))

	int	aY, aX, bY, bX, tY, tX;



//fprintf(stderr,"conv2d: %d x %d\n",2*bCtrX+1, 2*bCtrY+1);
#ifdef TEST
	if (bCtrX == 1 && bCtrY == 1){
		float *k;
		/* special case for a 3x3 kernel */
		/*upper right (0,0) pixel*/
		out[0] = a[0] *b[4] + a[1]*b[5] + a[aSizeX]*b[7] + a[aSizeX+1]*b[8];
		/*(aSizeX-1,0) pixel*/
		out[aSizeX-1] = a[aSizeX-2]*b[3] + a[aSizeX-1]*b[4] +
							 a[2*aSizeX-2]*b[6] + a[aSizeX-2]*b[7];
		/*lower left (0,aSizeY-1) pixel*/
		out[(aSizeY-1)*aSizeX] = a[(aSizeY-2)*aSizeX]*b[1] +
										 a[(aSizeY-2)*aSizeX+1]*b[2] +
										 a[(aSizeY-1)*aSizeX]*b[4] +
										 a[(aSizeY-1)*aSizeX+1]*b[5];
		/*lower right (aSizeX-1,aSizeY-1) pixel*/
		out[(aSizeY-1)*aSizeX] = a[(aSizeY-1)*aSizeX-2]*b[0] +
										 a[(aSizeY-1)*aSizeX-1]*b[1] +
										 a[(aSizeY)*aSizeX-2]*b[3] +
										 a[(aSizeY)*aSizeX-1]*b[4];
		/* do the top row except first and last columns*/
		for(aX=1; aX<aSizeX-1; ++aX)
			out[aX] = a[aX]*b[4] + 
						 a[aX+1]*b[5] + 
						 a[aX-1]*b[3] +
						 a[aX+aSizeX]*b[7]+
						 a[aX+aSizeX+1]*b[8]+
						 a[aX+aSizeX-1]*b[6];

		/* do bottom row except 1st and last columns*/
		for(aX=1; aX<aSizeX-1; ++aX)
			out[aX+(aSizeY-1)*aSizeX] = a[(aSizeY-1)*aSizeX+aX]*b[4] + 
											 	 a[(aSizeY-1)*aSizeX+aX+1]*b[5] + 
											 	 a[(aSizeY-1)*aSizeX+aX-1]*b[3] + 
											 	 a[(aSizeY-2)*aSizeX+aX]*b[1] + 
											 	 a[(aSizeY-2)*aSizeX+aX+1]*b[2] + 
											 	 a[(aSizeY-2)*aSizeX+aX-1]*b[0];
		/* do left column except 1st and last rows*/
		for(aY=1; aY<aSizeY-1; ++aY)
			out[aY*aSizeX] = a[aY*aSizeX]*b[4] +
								  a[aY*aSizeX+1]*b[5] +
								  a[(aY+1)*aSizeX]*b[7] +
								  a[(aY+1)*aSizeX+1]*b[8] +
								  a[(aY-1)*aSizeX]*b[1] +
								  a[(aY-1)*aSizeX+1]*b[2];

		/* do right column except 1st and last rows*/
		for(aY=1; aY<aSizeY-1; ++aY)
			out[aY*aSizeX+aSizeX-1] = a[aY*aSizeX+aSizeX-1]*b[4] +
											  a[aY*aSizeX+aSizeX-2]*b[3] +
											  a[(aY+1)*aSizeX+aSizeX-1]*b[7] +
											  a[(aY+1)*aSizeX+aSizeX-2]*b[6] +
											  a[(aY-1)*aSizeX+aSizeX-1]*b[1] +
											  a[(aY-1)*aSizeX+aSizeX-2]*b[0];
		/* now we can do everything except the side columns */
		for(aY=1; aY<aSizeY-1; ++aY)
			for(aX=1; aX<aSizeX-1; ++aX)
				out[aX + aY*aSizeX] = a[aX + aY*aSizeX]*b[4];
											 /*
											 a[aX+1 + aY*aSizeX]*b[5] +
											 a[aX-1 + aY*aSizeX]*b[3] +
											 a[aX + (aY+1)*aSizeX]*b[7] +
											 a[aX+1 + (aY+1)*aSizeX]*b[8] +
											 a[aX-1 + (aY+1)*aSizeX]*b[6] +
											 a[aX + (aY-1)*aSizeX]*b[1] +
											 a[aX+1 + (aY-1)*aSizeX]*b[2] +
											 a[aX-1 + (aY-1)*aSizeX]*b[0];
											 */
	}else
#endif
	{
/* for other size b's, place the array a in tmp with bCtrX zeros on 
	left and right to avoid need to test for array bounds*/
		float *tmp;
		/* other kernel sizes */
		if (pfgConvMatrix  == NULL || 
					(igConvMatrixSize < aSizeY * (aSizeX+2*bCtrX))){
			pfgConvMatrix = (float *)
				pvIrlRealloc(pfgConvMatrix,
							aSizeY * (aSizeX+2*bCtrX)*sizeof(float), 
							"conv2d:convmatrix");
			igConvMatrixSize = aSizeY * (aSizeX+2*bCtrX);
		}
		tmp = pfgConvMatrix;

		for (tY=0; tY<aSizeY; tY++)
		for (tX=0; tX<bCtrX; tX++) {
			tmp[t_ptr(tY, tX)] = 0.0;
			tmp[t_ptr(tY, tX+aSizeX+bCtrX)] = 0.0;
		}

		for (aY=0; aY<aSizeY; aY++)
		for (aX=0; aX<aSizeX; aX++) {
			tmp[t_ptr(aY, aX+bCtrX)] = a[a_ptr(aY, aX)];
			out[a_ptr(aY, aX)] = 0.0;
		}

		for (bY=(-bCtrY); bY<=bCtrY; bY++)
		for (aY=0; aY<aSizeY; aY++)
			if (aY-bY>=0 && aY-bY<=aSizeY-1)
				for (bX=(-bCtrX); bX<=bCtrX; bX++)
				for (aX=0; aX<aSizeX; aX++)
					out[a_ptr(aY, aX)] += b[b_ptr(bY+bCtrY, bX+bCtrX)] * 
													tmp[t_ptr(aY-bY, aX-bX+bCtrX)];

	}
#undef a_ptr
#undef b_ptr
#undef t_ptr
}
