/*
	krnlprj3d.c

	$Id: krnlprj3d.c 96 2010-01-22 05:05:40Z binhe $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <limits.h>

#include <mip/realft.h>
#include <mip/errdefs.h>
#include <mip/printmsg.h>
#include <mip/miputil.h>

#include "irlprivate.h"
#include "krnlio.h"
#include "krnlprj.h"

#ifdef WIN32
#define hypot _hypot
#endif

#ifndef MIN
#define MIN(a,b) (((a) > (b)) ? (b) : (a))
#endif

/* whine if difference in requested pixel kernel and user pixel size is
 * greater than the following.
*/
#define KRNL_PIXSIZE_DIFF_THRESH 0.02 
#define INDX3D(ix,iy,iz,nx,ny,nz) ((ix)+(nx)*((iz)+(nz)*(iy)))
#define INDX2D(ix,iy,nx,ny) ((ix)+(nx)*(iy))
static int Debug = 0;

#define vZeroVec(size, vec) set_float((vec), (size), 0.0)

static float *pfReorder3dKrnl(float *pfKrnlIn, int iInNoffs, int iInNdepths,
	int iOutNoffs, int iOutNslices, int iOutNdepths)
/* the input kernel is in the format where the fastest varying two dimensions
	are specified by noffs and represent offsets on the face of the detector.
	In these dimensions, the maximum is at 0,0. The third dimension, represents 
	the distance from the face of the collimator. The center in this direction
	is at iInNdepths/2-1 and the closest distance to the collimator is the last
	slice in the image (note the output is in the same format as ksimind
	and should not be flipped. The output image
	has size in the first 2 dimensions of iOutNoffs (direction on detector
	perpendicular to radius of rotation) by iOutNslices (direction on
	detector parallel to radius of rotation). The slices are again parallel
	to the number of depths. The major difference is that the ouptput is
	in the wraparound format needed for Fourier transforms. That is,
	we start with (0,0,0), go to positive values of the coordinate
	each dimension until we get to half the dimension, then go to the most 
	negative value of the coordinate. For the direction perpendicular
	to the collimator, we want the most negative coordinate to be the
	closest to the collimator. Note that the input and output sizes
	can be different. If the output size is bigger, then zeros are
	substituted appropriately.
*/
{
	float *pfKrnlOut;
	int ixout, iyout, izout, iyin, ixin, izin;
	int inzdim=iInNdepths;
	int inydim=iInNoffs;
	int inxdim=iInNoffs;
	int outzdim=iOutNdepths;
	int outydim=iOutNslices;
	int outxdim=iOutNoffs;
	int xmin, ymin;

	pfKrnlOut = vector(iOutNdepths*iOutNoffs*iOutNslices, 
						"Reorder3dKrnl:KrnlOut");
	vZeroVec(iOutNdepths*iOutNoffs*iOutNslices,pfKrnlOut);
	/* iz is loop over distances from collimator*/

	/* compute upper limit of number of pixels to copy in x and y dimension
	   these are the half limits (i.e. number on each size of zero, so
		the output size is 1/2 of that requested since this array contains
		both sides of the origin. The input array only contains the positive side
	*/
	xmin = MIN(inxdim,outxdim/2);
	ymin = MIN(inydim,outydim/2);
	vPrintMsg(6, "xmin,ymin=%d,%d\n",xmin,ymin);
	for(izin=0; izin < inzdim; ++izin){
		// The following takes care of shifting to the center of the new
		// Array and flipping so the detector is at small indices
		izout = outzdim/2 - izin + inzdim/2 - 1;
		// Check to make sure that the output coordinate isn't out of the
		// array
		if (izout < 0 || izout >= outzdim) continue;
		// swap the axes into Fourier ordering
		izout += izout < outzdim/2 ? outzdim/2 : -outzdim/2;
		//printf("%d -> %d (%d %d)\n",izin,izout,inzdim,outzdim);
		for(iyin=0; iyin < ymin; iyin++){
			iyout=iyin;
			for(ixin=0; ixin < xmin; ixin++){
				ixout=ixin;
				pfKrnlOut [ixout+outxdim*(iyout + outydim*izout)] =
						pfKrnlIn[ixin + inxdim*(iyin + inydim*izin)];
			}
			for(ixin=1; ixin < xmin; ixin++){
				ixout=outxdim-ixin;
				pfKrnlOut [ixout+outxdim*(iyout + outydim*izout)] =
						pfKrnlIn[ixin + inxdim*(iyin + inydim*izin)];
			}
		}
		for(iyin=1; iyin < ymin; iyin++){
			iyout=outydim - iyin;
			for(ixin=0; ixin < xmin; ixin++){
				ixout=ixin;
				pfKrnlOut [ixout+outxdim*(iyout + outydim*izout)] =
						pfKrnlIn[ixin + inxdim*(iyin + inydim*izin)];
			}
			for(ixin=1; ixin < xmin; ixin++){
				ixout=outxdim-ixin;
				pfKrnlOut [ixout+outxdim*(iyout + outydim*izout)] =
						pfKrnlIn[ixin + inxdim*(iyin + inydim*izin)];
			}
		}
	}
	return pfKrnlOut;
}

float *pfInterpPsf(float *pfRadialPsf,
						int iInNx, int iInNy, int iOutNx, int iOutNy, int iOutNz)
{
	float *pfPsf;
	int ix, iy, iz, ir, iYin;
	float fVal;
	double r, x, z, frac;

	pfPsf = vector(iOutNx*iOutNy*iOutNz, "Psf in InterpPsf");
	for(iy = 0; iy < iOutNy; iy++){
		if (iy < iOutNy/2)
			iYin = iInNy/2 + iy;
		else
			iYin = iy +iInNy/2 - iOutNy;
		if (iYin >= iInNy || iYin < 0) continue;
		for(iz = 0; iz < iOutNz; ++iz){
			if (iz > iOutNz/2)
				z = (double)iz - iOutNz;
			else
				z = (double)iz;
			/*
			z = (double)iz - (double)iOutNz/2.0;
			*/
			for(ix = 0; ix < iOutNx; ix++){
				if (ix > iOutNx/2)
					x = (double)ix - iOutNx;
				else
					x = (double)ix;
				/*
				x = (double)ix - (double)iOutNx/2.0;
				*/
				r = hypot(x,z);
				ir = (int)r;
				frac = fmod((double)r, 1.0);
				/*
				if (iz == 0 && ix == 0)
					fprintf(stderr,"%d %d %d-> %d %d %.4f %.4f\n",
						ix,iy, iz, iYin,ir, r,frac);
				*/
				if (ir >= iInNx){
					fVal = 0;
				}else if (ir == iInNx - 1){
					fVal = pfRadialPsf[ir + iYin*iInNx]*(1-frac);
				}else{
					fVal = pfRadialPsf[ir + iYin*iInNx]*(1-frac) + 
							frac*pfRadialPsf[ir+1+iYin*iInNx];
				}
				pfPsf[ix+iOutNx*(iy + iOutNy*iz)] = fVal;
			}
		}
	}
	return(pfPsf);
}

static float *pfPadHalf(float *pfKrnl, int nx, int ny)
{
	int ix, iy, nyshift, nxshift, i;
	float *pfKrnl2;

	pfKrnl2 = vector(nx*ny*8, "pfKrnl2 in pfPadHalf");
	nyshift = ny/2;
	nxshift = nx;
	for(i=0; i<nx*ny * 8; ++i)
		pfKrnl2[i] = 0.0;

	for(iy=0; iy < ny; ++iy){
		for(ix=0; ix < nx; ++ix)
			pfKrnl2[nx + nxshift + ix + (iy+nyshift)*nx*4] = 
				pfKrnl[ix + iy*nx];
		for(ix=1; ix < nx; ++ix)
			pfKrnl2[nxshift + ix + (iy+nyshift)*nx*4] = 
				pfKrnl[nx- ix + iy*nx];
	}
	free_vector(pfKrnl);
	return(pfKrnl2);
}

float *pfSwapQuadrants(float *pfKrnl, int nx, int ny)
{
	/* this routine swaps 1<->3 and 2<->4. The center is assumed to be
		at (nx/2,(ny/2) which means that (nx/2,ny/2 -> (0,0) and 
			(0,0) -> (nx-nx/2, ny-ny/2). For example, if nx=ny=64, (32,32)->(0,0)
			and (0,0) -> (32,32), for nx=ny=65, (32,32)->(0,0) and (0,0)->(33,33)*/
	int ix, iy, nxshift, nyshift;
	float *pfSwap;

	nxshift = nx/2;
	nyshift = ny/2;
	pfSwap = vector(nx*ny, "pfSwap in pfSwapQuadrants");

	/*outer loop for new quadrants 1 and 2*/
	for(iy=0; iy < ny - nyshift; ++iy){
		/*compute new quadrant 2*/
		for(ix=0; ix < nx - nxshift; ++ix)
			pfSwap[ix + iy*nx] = pfKrnl[ix+nxshift + (iy+nyshift)*nx];
		/*new quadrant 1*/
		for(ix = nx-nxshift; ix < nx; ++ix)
			pfSwap[ix + iy*nx] = pfKrnl[ix-nx+nxshift + (iy+nyshift)*nx];
	}
	/*outer loop for new quadrants 1 and 2*/
	for(iy=ny-nyshift; iy < ny; ++iy){
		/*compute new quadrant 3*/
		for(ix=0; ix < nx - nxshift; ++ix)
			pfSwap[ix + iy*nx] = pfKrnl[ix+nxshift + (iy-ny+nyshift)*nx];
		/*new quadrant 4*/
		for(ix = nx-nxshift; ix < nx; ++ix)
			pfSwap[ix + iy*nx] = pfKrnl[ix-nx+nxshift + (iy-ny+nyshift)*nx];
	}
	free_vector(pfKrnl);
	return(pfSwap);
}

static float *pfPad(float *pfKrnl, int nx, int ny)
{
	int ix, iy, nyshift, nxshift, i;
	float *pfKrnl2;


	pfKrnl2 = vector(nx*ny*4, "pfKrnl2 in pfPad");
	nyshift = ny/2;
	nxshift = nx/2;
	for(i=0; i<nx*ny * 4; ++i)
		pfKrnl2[i] = 0.0;

	for(iy=0; iy < ny; ++iy){
		for(ix=0; ix < nx; ++ix)
			pfKrnl2[nxshift + ix + (iy+nyshift)*nx*2] = 
				pfKrnl[ix + iy*nx];
	}
	free_vector(pfKrnl);
	return(pfKrnl2);
}

static int iPowerTwo(int iNum)
{
	int i;

	for(i=1; i<iNum; i *= 2)
		;
	return (i);
}


static float *pfUnPad(float *pfKrnl, int nx, int ny)
{
	int ix, iy, nyshift, nxshift, i;
	float *pfKrnl2;

	pfKrnl2 = vector(nx*ny, "pfKrnl2 in pfUnPad");
	nyshift = ny/2;
	nxshift = nx/2;
	for(i=0; i<nx*ny; ++i)
		pfKrnl2[i] = 0.0;

	for(iy=0; iy < ny; ++iy){
		for(ix=0; ix < nx; ++ix)
			pfKrnl2[ix + iy*nx] = pfKrnl[nxshift + ix + (iy+nyshift)*nx*2];
				pfKrnl[ix + iy*nx];
	}
	free_vector(pfKrnl);
	return(pfKrnl2);
}

static void vUnPadTo(float *pfKrnl, int nx, int ny, float *pfKrnl2)
/* same as pfUnPad, but instead of unpadding "in place", puts the result in 
	pfKrnl2*/
{
	int ix, iy, nyshift, nxshift, i;

	nyshift = ny/2;
	nxshift = nx/2;
	for(i=0; i<nx*ny; ++i)
		pfKrnl2[i] = 0.0;

	for(iy=0; iy < ny; ++iy){
		for(ix=0; ix < nx; ++ix)
			pfKrnl2[ix + iy*nx] = pfKrnl[nxshift + ix + (iy+nyshift)*nx*2];
				pfKrnl[ix + iy*nx];
	}
}

static void vCopy3dArray(float *pfInPix, int nx, int ny, int nz, 
	float *pfOutPix, int outnx, int outny, int outnz)
/* copies the image with size nx, ny, nz into an image with size outnx, outny,
	outnz. Note that this can be used for padding or unpadding as unaccessed
	parts of the output array are set to zero.
*/
{
	int ix, iy, iz;
	int endx, endy, endz;

	endx = MIN(nx,outnx);
	endy = MIN(ny,outny);
	endz = MIN(nz, outnz);

	for (iz=0; iz< endz; iz++){
		for(iy=0; iy < endy; ++iy){
			for(ix=0; ix < endx; ++ix){
				pfOutPix[ix + outnx * (iy + iz*outny)] = 
					pfInPix[ix + nx * (iy + iz*ny)];
			}
			/* zero the rest of the row */
			for(ix = endx; ix < outnx; ++ix)
				pfOutPix[ix + outnx * (iy + iz*outny)] = 0.0;
		}
		/* zero the rest of the column*/
			for(iy = endy*outnx; iy < outny*outnx; ++iy)
				pfOutPix[iz*outnx*outny + iy] = 0.0;
	}
	/* zero unused slices*/
	for(iz=endz*outnx*outny; iz < outnx*outny*outnz; ++iz)
		pfOutPix[iz] = 0.0;
}

static void ErrAbort(char *pchRoutine, char *pchErrMsg)
{
	vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, pchRoutine,
		pchErrMsg);
}
			
typedef struct {
	int iNumKrnlOffs;
	int iNumKrnlDepths;
	int iNumKrnlSlices;
	int iKrnlCollapseFac; /* djk:  coarse-grid collapse fac for ESSE */
	float fMu0;
	float fMuAvg;
	float fPixSize;
	float fScatFac;
	float fPriFac;
	float fCntrPixFac;
	enum {GEOM_2D, GEOM_3D} eGeom2D_3D;
	int iNumExpTerms;
	float *pfKrnlFT;
	float *pfKMuFT;
	float *pfKMu2FT;
} SCATTER_PARMS;

static struct {
	/*Working storage used in CalcEffSource*/
	int iTmpSize;
	float *pfTmp;
	float *pfActFT;
	/*Storage used in CalcEffSource and AttenProject to store summed attn coefs*/
	int iSumAtnSize;
	float *pfSumAtn;
} sGTmpStore = {0,NULL,NULL,0,NULL};

void *pvDeinitScatterPrj(void *pvScatterData)
/* frees temporary memory used by scatter projector */
{
	SCATTER_PARMS *pvScatterParms = (SCATTER_PARMS *)pvScatterData;

	if (pvScatterData == NULL)return (NULL);
	if (pvScatterParms->pfKrnlFT != NULL)free_vector(pvScatterParms->pfKrnlFT);

	if (pvScatterParms->pfKMuFT != NULL)free_vector(pvScatterParms->pfKMuFT);

	if (pvScatterParms->pfKMu2FT != NULL)free_vector(pvScatterParms->pfKMu2FT);

	sGTmpStore.iTmpSize=0;
	if (sGTmpStore.pfTmp != NULL)free_vector(sGTmpStore.pfTmp);
	sGTmpStore.pfTmp=NULL;

	if (sGTmpStore.pfActFT != NULL)free_vector(sGTmpStore.pfActFT);
	sGTmpStore.pfActFT=NULL;

	sGTmpStore.iSumAtnSize=0;
	if(sGTmpStore.pfSumAtn != NULL)free_vector(sGTmpStore.pfSumAtn);
	sGTmpStore.pfSumAtn=NULL;
	IrlFree(pvScatterData);
	return (NULL);
}

static void *pvInitScatterPrj(KrnlData_t *psKrnl,
						float pixsize,
						int iGeom2D_3D,
						int iNumExpTerms,
						int iNPix,
						int iNslices,
						int iCollapseFac,
						float fXYpadfac,
						float fZpadfac)
/* this routine reads in the kernel and mu images, puts them in the
	format needed for doing the projection, and takes the relevant FTs.
	The resulting FT images are stored in the elements of the scatter
	parmeters structure pfKrnlFT, pfKMuFT, pfMu2FT. It also stuffs
	the necessary parameters into the scatter parameters structure.
	These parameters and images are used by CalcEffScatterSource to compute
	the effective scatter source given an activity map and attenuation
	distribution. They are returned as a void pointer to the structure.
*/
{
	float *pfKrnl, *pfMu;
	float fKrnlSum, fKrnlMuSum;
	float *pfKrnl3d, *pfMu3d, *pfKMu3d, *pfKMuSq3d;
	int noffs, ndepths, nx, ny, i,nz;
	SCATTER_PARMS *psScatParms;

	if (psKrnl->fScatFac == 0.0) return NULL;
	/* save parameters for use by other related routines*/
	psScatParms = pvIrlMalloc(sizeof(SCATTER_PARMS),
			"pvInitScatterPrj:ScatParms");
	if (psScatParms == NULL)
		vErrorHandler(ECLASS_FATAL, ETYPE_MALLOC, "InitScatterPrj", "ScatParms");
	


	psScatParms->fMu0 = psKrnl->fKrnlMu0;
	psScatParms->fScatFac = psKrnl->fScatFac;
	psScatParms->fPixSize = pixsize;
	psScatParms->fPriFac = psKrnl->fPriFac;
	psScatParms->fCntrPixFac = psKrnl->fCntrPixFac;
	psScatParms->iNumExpTerms=iNumExpTerms;
	psScatParms->iKrnlCollapseFac=iCollapseFac;
	noffs = psKrnl->iNumXYoffs;
	ndepths = psKrnl->iNumDists;
	nx=iNPix / iCollapseFac;
	ny=iNPix / iCollapseFac;
	nz = iNslices / iCollapseFac;		/* kernel only needs as many slices
													as we are convolving with*/
	nx = (int)ceil(fXYpadfac*(float)(nx));
	ny = (int)ceil(fXYpadfac*(float)(ny));
	nz = (int)ceil(fZpadfac*(float)nz);
	vPrintMsg(8,"padded kernel image size=(%d,%d,%d)\n",nx,ny,nz);
	psScatParms->iNumKrnlOffs = nx;
	psScatParms->iNumKrnlDepths = ny;
	psScatParms->iNumKrnlSlices = nz;
	pfKrnl = psKrnl->pfKrnl;
	pfMu = psKrnl->pfKrnlMu;
	if (iNumExpTerms <1 || iNumExpTerms > 3){
		ErrAbort("InitScatterPrj",
			"number of terms in the exponential must be 1, 2 or 3\n");
	}
	if (iGeom2D_3D == 3){
		/*generate full 3d psf from the 1 octant version*/
		psScatParms->eGeom2D_3D = GEOM_3D;
		pfKrnl3d = pfReorder3dKrnl(pfKrnl,noffs,ndepths,nx,nz,ny);
		free_vector(pfKrnl);
		psKrnl->pfKrnl = NULL;

		pfMu3d = pfReorder3dKrnl(pfMu,noffs,ndepths,nx,nz,ny);
		free_vector(pfMu);
		psKrnl->pfKrnlMu = NULL;

		if (Debug){
#ifdef DEBUG
			writeimage("krnl3d.out.im",nx, nz, ny, pfKrnl3d);
			writeimage("krnlmu3d.out.im",nx, nz, ny, pfMu3d);
			pfKrnl = vector(nx*ny,"pfKrnl for Debugging");
			for(iy=0; iy<ny; ++iy)
				for(iz=0; iz<nz; ++iz)
					for(ix=0; ix<nx; ++ix)
						pfKrnl[INDX2D(ix,iy,nx,ny)] += 
							pfKrnl3d[INDX3D(ix,iy,iz,nx,ny,nz)];
			
			pfKrnl = pfSwapQuadrants(pfKrnl, nx, ny);
			writeimage("krnl2d.out.im",nx,ny,1,pfKrnl);
			free_vector(pfKrnl);
#endif /*DEBUG*/
		}
	}else{
		psScatParms->eGeom2D_3D = GEOM_2D;
		/*generate 2d psf*/
			ErrAbort("InitScatterPrj","2D not supported yet, sorry");
	}
	
	psScatParms->pfKrnlFT = pfKrnl3d;
	/* to compute the effective scatter source we need the product of
	the krnl times mu and mu^2. compute these here*/
	if (iNumExpTerms > 1){
		/* compute weighted average kmu */
		fKrnlSum=fKrnlMuSum=0.0;
		for(i=0; i<nx*ny*nz; ++i){
			fKrnlSum += pfKrnl3d[i];
			fKrnlMuSum += pfKrnl3d[i]*pfMu3d[i];
		}
		if (fKrnlSum == 0.0){
			vErrorHandler(ECLASS_WARN, ETYPE_ILLEGAL_VALUE,
				"IntScatterPrj","Sum of Krnl is 0!");
			psScatParms->fMuAvg = psKrnl->fKrnlMu0;
		}else{
			psScatParms->fMuAvg = fKrnlMuSum/fKrnlSum;
			vPrintMsg(5,"average mu for scatter kernel=%.4g\n",
							psScatParms->fMuAvg);
			vPrintMsg(7,"fKrnlSum=%f, fKrnlMuSum=%f\n",fKrnlSum,fKrnlMuSum);
		}
		pfKMu3d = vector(nx*ny*nz, "pfKMu3d");
		psScatParms->pfKMuFT = pfKMu3d;
		for(i=0; i<nx*ny*nz; ++i){
			pfMu3d[i] -= psScatParms->fMuAvg;
			pfKMu3d[i] = pfMu3d[i]*pfKrnl3d[i];
		}
	}else {
		psScatParms->pfKMuFT = NULL;
	}
	if (iNumExpTerms > 2){
		pfKMuSq3d = vector(nx*ny*nz, "pfKMu2");
		psScatParms->pfKMu2FT = pfKMuSq3d;
		for(i=0; i<nx*ny*nz; ++i){
			pfKMuSq3d[i] = pfMu3d[i]*pfKMu3d[i];
		}
	}else {
		psScatParms->pfKMu2FT = NULL;
	}

	vPrintMsg(2,"computing kernel-related FTs\n");
	/* note that the matrices are reorderd so the slicing is
	 	perp to y axis, so dimension in z and y dimension are switched
	*/
	realft3d(pfKrnl3d, nx, nz, ny);
	if (iNumExpTerms > 1)
		realft3d(pfKMu3d, nx, nz, ny);
	if (iNumExpTerms > 2)
		realft3d(pfKMuSq3d, nx, nz, ny);

	free_vector(pfMu3d);
	return (void *)psScatParms;
}

static void CollapseByFac(float *pfImage, int iFac, int nx, int ny, int nz,
			int *piNewnx, int *piNewny, int *piNewnz)
/* Collapses 3D image by factor iFac in all 3 dimensions.  Collapsing
is performed in-place, so that pfImage is replaced by the collapsed version */
{
	int i,j,k;
	int ii,jj,kk;
	int ni,nj,nk;
	float tmp;

	if (nx%iFac != 0 || ny%iFac != 0 || nz%iFac != 0)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "CollapseByFac",
	"factor does not evenly divide image dimensions \n     iFac=%d, nx=%d, ny=%d, nz=%d",
		iFac, nx, ny, nz);

	*piNewnx = nx / iFac;
	*piNewny = ny / iFac;
	*piNewnz = nz / iFac;

	for (kk=0; kk<*piNewnz; ++kk){
	for (jj=0; jj<*piNewny; ++jj){
	for (ii=0; ii<*piNewnx; ++ii){

		tmp = 0.0;
		for (k=0; k<iFac; ++k){
			nk = kk * iFac + k;
		for (j=0; j<iFac; ++j){
			nj = jj * iFac + j;
		for (i=0; i<iFac; ++i){
			ni = ii * iFac + i;

			tmp += pfImage[ni + nx * (nj + ny * nk)];

		}
		}
		}

		pfImage[ii + *piNewnx * (jj + *piNewny * kk)] = tmp;

	}
	}
	}

	return;
}

static void ExpandByFac(float *pfImage, int iFac, int nx, int ny, int nz,
			int *piNewnx, int *piNewny, int *piNewnz)
/* Expands 3D image by factor iFac in all 3 dimensions.  Expansion
is performed in-place, so that pfImage is replaced by the expanded version */
/* Expansion is followed by (iFac+1)-point smoothing in each direction */
{
	int i,j,k;
	int ii,jj,kk;
	int ni,nj,nk;
	float tmp;

	*piNewnx = nx * iFac;
	*piNewny = ny * iFac;
	*piNewnz = nz * iFac;

	for (kk=nz-1; kk>=0; --kk){
	for (jj=ny-1; jj>=0; --jj){
	for (ii=nx-1; ii>=0; --ii){

		tmp = pfImage[ii + nx * (jj + ny * kk)] / iFac / iFac / iFac;

		for (k=0; k<iFac; ++k){
			nk = kk * iFac + k;
		for (j=0; j<iFac; ++j){
			nj = jj * iFac + j;
		for (i=0; i<iFac; ++i){
			ni = ii * iFac + i;

			pfImage[ni + *piNewnx * (nj + *piNewny * nk)] = tmp;

		}
		}
		}
	}
	}
	}

/* now smooth the expanded image by an (2*jj+1)-point rect window */
/* (this reduces 'blockiness' due to large pixel size)            */
	jj=iFac/2;

	for (k=jj; k<*piNewnz-jj; ++k){
	for (j=0; j<*piNewny; ++j){
	for (i=0; i<*piNewnx; ++i){
	 tmp=0.0;
	 for (kk=k-jj; kk<=k+jj; ++kk){
	  tmp+=pfImage[i+*piNewnx*(j+*piNewny*kk)];
	 }
	 pfImage[i+*piNewnx*(j+*piNewny*k)]=tmp/(2*jj+1);
	}
	}
	}

	for (j=jj; j<*piNewny-jj; ++j){
	for (k=0; k<*piNewnz; ++k){
	for (i=0; i<*piNewnx; ++i){
	 tmp=0.0;
	 for (kk=j-jj; kk<=j+jj; ++kk){
	  tmp+=pfImage[i+*piNewnx*(kk+*piNewny*k)];
	 }
	 pfImage[i+*piNewnx*(j+*piNewny*k)]=tmp/(2*jj+1);
	}
	}
	}

	for (i=jj; i<*piNewnx-jj; ++i){
	for (j=0; j<*piNewny; ++j){
	for (k=0; k<*piNewnz; ++k){
	 tmp=0.0;
	 for (kk=i-jj; kk<=i+jj; ++kk){
	  tmp+=pfImage[kk+*piNewnx*(j+*piNewny*k)];
	 }
	 pfImage[i+*piNewnx*(j+*piNewny*k)]=tmp/(2*jj+1);
	}
	}
	}

	return;
}


void CalcEffScatterSource(void *pvScatterData,
							float *pfAct, 
							float *pfMap, 
							int nx, 
							int ny, 
							int nz, 
							float *pfEffSrc)
/* calculates the effective scatter source using the krnl data saved
by InitScatterPrj. The returned source has the same size as the input
map and activity */
{
	float fAtnCnvFac, fDepth;
	float *pfLActFT;
	float *pfLTmp;
	float *pfLSumAtn;
	float fMu;
	int i, ix, iy, iz, iNumKrnlPix, iNumPix, iNumSlicePix;
	int iMaxNumPix;
	int iKrnlNx, iKrnlNy, iKrnlNz;
	int nxnew,nynew,nznew,mx,my,mz;
	int iNumExpTerms, iCllpseFac;
	SCATTER_PARMS *psScatParms=pvScatterData;

	if (pvScatterData == NULL){
		set_float(pfEffSrc,nx*ny*nz,0.0);
		return;
	}

	if (psScatParms->fScatFac == 0.0){
		set_float(pfEffSrc,nx*ny*nz,0.0);
		return;
	}

	fAtnCnvFac = 1.0/(psScatParms->fMu0);
	iNumExpTerms = psScatParms->iNumExpTerms;

	if (psScatParms->pfKrnlFT == NULL || 
		(iNumExpTerms > 1 && psScatParms->pfKMuFT == NULL) || 
			(iNumExpTerms > 2 && psScatParms->pfKMu2FT == NULL))
		ErrAbort("CalcEffScatterSource",
			"Must initialize using InitScatterPrj first");
	
	iCllpseFac = psScatParms->iKrnlCollapseFac;
	vPrintMsg(8,"Collapse Factor=%d\n",iCllpseFac);

	if(nx > psScatParms->iNumKrnlOffs*iCllpseFac || 
		ny > psScatParms->iNumKrnlDepths*iCllpseFac ||
		nz > psScatParms->iNumKrnlSlices*iCllpseFac
		)
		ErrAbort("CalcEffScatterSource",
			"Collapsed activity image can not be larger than kernel image");
	
	iNumPix = nx*ny*nz;
	iNumSlicePix = nx*nz;
	iKrnlNx = psScatParms->iNumKrnlOffs ;
	iKrnlNy = psScatParms->iNumKrnlDepths; 
	iKrnlNz = psScatParms->iNumKrnlSlices;
	iMaxNumPix=iNumKrnlPix = iKrnlNx*iKrnlNy*iKrnlNz;
	if (iNumPix > iMaxNumPix) iMaxNumPix=iNumPix;

	/*Allocate working storage if it hasn't already been allocated*/
	if (sGTmpStore.pfTmp != NULL && sGTmpStore.iTmpSize < iMaxNumPix){
		free_vector(sGTmpStore.pfTmp);
		vPrintMsg(3,"CalcEffScatterSource: warning, change in image size\n");
		sGTmpStore.pfTmp = NULL;
		sGTmpStore.iTmpSize = 0;
		free_vector(sGTmpStore.pfActFT);
	}
	if (sGTmpStore.pfTmp == NULL){
		sGTmpStore.pfTmp = vector(iMaxNumPix, "CalcEffScatterSource: pfTmp");
		sGTmpStore.pfActFT = vector(iMaxNumPix, "pfActFT");
		sGTmpStore.iTmpSize = iMaxNumPix;
	}

	if (sGTmpStore.pfSumAtn != NULL && sGTmpStore.iSumAtnSize  < iNumSlicePix){
		free_vector(sGTmpStore.pfSumAtn);
		sGTmpStore.pfSumAtn = NULL;
		sGTmpStore.iSumAtnSize = 0;
	}
	if (sGTmpStore.pfSumAtn == NULL){
		sGTmpStore.pfSumAtn = 
			vector(iNumSlicePix, "CalcEffScatterSource: pfSumAtn");
		sGTmpStore.iSumAtnSize = iNumSlicePix;
	}
	pfLTmp=sGTmpStore.pfTmp;
	pfLActFT=sGTmpStore.pfActFT;
	pfLSumAtn = sGTmpStore.pfSumAtn;

/* djk:  To retain quality of pfAct, must make copy before collapsing */
	memcpy(pfLTmp,pfAct,iNumPix*sizeof(float));

/* Collapse activity image here by iKrnlCollapseFactor */
	if (iCllpseFac != 1){
		CollapseByFac(pfLTmp, iCllpseFac, nx, nz, ny,
			&nxnew, &nznew, &nynew);
	} else {
		nxnew=nx;
		nynew=ny;
		nznew=nz;
	}

	vCopy3dArray(pfLTmp, nxnew, nznew, nynew, 
			pfLActFT, iKrnlNx, iKrnlNz, iKrnlNy);
#ifdef DEBUG
	if (Debug) 
		writeimage("act.collandpad.im", iKrnlNx, iKrnlNz, iKrnlNy, pfLActFT);
#endif
	realft3d(pfLActFT, iKrnlNx, iKrnlNz, iKrnlNy);

	/* zeroth order term in taylor series */
	vPrintMsg(8,"krnl=(%d,%d,%d), image=(%d,%d,%d), collapsed image=(%d,%d,%d)\n",
		iKrnlNx, iKrnlNy, iKrnlNz, nx, ny, nz, nxnew, nynew, nznew);
	multrealft3d(iKrnlNx, iKrnlNz, iKrnlNy, psScatParms->pfKrnlFT, 
					pfLActFT, pfLTmp);
	realinvft3d(pfLTmp, iKrnlNx, iKrnlNz, iKrnlNy);

	if (Debug)
		vPrintMsg(9,"sum of tmp0 before vCopy=%f\n", 
			sum_float(pfLTmp,iKrnlNx*iKrnlNy*iKrnlNz));

	vCopy3dArray(pfLTmp,iKrnlNx,iKrnlNz,iKrnlNy,pfEffSrc,nxnew,nznew,nynew);
	if (iCllpseFac != 1){
		ExpandByFac(pfEffSrc, iCllpseFac, nxnew, nznew, nynew,
		&mx, &mz, &my);
		if (mx != nx || my != ny || mz != nz)
		 ErrAbort("CalcEffScatterSource","Expanded sizes did not match original");
	}
#ifdef DEBUG
	if (Debug){ 
		vPrintMsg(9,"sum of tmp0=%f, sum of effsrc0=%f\n", 
			sum_float(pfLTmp, iKrnlNx*iKrnlNy*iKrnlNz),
			sum_float(pfEffSrc, nx*ny*nz));
		writeimage("effsrc0.im", nx,nz,ny,pfEffSrc);
	}
#endif

	if (iNumExpTerms > 1){
		/*do the first order term in the taylor series and add it to the zeroth
		order term*/
		multrealft3d(iKrnlNx, iKrnlNz, iKrnlNy, psScatParms->pfKMuFT, 
				pfLActFT, pfLTmp);
		realinvft3d(pfLTmp, iKrnlNx, iKrnlNz, iKrnlNy);
		vCopy3dArray(pfLTmp,iKrnlNx, iKrnlNz, iKrnlNy, pfLTmp, nxnew, nznew, nynew);
		if (iCllpseFac != 1){
			ExpandByFac(pfLTmp, iCllpseFac, nxnew, nznew, nynew,
			&mx, &mz, &my);
			if (mx != nx || my != ny || mz != nz)
			 ErrAbort("CalcEffScatterSource","Expanded sizes did not match original");
		}

		vZeroVec(iNumSlicePix,pfLSumAtn);
		for(iy=0; iy < ny; ++iy){
			for(iz = 0; iz < nz; ++iz){
				for(ix=0; ix < nx; ++ix){
					pfLSumAtn[INDX2D(ix,iz,nx,nz)] += 
						0.5*pfMap[INDX3D(ix,iy,iz,nx,ny,nz)];
				}
			}
			for(iz = 0; iz < nz; ++iz){
				for(ix=0; ix < nx; ++ix){
					fDepth = pfLSumAtn[INDX2D(ix,iz,nx,nz)] *fAtnCnvFac;
					pfEffSrc[INDX3D(ix,iy,iz,nx,ny,nz)] -= 
						fDepth * pfLTmp[INDX3D(ix,iy,iz,nx,ny,nz)];
					pfLTmp[INDX3D(ix,iy,iz,nx,ny,nz)]*= fDepth;
				}
			}
			for(iz = 0; iz < nz; ++iz){
				for(ix=0; ix < nx; ++ix){
					pfLSumAtn[INDX2D(ix,iz,nx,nz)] += 
						0.5*pfMap[INDX3D(ix,iy,iz,nx,ny,nz)];
				}
			}
		}
#ifdef DEBUG
		if (Debug){
			writeimage("effsrc1.im", nx,nz,ny,pfEffSrc);
		}
#endif
	}

	if (iNumExpTerms > 2){
		/* now do the second order term in the Taylor series  and add it to
			the 0th and 1st order terms */
		multrealft3d(iKrnlNx, iKrnlNz, iKrnlNy, psScatParms->pfKMu2FT, 
				pfLActFT, pfLTmp);
		realinvft3d(pfLTmp, iKrnlNx, iKrnlNz, iKrnlNy);
		vCopy3dArray(pfLTmp, iKrnlNx, iKrnlNz, iKrnlNy, pfLTmp, nxnew, nznew, nynew);
		if (iCllpseFac != 1){
			ExpandByFac(pfLTmp, iCllpseFac, nxnew, nznew, nynew,
			&mx, &mz, &my);
			if (mx != nx || my != ny || mz != nz)
			 ErrAbort("CalcEffScatterSource",
			 	"Expanded sizes did not match original");
		}

		vZeroVec(iNumSlicePix,pfLSumAtn);
		for(iy=0; iy < ny; ++iy){
			for(iz = 0; iz < nz; ++iz){
				for(ix=0; ix < nx; ++ix){
					pfLSumAtn[INDX2D(ix,iz,nx,nz)] += 
						0.5*pfMap[INDX3D(ix,iy,iz,nx,ny,nz)];
				}
			}
			for(iz = 0; iz < nz; ++iz){
				for(ix=0; ix < nx; ++ix){
					fDepth = pfLSumAtn[INDX2D(ix,iz,nx,nz)] *fAtnCnvFac;
					pfEffSrc[INDX3D(ix,iy,iz,nx,ny,nz)] += 
						0.5 * fDepth * fDepth * pfLTmp[INDX3D(ix,iy,iz,nx,ny,nz)];
					pfLTmp[INDX3D(ix,iy,iz,nx,ny,nz)] *= fDepth*fDepth;
				}
			}
			for(iz = 0; iz < nz; ++iz){
				for(ix=0; ix < nx; ++ix){
					pfLSumAtn[INDX2D(ix,iz,nx,nz)] += 
						0.5*pfMap[INDX3D(ix,iy,iz,nx,ny,nz)];
				}
			}
		}
#ifdef DEBUG
		if (Debug){
			writeimage("effsrc2.im", nx,nz,ny,pfEffSrc);
		}
#endif
	}

	/* since we factored out an exp(-fMuAvg*d), and the attenuated
		projection only accounts for exp(-fMu0*d), do an attenuated
		projection to account for this difference*/
	fMu = psScatParms->fMu0 - psScatParms->fMuAvg;
	vZeroVec(iNumSlicePix,pfLSumAtn);
	for(iy=0; iy < ny; ++iy){
		for(iz = 0; iz < nz; ++iz){
			for(ix=0; ix < nx; ++ix){
				pfLSumAtn[INDX2D(ix,iz,nx,nz)] += 
					0.5*pfMap[INDX3D(ix,iy,iz,nx,ny,nz)];
			}
		}
		for(iz = 0; iz < nz; ++iz){
			for(ix=0; ix < nx; ++ix){
				fDepth = pfLSumAtn[INDX2D(ix,iz,nx,nz)] *fAtnCnvFac;
				pfEffSrc[INDX3D(ix,iy,iz,nx,ny,nz)] *= exp(fDepth*fMu);
			}
		}
		for(iz = 0; iz < nz; ++iz){
			for(ix=0; ix < nx; ++ix){
				pfLSumAtn[INDX2D(ix,iz,nx,nz)] += 
					0.5*pfMap[INDX3D(ix,iy,iz,nx,ny,nz)];
			}
		}
	}

	/* now add a scaled version of the activity distribution to account
		for "peak flattening" when using a collapsed kernel, if required */
	if (psScatParms->fCntrPixFac > 0){
	 vPrintMsg(4,"Using primary factor adjustment for SRF modeling: %f\n",
	 	psScatParms->fCntrPixFac);
	 for(i=0; i < iNumPix; ++i)
		pfEffSrc[i] += pfAct[i] * psScatParms->fCntrPixFac;
	}


	/* Now multiply each point in the effective scatter source by the density at
	 	the point. We can get the density from the attenuation map by dividing by
		the attenuation coefficient of water in units of per pixel*/
	/* fAtnCnvFac is already 1/mu0 in cm. Covert to inv. pixels */
	fAtnCnvFac /= psScatParms->fPixSize;
	vPrintMsg(9,"atncnvfac=%f before scaling source\n", fAtnCnvFac);
	vPrintMsg(9,"mu0=%.5f \n", psScatParms->fMu0);
	vPrintMsg(9,"pixsize=%.5f \n", psScatParms->fPixSize);
	for(i=0; i < iNumPix; ++i){
		pfEffSrc[i] *= pfMap[i] * fAtnCnvFac;
		if (pfEffSrc[i] < 0) pfEffSrc[i] = 0.0;
	}

#ifdef DEBUG
	if (Debug)
		writeimage("effsrc.out.im", nx, nz, ny, pfEffSrc);
#endif
}

void *pvInitESSEScatterModel(char *pchKrnlFile,/* name of file containing krnl*/
							char *pchParString, /*string containing 
															xypad zpad expterms geom*/
							float fPixsize, /*size of pixels*/
							int nPix,
							int nSlices,
							int iCollapseFac,
							int bDebugFlag)
{
	float fXYPadFac,fZPadFac;
	int iGeom, iNumExpTerms;
	void *pvScatData;
	KrnlData_t *psKrnl;

	Debug = bDebugFlag;

	vPrintMsg(2,"Initializing ESSE scatter model\n");
	vPrintMsg(9,"parstr=%s\n",pchParString);
	if (sscanf(pchParString,"%f %f %d %d",
				&fXYPadFac, &fZPadFac, &iNumExpTerms, &iGeom) <0){
		vPrintMsg(1,"Error: for ESSE, esse_parms  should contain:");
		vPrintMsg(1,"	 xypadfac, zpadfac, numexpterms\n");
		ErrAbort("InitESSEScatterModel","illegal srf par string");
	}
	vPrintMsg(8,"Collapse Factor=%d\n",iCollapseFac);
	vPrintMsg(9, "xypadfactor=%.2f, zpadfac=%.2f, use %d terms\n",
			fXYPadFac, fZPadFac, iNumExpTerms);
	vPrintMsg(9, "geom=%dd krnlfile=%s \n",
			iGeom, pchKrnlFile);
	if (iGeom!=2 && iGeom != 3)
		ErrAbort("InitESSEScatterModel","geometry must be 2 or 3");
	if (iNumExpTerms<1 || iNumExpTerms >3)
		ErrAbort("InitESSEScatterModel","numexpterms must be 1, 2, or 3");

	psKrnl = psReadKrnlFile(pchKrnlFile,1);
	if (psKrnl == NULL)
		ErrAbort("InitESSEScatterModel","Error reading kernel file");
	pvScatData = pvInitScatterPrj(psKrnl, fPixsize, 
											iGeom,iNumExpTerms,
											nPix,
											nSlices,iCollapseFac,
											fXYPadFac, fZPadFac);
	if (fabs(psKrnl->fPixSize - fPixsize * iCollapseFac) >
				KRNL_PIXSIZE_DIFF_THRESH){
		vErrorHandler(ECLASS_WARN, ETYPE_INFO, "vInitESSEScatterModel",
			"Pixel size after collapse=%.3g, Pixel size in kernel=%.3g\n",
				iCollapseFac*fPixsize, psKrnl->fPixSize);
	}
	vPrintMsg(2,"done initializing ESSE model\n");
	/* free the psKrnl structure and its data. Note that pvInitScatterPrj
		has allready freed pfKrnl and pfMu */
	psKrnl=psFreeKrnlData(psKrnl,0);
	return (pvScatData);
}

#ifdef STANDALONE

#include <sys/times.h>
#include <unistd.h>

#include "imgio.h"
void AttenProject(void *pvScatterData, float *pfAct, float *pfEffSrc, 
		float *pfMap, int nx, int ny, int nz, float *pfPrj)
{
	int iNumSlicePix;
	int ix, iy, iz;
	float *pfLSumAtn;
	SCATTER_PARMS *psScatParms=(SCATTER_PARMS *)pvScatterData;

	/* allocate temporary memory for new array if it hasn't been allocated 
		before */
	iNumSlicePix = nx*nz;
	if (sGTmpStore.pfSumAtn != NULL && sGTmpStore.iSumAtnSize  < iNumSlicePix){
		free_vector(sGTmpStore.pfSumAtn);
		sGTmpStore.pfSumAtn == NULL;
		sGTmpStore.iSumAtnSize = 0;
	}
	if (sGTmpStore.pfSumAtn == NULL){
		sGTmpStore.pfSumAtn = 
			vector(iNumSlicePix, "AttenProject: pfSumAtn");
			sGTmpStore.iSumAtnSize = iNumSlicePix;
	}
	pfLSumAtn=sGTmpStore.pfSumAtn;

	for(ix = 0; ix < iNumSlicePix; ++ix){
			pfLSumAtn[ix] = 0.0;
		pfPrj[ix] = 0.0;
	}

	for(iy=0; iy < ny; ++iy){
		for(iz = 0; iz < nz; ++iz){
			for(ix=0; ix < nx; ++ix){
				pfLSumAtn[INDX2D(ix,iz,nx,nz)] += 
					0.5*pfMap[INDX3D(ix,iy,iz,nx,ny,nz)];
			}
		}
		for(iz = 0; iz < nz; ++iz){
			for(ix=0; ix < nx; ++ix){
				pfPrj[INDX2D(ix,iz,nx,nz)] += exp(-pfLSumAtn[INDX2D(ix,iz,nx,nz)])*
						(
							
								psScatParms->fPriFac * 
									pfAct[INDX3D(ix,iy,iz,nx,ny,nz)]
								+ psScatParms->fScatFac *
									pfEffSrc[INDX3D(ix,iy,iz,nx,ny,nz)]
						);
			}
		}
		for(iz = 0; iz < nz; ++iz){
			for(ix=0; ix < nx; ++ix){
				pfLSumAtn[INDX2D(ix,iz,nx,nz)] += 
					0.5*pfMap[INDX3D(ix,iy,iz,nx,ny,nz)];
			}
		}
	}
}

void vIllegalUsageAbort()
{
	fprintf(stderr, "usage: [-d] [-2] krnlprj krnl.par line mu0 pixsize clpsfac act map out\n");
		exit(1);
}

float *local_reorder(int zdim, int ydim, int xdim, float *Matrix, float *TempMatrix){
  int i,j,k;
  if (TempMatrix == NULL)
  	TempMatrix=vector(xdim*ydim*zdim,"reorder: TempMatrix");
  for(j=0;j<ydim;j++)
	 for(k=0;k<zdim;k++)
		for(i=0;i<xdim;i++)
		  TempMatrix[j*zdim*xdim +k*xdim + i] = Matrix[k*ydim*xdim + j*xdim + i];
  return (TempMatrix);
}

main(int argc, char **argv)
{
	void *pvScatterData;
	float *pfAct, *pfMap;
	float *pfPrj, *pfEffSrc;
	float *pfTmp;
	float mu0, pixsize;
	char *pchKrnlPars, *pchAct, *pchMap, *pchOut;
	char *esse_pars="1 1 3 3"; /*xypad zpad numexp dim*/
	char *pch;
	int iKrnlLine, iCollapseFac=1;
	int nx, ny, nz, nx_act, ny_act, nz_act;
	int i;
	int iNumExpTerms;
	int iGeom2D_3D = 3;
	int nxshift, nyshift;
	FILE *fpOut;
	struct tms sTimes;
	clock_t iLastTime, iTime;
	float fSecsPerTick;

	vSetMsgLevel(8);
	for(argc--, ++argv ;argc > 0 && **argv == '-'; argc--, ++argv){
		for(pch = *argv+1; *pch != '\0'; ++pch){
			switch (*pch){
				case 'd': 
					Debug=1;
					break;
				case '2':
					iGeom2D_3D = 2;
					ErrAbort("krnlprj","Sorry, 2D is not supported yet\n");
					break;
				default:
					fprintf(stderr,"krnlprj: illegal option: %c\n",*pch);
					vIllegalUsageAbort();
					break;
			}
		}
	}
					
	if (argc != 8) {
		vIllegalUsageAbort();
	}
	pchKrnlPars=*argv;
	iKrnlLine=atoi(*++argv);
	mu0 = atof(*++argv);
	if (mu0 <= 0)ErrAbort("krnlprj3d", "mu0 must be > 0");
	pixsize = atof(*++argv);
	if (pixsize <= 0)ErrAbort("krnlprj3d", "pixelsize must be > 0");
	fprintf(stderr,"pixsize=%f\n",pixsize);
	iCollapseFac = atoi(*++argv);
	pchAct = *++argv;
	pchMap = *++argv;
	pchOut = *++argv;

	pfTmp = readimage3d(pchAct, &nx_act, &ny_act, &nz_act);
	pfAct = local_reorder(nz_act, ny_act, nx_act, pfTmp, NULL);
	free_vector(pfTmp);
	pfTmp = readimage3d(pchMap, &nx, &ny, &nz);
	pfMap = local_reorder(nz, ny, nx, pfTmp, NULL);
	if (nx != nx_act || ny != ny_act || nz != nz_act)
		ErrAbort("krnlprj", 
			"Attenuation map and activity must be same size\n");

	fSecsPerTick = 1.0/((float)sysconf(_SC_CLK_TCK));
	times(&sTimes);
	iLastTime=sTimes.tms_utime;
	
   pvScatterData=pvInitESSEScatterModel(pchKrnlPars, iKrnlLine,
			  esse_pars,
	        pixsize, mu0, 
	         nz, iCollapseFac, Debug);
	times(&sTimes);
	fprintf(stderr,"cpu time to initialize scatter estimate=%.2fs\n",
		(sTimes.tms_utime - iLastTime)*fSecsPerTick);

	pfPrj = vector(nx*nz, "pfPrj");
	pfEffSrc = vector(nx*ny*nz, "pfEffSrc");
	/* compute the effective scatter source*/
	fprintf(stderr,"activity image size=(%d,%d,%d)\n",nx,ny,nz);

	times(&sTimes);
	iLastTime=sTimes.tms_utime;
	CalcEffScatterSource(pvScatterData,pfAct, pfMap, nx, ny, nz, pfEffSrc);
	times(&sTimes);
	fprintf(stderr,"cpu time to compute scatter source=%.2fs\n",
		(sTimes.tms_utime - iLastTime)*fSecsPerTick);
	iLastTime = sTimes.tms_utime;
	/* now do the attenuated projecton  of the effective source*/
	AttenProject(pvScatterData, pfAct, pfEffSrc, pfMap, nx, ny, nz, pfPrj);
	times(&sTimes);
	fprintf(stderr,"cpu time for attenuated projection=%.2fs\n",
		(sTimes.tms_utime - iLastTime)*fSecsPerTick);
	writeimage(pchOut,nx,nz,1,pfPrj);
	
	pvScatterData=pvDeinitScatterPrj(pvScatterData);
	free_vector(pfPrj);
	free_vector(pfAct);
	free_vector(pfMap);
	free_vector(pfEffSrc);
	times(&sTimes);
	fprintf(stderr,"total cpu time=%.2fs\n",
		sTimes.tms_utime*fSecsPerTick);
	exit(0);
}
#endif
