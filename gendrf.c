/*
	gendrf.c
	
  $Id: gendrf.c 109 2010-04-22 16:15:01Z frey $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

/*
 *	Based on Eric's original version
 *
 * Description
 *		gendrf computes a table containing the detector point source
 *		response function. Since the table is circulary symmetric, only
 *		the portion in the first quadrant is computed. The dpsrf is the
 *		average response for a collimator with round holes having, the
 *		specified Intrinsic resolution, binned into bins with the 
 *		specified size. The table starts at a distance pfZeroDist 
 *              from the face of the collimator and each row is separated 
 *	        by the specified distance.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <mip/miputil.h>
#include <mip/printmsg.h>
#include <mip/errdefs.h>
#include <mip/fftw3d.h>

#include "irlprivate.h"
#include "drftabio.h"

#ifdef WIN32
#define hypot _hypot
#define j1 _j1
#endif

#ifdef DEBUG
#include <mip/imgio.h>
#endif

static void GenDPSRF(int Type, double Rfan, double ColToSource,
		     double ColThick, double ColToImage, double HoleDiam,
		     double BinSize, double IntrinsicFWHM, int MatSize, 
		     float *pDPSRF,int *pAliasing)
     /* computes the detector point source response function */
     /* since it is symmetric, only 1/4 of the table is computed */
     /* calculate the detector gtf in the frequency domain, Take 2D FT to give
	the dpsrf */
     /* includes sampling width and intrisnic resolution */
     /* author: eric frey. 1/28/93*/
{
  
  double intarg, samwid, normsum, normfac;
  double alpha, beta, xfreq, yfreq, freq,delf, ftmp;
  double PsfZero; /* this is the spatial extent of the psf. The product of
		     the BinSize and MatSize should be less than this */
  double NewBinSize; /*new Bin size used to avoid arg of bessel function
		       from getting too high getting to big*/
  int BinRatio; /* ratio of BinSize/NewBinSize */
  int MyMatSize; /* matrix size used in computation to avoid extra work and
		    having a matrix that is too small for the psf to fit in*/
  int NewMatSize; 
  int ix,iy, ix1, iy1;
  float *pDrfFt;/*pointer to FT of the gdpsrf*/
  int nn[2];
  
  *pAliasing = 0;
  /*Fourier transform of the Gaussian intrinsic resolution*/
  intarg = M_PI*M_PI*IntrinsicFWHM*IntrinsicFWHM/(4.0*log(2.0));

  /*the gmtf for a round hole collimator is given by 4(j1(beta*v)/(beta*v))^2*/
  alpha = (ColToSource + ColToImage + ColThick)/ColThick;
  if(Type == CONE)
    alpha *= (Rfan + ColThick)/(Rfan - ColToSource);
  beta = M_PI * alpha * HoleDiam;
  
  /*dmtf is also blurred by sampling size. This is argument for sampling
    size sinc, the FT of the sampling function*/
  samwid = M_PI * BinSize;
  
  /* We want to make sure we properly sample the gmtf. The function
     (0.25*j0(arg)/arg)^2) has a node near arg=10. The next maximum is
     about 0.2% of the value at the center. We could get rid of more sampling
     error by requiring that this argument be larger, but arg=10 seems to
     be o.k. We choose the new pixel size to be an integral fraction of the
     sampling size so we don't need interpolation later. We try to make the
     value of arg at the NyQuist Frequency, 1/(2*BinSize), be close to 10.
     We won't allow pixels larger than the sampling size as this 
     too would require interpolation
     */
  BinRatio = 20.0*BinSize/beta + 0.5;
  BinRatio = BinRatio < 1 ? 1 : BinRatio;
  NewBinSize=BinSize/BinRatio;
  
  /* we also don't want the spatial doman grf to be larger than the sampling
     array. We can estimate the size of the grf from the acceptance angle of
     the collimator using the following and then choose the new matrix size 
     accordingly. PsfZero is where the gpsf goes to zero, we need a matrix
     twice as wide as this. */
  
  PsfZero = HoleDiam*(ColToSource + ColToImage + ColThick)/ColThick;
  vPrintMsg(8,"PsfZero = %f\n",PsfZero);
  if(Type != PARALLEL){
    PsfZero *= (Rfan + ColThick)/(Rfan - ColToSource);
    vPrintMsg(8,"ColToSource:%f\n",ColToSource);
    vPrintMsg(8,"factor = %f\n",(Rfan + ColThick)/(Rfan - ColToSource));
  }
  vPrintMsg(8,"PsfZerom = %f\n",PsfZero);  
  NewMatSize = 4.0 * PsfZero/NewBinSize;
  
  /* Make the new matrix size be a the next power of 2 >= the one needed */
  MyMatSize = 2;
/* Modified 4/2/98 by Dan Kadrmas
		limited maximum MyMatSize to be 1024 here; this may cause some aliasing
		when the response is wide (ie, near the focal line for fan beam), but
		that shouldn't introduce much bias.  */
  while ((MyMatSize < NewMatSize) && (MyMatSize < 1024)) 
    MyMatSize *= 2;
  
  /* allocate memory to compute the dtf. It would be more efficient if we
     calculated using a real 2-d transform, but I don't have one so use
     a complex transform and allocate 2 times as much memory
     */
  pDrfFt = vector(2*MyMatSize*MyMatSize ,"genDPSRF: pDrfFT");
  if (pDrfFt == NULL) {
    vErrorHandler(1,1,"GenDPSRF",
		  "Error alocating memory to compute DPSF table");
  }
  
  delf = 1.0 / ((double)MyMatSize * NewBinSize);/*freq. domain sample spacing*/
  
  /* value at zero frequency is 1 */	
  pDrfFt[0] = 1;
  
  /* zero everything else just to be safe */
  for (ix=1; ix<2*MyMatSize*MyMatSize; ix++)
    pDrfFt[ix] = 0.0;
  /* only calculate gtf in first quadrant and on positive x and y axes. Use
     symmetry to compute other points
     */
  for (ix=0; ix<=MyMatSize/2; ++ix) {
    xfreq = (double)ix*delf;
    for (iy=0; iy<=MyMatSize/2; ++iy) {
      if (ix == 0 && iy == 0) continue;
      yfreq = (double)iy*delf;
      freq = hypot(yfreq, xfreq);
      /* here is the computation of the gtf and the intrinsic */
      ftmp = 2.0*j1(beta*freq)/(beta*freq);
      ftmp = ftmp*ftmp*exp(-freq*freq*intarg); /* intrinsic */
      /* now add in the sampling */
      if (samwid > 0.0)
			ftmp *= sin(samwid*freq)/(freq*samwid);
      pDrfFt[2*(MyMatSize*iy+ix)] = ftmp;
      if (iy != 0) {
			/* reflect about y-axis */
			pDrfFt[2*(MyMatSize*(MyMatSize-iy)+ix)] = ftmp;
      }
      if (ix != 0) {
			/* reflect about x-axis */
			pDrfFt[2*(MyMatSize*iy+(MyMatSize-ix))] = ftmp;
      }
      if (ix !=0 && iy != 0) {
			/* reflect through origin */
			pDrfFt[2*(MyMatSize*(MyMatSize-iy)+(MyMatSize-ix))] = ftmp;
      }
    }
  }
  
  
  /* take the nverse fourier transform */
  nn[0] = nn[1] = MyMatSize;
  invfft3d(pDrfFt,MyMatSize, MyMatSize, 1);
  
  /* This is the largest element in the matrix we need when we extract */
  NewMatSize = MIN(MyMatSize/2, BinRatio*MatSize);
  
  /* Warn user that they didn't make their table big enough */
  *pAliasing = MyMatSize > 2*MatSize;
  
  /* make sure there are no negative real values */
  for (ix=0; ix<MyMatSize*MyMatSize; ix+=2)
    if (pDrfFt[ix] < 0) 
      pDrfFt[ix]=0;
  
  /* find the sum so we can normalize to 1 */
  /* sum is value along origin, plus 2 times sum of values along positive
     x and y axes, and 4 times values inside 1st quadrant */
  /* start w/value at origin */
  normsum = pDrfFt[0];
  
  /* add values along x and y axes */
  for (ix=BinRatio; ix<NewMatSize; ix+=BinRatio)
    normsum += 2.0 * (pDrfFt[2*ix] + pDrfFt[2*ix*MyMatSize]);
  
  /* now add values in 1st quadrant */
  for (ix=BinRatio; ix<NewMatSize; ix+=BinRatio)
    for (iy=BinRatio; iy<NewMatSize; iy+=BinRatio)
      normsum += 4.0*pDrfFt[2*(ix+iy*MyMatSize)];
  
  /* now normalise, extract every  BinRatio-th pixel in the first quadrant,
     and return it pDPSRF[] */
  normfac = 1.0 / normsum;
  for (ix=0; ix<MatSize*MatSize; ++ix)
    pDPSRF[ix] = 0.0;
  
  for(ix=0, ix1=0; ix<NewMatSize; ix+=BinRatio, ix1++)
    for (iy=0, iy1=0; iy<NewMatSize; iy+=BinRatio, iy1++)
      pDPSRF[ix1+iy1*MatSize] = pDrfFt[2*(ix+iy*MyMatSize)] * normfac;
  
  /* and free the memory used for the FT. It would be more efficient not to
     reallocate this every time, but it makes other things harder. */
  free_vector(pDrfFt);
}

void gendrf(int Type, float Rfan, float fColThick, float fHoleDiam, 
	    float fGap, float fIntrinsic, float fMaxErr, float fBinSize, 
		 float fDrfSpacing, float minDist, int NumDists, int NOffs, Drf_t *DrfTab)
{
  /*
	 Computes a table of DRFs starting with minDist and going every
	 fDrfSpacing for NumDist distances. NOffs is the maximum number
	 of offsets. The number of pixels in a given direction is
	 (2*Noffs -1).   For each distance, the pixel size is
	 specified by BinSize. The collimator is assumed to have round
	 holes with diameter fHoleDiam and is fColThick thick (all units in cm).
	 Type specifies the collimator type which can be PARALLEL, FAN, or CONE.
	 fMaxErr is the maximum fractional area in the DRF that it is o.k. to
	 truncate.
    */

#define ptr(y,x) ((y)*(2*MaxOffsX-1)+(x))
  int 	iDist, iXOffs, iYOffs, MaxOffsX, MaxOffsY;
  int 	bAliasing;
  float fColToSource; 
  float *pfQuarterDpsf = vector(NOffs*NOffs,"quarterDpsf:gendrf");
  double dSum = 0.0, dNormFac, factor;
  double dMinSum = 1.0 - fMaxErr;
  int	x_shift, y_shift;
  char psf_lsf[100];
  
  for (iDist=0, fColToSource = minDist;
       iDist < NumDists; 
       iDist++, fColToSource += fDrfSpacing) {
    /* calculate the dpsf for this plane */
    
    if (fColToSource <= 0.0){
      for(iXOffs = 1; iXOffs < NOffs*NOffs; ++iXOffs)
	pfQuarterDpsf[iXOffs] = 0.0;
      pfQuarterDpsf[0] = 1.0;
    }
    else{
      GenDPSRF(Type,(double)Rfan,(double)fColToSource, (double)fColThick, 
	       (double)fGap, (double)fHoleDiam, (double)fBinSize, 
	       (double)fIntrinsic, NOffs, pfQuarterDpsf, &bAliasing);
    }
	 if (fabs(Rfan - fColToSource) < 1e-3) 
	 	factor = 1000;
	 else
		 factor = (Rfan + fColThick)/(Rfan - fColToSource);

    vPrintMsg(8," iDist= %d,size of drf type =%d\n",iDist, sizeof(Drf_t));
    
    if (bAliasing)
      vErrorHandler(2,4,"gendrf",
		  "Possible aliasing in DRF table computation at distance %.2f\n",fColToSource);
    
    if (fMaxErr <= 0.0) {
      MaxOffsX= MaxOffsY = NOffs;
      dNormFac = 1.0; 
    }
    else {
      /* 
	 find the radius required to give the specified error. Note:
	 The drf for this plane already sums to 1.0. Use only the 
	 position along the x-axis and do an integral in polar coords
	 based on the area of the rings times the value for the ring.
	 The center point in the drf is a disk with radius 0.5. 
	 */
      
      if(Type != FAN){
	int iXTempOffs = 0, iYTempOffs=1, previXOff =1;
	dSum = pfQuarterDpsf[0]; 
	iYOffs =1;
	
	for (iXOffs = 1; iXOffs < NOffs; iXOffs++) {
	  if (dSum > dMinSum) 
	    break;
	  dSum += 2*pfQuarterDpsf[iXOffs]; 
	  iYOffs = iXOffs;
	  for(iYTempOffs = 1 ;iYTempOffs <iYOffs; iYTempOffs++)
	    for(iXTempOffs =(previXOff+1); iXTempOffs<= iXOffs; iXTempOffs++)
	      dSum += 4*pfQuarterDpsf[iYTempOffs*NOffs + iXTempOffs];
	  for(iXTempOffs = 1;iXTempOffs <= iXOffs; iXTempOffs++)
	    dSum += 4*pfQuarterDpsf[iYOffs*NOffs + iXTempOffs];
	  dSum += 2*pfQuarterDpsf[iYOffs*NOffs];
	  iYOffs++;
	  previXOff = iXOffs;
	}
      }
      
      /*for fan, the response is not different for the x and y dimenssions so we
	grow a rectangular region from the center until the sum has been 
	reached. This way we get a bounding box in the form of a rectangle 
	whose x/y side ratio = the factor F+L/F-Z. 
	We do Not use an analytical ring formula*/
      
      if(Type == FAN){
	int iXTempOffs = 0, iYTempOffs=1, previXOff =1;
	dSum = pfQuarterDpsf[0];
	iYOffs= 1;
	for (iXOffs=1; iXOffs < NOffs; iXOffs++){
	  if (dSum > dMinSum)
	    break;
	  dSum += 2*pfQuarterDpsf[iXOffs];
	  if((float)iXOffs/factor - (float)iYOffs >= 0){
	    for(iYTempOffs = 1 ;iYTempOffs <iYOffs; iYTempOffs++)
	      for(iXTempOffs =(previXOff+1); iXTempOffs<= iXOffs; iXTempOffs++)
		dSum += 4*pfQuarterDpsf[iYTempOffs*NOffs + iXTempOffs];
	    for(iXTempOffs = 1;iXTempOffs <= iXOffs; iXTempOffs++)
	      dSum += 4*pfQuarterDpsf[iYOffs*NOffs + iXTempOffs];
	    dSum += 2*pfQuarterDpsf[iYOffs*NOffs];
	    iYOffs++;
	    previXOff = iXOffs;
	  }
	}  
      }
      MaxOffsX = iXOffs; /* only need to use iXOffs elements */
      if((Type == FAN) && (MaxOffsX > (int)factor)){
	MaxOffsY = (int)((float)iXOffs/factor);
	if(MaxOffsX > MaxOffsY)
	  MaxOffsY +=1;  /*just  one extra to prevent trucation in y*/
      }
      else
	MaxOffsY = iXOffs;
      vPrintMsg(8,"MaxOffsX MaxOffsY : = %d\t%d\n",MaxOffsX,MaxOffsY);
      /*
	renormalize so iXOffs elements have proper area. First compute
	actual area. Summing is a little tricky since we only have one
	quarter of the actual array 
	*/
      
      dSum = pfQuarterDpsf[0];
      if((Type == PARALLEL)||(Type == CONE)){
	for (iXOffs=1; iXOffs < MaxOffsX; ++iXOffs) {
	  /* sum along X axis, multiply by 2 to account for negative X axis */
	  dSum += pfQuarterDpsf[iXOffs]*2.0;
	  
	  /* sum along Y axis, multiply by 2 to account for negative Y axis */
	  dSum += pfQuarterDpsf[iXOffs * NOffs]*2.0;
	  
	  /* Now sum along colums excluding points on the x and y axis. 
	     Multiply these by four to account for the other 3 quadrants */
	  for (iYOffs=1; iYOffs < MaxOffsY; ++iYOffs)
	    dSum += 4.0*pfQuarterDpsf[iXOffs + iYOffs * NOffs];
	}
      }
      /*for fan type or asymmetric fan type you should sum on the
	y axis only till the margin of error is reached*/ 
      if(Type == FAN){
	/* sum along Xaxis */
	for (iXOffs = 1; iXOffs < MaxOffsX; ++iXOffs){
	  dSum += pfQuarterDpsf[iXOffs]*2.0;
	  /*sum along yaxis*/
	  if(iXOffs < MaxOffsY)
	    dSum += pfQuarterDpsf[iXOffs * NOffs] * 2.0;	
	  /*sum along columns excluding points on the x and y axis.
	    Multiply these by four to account for the other 3 quadrants */
	  for(iYOffs=1; iYOffs < MaxOffsY; ++iYOffs) 
	    dSum += 4.0*pfQuarterDpsf[iXOffs + iYOffs * NOffs];	
	}
      }
      dNormFac = 1.0 / dSum;
#ifdef DEBUG
      fprintf(stderr,"dNormFac = %f\n",dNormFac);
#endif
    }
    
    DrfTab[iDist].MaxOffsX = MaxOffsX;
    DrfTab[iDist].MaxOffsY = MaxOffsY;
    
    /* full psf */
    sprintf(psf_lsf,"DrfTab[%d].psf : gendrf",iDist);
    DrfTab[iDist].Psf = vector((2*MaxOffsX-1)*(2*MaxOffsY-1),psf_lsf);
	 DrfTab[iDist].bFFTdone=0;
    vPrintMsg(9, "address of psrf = %x\n",DrfTab[iDist].Psf);
    x_shift = MaxOffsX - 1;
    y_shift = MaxOffsY - 1;

    DrfTab[iDist].Psf[ptr(y_shift+0, x_shift+0)] = pfQuarterDpsf[0]*dNormFac;
    for (iYOffs =1; iYOffs < MaxOffsY; ++iYOffs) {
      DrfTab[iDist].Psf[ptr(y_shift+iYOffs, x_shift+0)] =
			DrfTab[iDist].Psf[ptr(y_shift-iYOffs, x_shift+0)] =
			pfQuarterDpsf[NOffs*iYOffs] * dNormFac;
    } 
    for (iXOffs=1; iXOffs < MaxOffsX; ++iXOffs)  
      DrfTab[iDist].Psf[ptr(y_shift+0, x_shift+iXOffs)] = 
			DrfTab[iDist].Psf[ptr(y_shift+0, x_shift-iXOffs)] =
				pfQuarterDpsf[iXOffs] * dNormFac;
    
    for (iYOffs =1; iYOffs < MaxOffsY; ++iYOffs) {
      for (iXOffs=1; iXOffs < MaxOffsX; ++iXOffs)
			DrfTab[iDist].Psf[ptr(y_shift + iYOffs, x_shift + iXOffs)] = 
	  		DrfTab[iDist].Psf[ptr(y_shift - iYOffs, x_shift + iXOffs)] = 
	  		DrfTab[iDist].Psf[ptr(y_shift + iYOffs, x_shift - iXOffs)] = 
	  		DrfTab[iDist].Psf[ptr(y_shift - iYOffs, x_shift - iXOffs)] = 
				pfQuarterDpsf[NOffs*iYOffs + iXOffs] * dNormFac;
    }
#ifdef DEBUG
    if((iDist==63)||(iDist==79)||(iDist==95)||(iDist==127)||(iDist==191)||
       (iDist==255)){
      sprintf(psf_lsf,"psf%d.im",iDist);
      writeimage(psf_lsf, 2*MaxOffsY-1, 2*MaxOffsX-1, 1, DrfTab[iDist].Psf);
    }
#endif
#ifdef DEBUG
    if((iDist==63)||(iDist==79)||(iDist==95)||(iDist==127)||(iDist==191)||
       (iDist==255)){
    }
#endif
  }
  free_vector(pfQuarterDpsf);
  
#undef ptr
}

DrfTab_t *psFreeDrf(DrfTab_t *psDrfTab)
{
	int iDist;
	Drf_t *psDrfs;

	if (psDrfTab != NULL){
		if (psDrfTab->psBckDrfTab != NULL){
			psDrfTab->psBckDrfTab = psFreeDrf(psDrfTab->psBckDrfTab);
			IrlFree(psDrfTab->psBckDrfTab);
		}
		psDrfs=psDrfTab->psDrfs;
		if (psDrfs != NULL){
			/* loop over each distance and free the data for it*/
			for(iDist = 0; iDist < psDrfTab->iNumDistances; ++iDist){
				free_vector(psDrfs[iDist].Psf);
			}
			IrlFree(psDrfs);
		}
		IrlFree(psDrfTab);
	}
	return NULL;
}
