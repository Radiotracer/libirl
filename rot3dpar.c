/*
	rot3d.par
	
	$Id: rot3dpar.c 58 2006-01-30 21:06:01Z mjs $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <math.h>

#include <mip/miputil.h>
#include <mip/errdefs.h>

#include "irlprivate.h"

#define EPS 0.00001

void rot3dpar (float fAngle, float *OldMatrix, float *NewMatrix,
	     int PB_FLAG, int NumSlices, int NumPixels, int NumRotPixs)
{
#define o_ptr(OldY,Slice,OldX) ((OldY)*NumSlices*OldDim+(Slice)*OldDim+(OldX))
#define n_ptr(NewY,Slice,NewX) ((NewY)*NumSlices*NewDim+(Slice)*NewDim+(NewX))
/*The followign macros are for rotation about the z axis */
/*#define o_ptr(OldY,Slice,OldX) ((Slice)*OldDim*OldDim + (OldY)*OldDim + (OldX))
#define n_ptr(NewY,Slice,NewX) ((Slice)*NewDim*NewDim + (NewY)*NewDim + (NewX))
*/
#ifdef FAST
#define f_ptr(NewY,Slice,NewX) ((NewY)*NumSlices*(NewDim+2)+(Slice)*(NewDim+2)+(NewX))
#endif
  
  int	Slice, NewDim,OldDim,NewX, NewY, OldX, OldY;
  float FloatX, FloatY, X, Y;
  float z0, z1, z2, z3, z01, z23, z0123;
  float NewAxisX, NewAxisY, OldAxisX, OldAxisY;
  int   z0out, z1out, z2out, z3out; /* flags indicating out of old matrix */
  int	HeadPixel, TailPixel, NumPix;
  float Radius;
  float Sin, Cos;
  
  
  /***************************************************************************
    rotation formula:
    
    FloatX = (NewX-NewAxisX)*Cos - (NewY-NewAxisY)*Sin + OldAxisX;
    FloatY = (NewX-NewAxisX)*Sin + (NewY-NewAxisY)*Cos + OldAxisY;
    *************************************************************************/
  /*fprintf(stderr,"PB_FLAG = %d\n",PB_FLAG);*/
  
  
  Sin=(float)sin(fAngle);
  Cos=(float)cos(fAngle);
  /*
  fprintf(stderr,"angle=%.1f, sin=%.4f, cos=%.4f\n", fAngle*180.0/M_PI,
		  Sin, Cos);
  */
  if(PB_FLAG == 0){
    OldDim = NumPixels;
    NewDim = NumRotPixs;
  }else{
    OldDim = NumRotPixs;
    NewDim = NumRotPixs;
  }
  set_float(NewMatrix,NewDim*NewDim*NumSlices,0.0);
  
#ifdef FAST
  NewDim -= 2;
  vPrintMsg(9,"FAST defined\n");
  /*assume the size of the NewMatrix is (NewDim+2)*(NewDim+2)*NumSlices, i.e.
    each slice has 4 edges(1 pixel wide per edge) on top, bottom, left, right,
    which will not be used in the rotation computation. */
  
  /* set 4 corners to zero */
  for (Slice=0; Slice<NumSlices; Slice++)
    NewMatrix[f_ptr(0, Slice, 0)] =
      NewMatrix[f_ptr(NewDim+1, Slice, NewDim+1)] =
      NewMatrix[f_ptr(0, Slice, NewDim+1)] =
      NewMatrix[f_ptr(NewDim+1, Slice, 0)] = 0.0;
  
  /* set 4 edges (excluding corners) to zero */
  for (Slice=0; Slice<NumSlices; Slice++)
    for (i=1; i<=NewDim; i++)
      NewMatrix[f_ptr(0,Slice,i)] = NewMatrix[f_ptr(NewDim+1,Slice,i)] =
	NewMatrix[f_ptr(i,Slice,0)] = NewMatrix[f_ptr(i,Slice,NewDim+1)]
	= 0.0;
#endif
  
  NewAxisX = NewAxisY = 0.5*(NewDim-1);
  OldAxisX = OldAxisY = 0.5*(OldDim-1);
  Radius = 0.5*NewDim;
  
  for (NewY=0; NewY < NewDim; NewY++) {
    Y = NewY - NewAxisY;
    /*compute head pixel and tail pixel that are going to be used for each
      row, assuming the image grid is at least 1 pixel wider than actually 
      needed on 4 sides, --> no need to do testing for out of bounds. */
    
    NumPix = (int)(sqrt(Radius*Radius-(Radius-NewY-.5)*(Radius-NewY-.5))+.5);
    NumPix *= 2;
    HeadPixel = (int)(0.5*(NewDim-NumPix));
    TailPixel = HeadPixel + NumPix - 1;
    
    /* pixels (0 to head-1, tail+1 to NewDim-1) are set to zero */
    
#ifdef FAST
    for (Slice = 0; Slice< NumSlices; Slice++)
      for (NewX=0; NewX < HeadPixel; NewX++) {
	NewMatrix[f_ptr(NewY+1, Slice, NewX+1)] = 0.0;
	NewMatrix[f_ptr(NewY+1, Slice, NewX+1+TailPixel+1)] = 0.0;
      }
#else
    for (Slice=0; Slice < NumSlices; Slice++)
      for (NewX=0; NewX < HeadPixel; NewX++) {
	NewMatrix[n_ptr(NewY, Slice, NewX)] = 0.0;
	NewMatrix[n_ptr(NewY, Slice, NewX+TailPixel+1)] = 0.0;
      }
#endif
    
    for (NewX = HeadPixel; NewX <= TailPixel; NewX++) {
      X = NewX - NewAxisX;
      
	FloatX = X*Cos - Y*Sin + OldAxisX;  /*rotation about center of image*/
	FloatY = X*Sin + Y*Cos + OldAxisY;
      OldX = (int)FloatX;
      OldY = (int)FloatY;
#ifdef FAST
      for (Slice=0; Slice < NumSlices; Slice++) {
	
	  z0 = OldMatrix[o_ptr(OldY, Slice, OldX)];
	  z1 = OldMatrix[o_ptr(OldY+1, Slice, OldX)];
	  z2 = OldMatrix[o_ptr(OldY, Slice, OldX+1)];
	  z3 = OldMatrix[o_ptr(OldY+1, Slice, OldX+1)];
	}
	
	z01 = z0 + (FloatY - OldY) * (z1 - z0);
	z23 = z2 + (FloatY - OldY) * (z3 - z2);
	z0123 = z01 + (FloatX - OldX) * (z23 - z01);
	
	NewMatrix[f_ptr(NewY+1, Slice, NewX+1)] = z0123;
      }
#else
      
      z0out = z1out = z2out = z3out = 0;
      
      if (OldX < 0 || OldX >= OldDim) {
	z0out = 1; 
	z1out = 1;
      }
      if (OldY < 0 || OldY >= OldDim) {
	z0out = 1; 
	z2out = 1;
      }
      if (OldX+1<0 || OldX+1>=OldDim) {
	z2out = 1; 
	z3out = 1;
      }
      if (OldY+1<0 || OldY+1>=OldDim) {
	z1out = 1; 
	z3out = 1;
      }
      
      for (Slice=0; Slice<NumSlices; Slice++) {
	z0 = z1 = z2 = z3 = 0.0;
	
	if (!z0out){
	   z0 = OldMatrix[o_ptr(OldY, Slice, OldX)];
	}
	if (!z1out){
	    z1 = OldMatrix[o_ptr(OldY+1, Slice, OldX)];
	}
	if (!z2out){
	    z2 = OldMatrix[o_ptr(OldY, Slice, OldX+1)];
	}
	if (!z3out){
	    z3 = OldMatrix[o_ptr(OldY+1, Slice, OldX+1)];
	}
	
	z01 = z0 + (FloatY - OldY) * (z1 - z0);
	z23 = z2 + (FloatY - OldY) * (z3 - z2);
	if(FloatY < 0){
	  z01 = z0;
	  z23 = z2;
	}  
	z0123 = z01 + (FloatX - OldX) * (z23 - z01);
	if(FloatX < 0){
	  z0123 = z01;
	}
	NewMatrix[n_ptr(NewY, Slice, NewX)] = z0123;
      }
#endif
    }
  }
#undef o_ptrj
#undef n_ptr
#ifdef FAST
#undef f_ptr
#endif
}

