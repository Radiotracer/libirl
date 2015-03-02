/*
	irlprivate.h

	$Id: irlprivate.h 127 2011-10-05 20:18:28Z frey $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#ifndef IRLPRIVATE_H
#define IRLPRIVATE_H

#include "irl.h" // This is the header for the public API
#include "drftabio.h" // Defines DRF Type

#ifndef TRUE
#define TRUE	1
#endif

#ifndef FALSE
#define FALSE	0
#endif

/* note that only parallel is supported in S.1, but the other constants are
	defined because they are used in drf.h
*/
#define PARALLEL           0
#define CONE               1
#define FAN                2
#define HFAN					3

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif
#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

#define MAX_COLLAPSE_FAC	4

#define SAVE_PIXEL_FORMAT_FLOAT 1
#define SAVE_PIXEL_FORMAT_SHORT 2

typedef struct {
	/* general parameters that are common to all recon and
		projection programs. They are divided into those that are input
		by the user and those that are computed*/
	/* user supplied parameters (via parameter file or command line */
	int Circular;          /* true for circular reconstruction */
	int NumPixels;         /* # of image pixels per row/col, NDIM(U) */
	int SliceStart, SliceEnd; /* starting and ending slices for reconstruction*/
	/* note that in the current version PixelWidth and BinWidth and 
			SliceThickness should probably be the same. If they aren't then
			the projection data is collapsed or interpolated to make them
			the same */
	float PixelWidth;      /* image pixel width in units of cm */
	float BinWidth;        /* projection bin width in units of cm */
	float SliceThickness;  /* image slice thickness in units of cm */

	/* parameters needed to fill views structures */
	/* parameters for attenuation map */
	float fAtnScaleFac;		  /* factor to use to scale the attenuation map
										  to units of attenuation per pixel at the
										  primary photon energy*/
	/* parameters for scatter modeling */
	int iNumSrfIterations; /* Number of iterations to do SRF Modeling */
	float fScatEstFac;				  /* factor to be used to multiply ScatEst */
	int iPrjModel; /* effects modeled in projector */
	int iBckModel; /* effects modeled in backprojector */
	int iModel; /* this is computed as the bitwise or of the two above */

	/* parameters that are required for reconstruction only */
	char *pchNormImageBase; /* if NULL, the normalization matrix is stored
										in memory. Otherwise, it is stored on disk.
									*/
	int NumAngPerSubset;       /* Number of angles in a subset*/
	int NumIterations;
	int StartIteration; /*starting iteration number (designed for restarting
								using initial estimate obtained from a previous
								run */
	int SaveInterval; /* distance between iterations to save. Used only if
									IterSaveString isn't specified*/
	char *pchIterSaveString; /* string specifying which iterations to change
								the save interval on and what the save interval
								should be */

	/* parameters that are determined from images (if .im is used) or user
		input (if .im is not used) */
	int NumBins;           /* # of projection bins per view, KDIMU in RECLBL*/
	int NumAngles;         /* usually = ZDim of the projection image */
		
	/* parameters for GRF modeling */
	float fHoleDiam;
	float fHoleLen;
	float fBackToDet;
	float fIntrinsicFWHM;
	float fMaxFracErr;
	/* parameters that are used to reduce edge effects at bottom and top of
		reconstructed stack of images */
	int iAxialPadLength, iAxialAvgLength;
	/* Parameters that are computed from user inputs, ...*/
	int  NumRotPixs;        /* Side length of transaxial images used in
										  rotation, projection and backprojection */
	int NumSlices;         /* # of slices to be reconstructed or projected*/
} Parms_t;

/* Masks that can be used to test iBckModel and iPrjModel to see which
	effects are modeled */
#define MODEL_ATN 1
#define MODEL_DRF 2
#define MODEL_SRF 4


typedef struct {
	int iSrfCollapseFac;
	int iSrfParmLineNum;
	void *pvScatterData;
} Srf_t;

typedef struct {
	float fHoleDiam;
	float fHoleLen;
	float fBackToDet;
	float fIntrinsicFWHM;
	float fSensitivity;
	int bFFTconvolve;
	int bUseGrfInBck;
	char *pchDrfTabFname;
	DrfTab_t *psDrfTab;
} Detector_t;

typedef struct {
	int iAngle; /* projection number */
	float Angle;            /* relative to the x-axis. radian/degree */
	float Sin;              /* sin(Angle) */
	float Cos;              /* cos(Angle) */
	float Center;           /* of the projection of the CoR, AXISU */
	float Left;					/* location in modified projections where start of
										measured projection bin projects */
	float Right;				/* location in modified projections where end of
										measured projection bin projects.*/
	float CFCR;             /* of the CoR to collimator face, CFCR (cm) */
	float ViewTime;			/* Time per projection view. This is used
										to scale the projections after projection,
										so the units affect the reconstructed image
										and must be consistent for all projection
										views */
	float Sensitivity;		/* Collimator Sensitivity. Units will affect
										units of reconstructed image */
	DrfTab_t *psDrfTab;
	
} View_t;

typedef struct {
	int bNormInMemory;/* true if we are storing norm matrices in memory*/
	int bSaveNormImages; /* true if we are saving the norm matrices */
	int iNumSlices, iNumPix;/* duplicated from psParms*/
	char *pchNormImageBase;/* base name for file to store norm images*/
	char *pchNormImageFileName;/* actual file name that is used to open
											disk versions of norm images*/
	float *pfNormImages; /*location of all images if we are storing them
									in memory*/
} NormData_t;


/* PROTOTYPES FOR INTERNAL FUNCTIONS */

/* computeorbit.c */

int ComputeOrbit(
	float *pfPixels,
	IrlParms_t *psParms,
	float fMargin,
	float *pfMargins,
	float fThresh,
	int OrbitType,
	PrjView_t *psViews
	);

/* conv2d.c */

void conv2d(
	float *a,
	int aSizeX,
	int aSizeY,
	float *b,
	int bCtrX,
	int bCtrY,
	float *out
	);

/* conv2dPsfFFT.c */

void conv2dPsfFFT(
	float *pfIn,
	int iInXdim,
	int iInYdim,
	float *pfPsfFT,
	int iPsfXdim,
	int iPsfYdim,
	float *pfOut
	);

/* conv2d_ypad.c */

void conv2d_ypad(
	float *a,
	int aSizeX,
	int aSizeY,
	float *b,
	int bCtrX,
	int bCtrY,
	int iPadY,
	int iAvgY,
	float *out
	);

/* convolve.c */

void ConvolveDrf(
	float *pfIn,
	int inX,
	int inY,
	Drf_t *psDrf,
	int bFFTConvolve,
	float *pfOut
	);

/* convolve_ypad.c */

void ConvolveDrfYpad(
	float *pfIn,
	int inX,
	int inY,
	Drf_t *psDrf,
	int bFFTConvolve,
	int iYpad,
	int iYavg,
	float *pfOut
	);

/* fastrotate.c */

void vShearsRotateCopy(
	float fAngle,
	float *pfOldMatrix,
	float *pfRotMatrix,
	int PB_FLAG,
	int NumSlices,
	int NumPixels
	);

void vShearsRotate(
	float fAngle,
	float *pfOldMatrix,
	float *pfRotMatrix,
	int PB_FLAG,
	int NumSlices,
	int NumPixels
	);

void vInitFastRotate(
  int NumPixels
  );
 
void vFreeFastRotate(void);

/* gendrf.c */

void gendrf(
	int Type,
	float Rfan,
	float fColThick,
	float fHoleDiam,
	float fGap,
	float fIntrinsic,
	float fMaxErr,
	float fBinSize,
	float fDrfSpacing,
	float minDist,
	int NumDists,
	int NOffs,
	Drf_t *DrfTab
	);

DrfTab_t *psFreeDrf(
	DrfTab_t *psDrfTab
	);

/* initial_est.c */

float *pfComputeInitialEstimate(
	Parms_t *psParms,
	int bUseContourSupport,
	float fAtnMapThresh,
	View_t *Views,
	float *pfPrjImage,
	float *pfAtnMap,
	float *pfInitImage
	);

/* IrlSetup.c */

void vSetupParms(
	IrlParms_t *psIrlParms,
	Options_t *psOptions,
	float *pfPrjImage,
	float *pfAtnMap,
	float *pfReconImage,
	float *pfScatterEstimate,
	char *pchScatKrnlFname,
	Parms_t *psParms
	);

void FreeParms(
	Parms_t *psParms
	);

View_t *psSetupViews(
	IrlParms_t *psIrlParms,
	Options_t *psOptions,
	PrjView_t *psPrjViews
	);

Srf_t *psSrfSetup(
	Parms_t *psParms,
	int iSrfCollapseFac,
	char *pchScatKrnlFname
	);

Srf_t *psFreeSrf(
	Parms_t *psParms,
	Srf_t *psSrf
	);



DrfTab_t *psDrfSetup(float fMinCFCR, float fMaxCFCR, float fPixelWidth, 
							float fBinWidth, int iReconMatSize, float fMaxFracErr, 
							Detector_t *psDetectorParms);

void vSetupDrfTables(Parms_t *psParms, Options_t *psOptions, 
					View_t *psViews, char *pchDrfTabFname);

void vFreeExtraPrjViewData(void);

void vFreeDetectors(void);

/* krnlprj3d.c */

float *pfInterpPsf(
	float *pfRadialPsf,
	int iInNx,
	int iInNy,
	int iOutNx,
	int iOutNy,
	int iOutNz
	);

float *pfSwapQuadrants(
	float *pfKrnl,
	int nx,
	int ny
	);

void *pvDeinitScatterPrj(
	void *pvScatterData
	);

void CalcEffScatterSource(
	void *pvScatterData,
	float *pfAct,
	float *pfMap,
	int nx,
	int ny,
	int nz,
	float *pfEffSrc
	);

void *pvInitESSEScatterModel(
	char *pchKrnlFile,
	char *pchParString,
	float fPixsize,
	int nPix,
	int nSlices,
	int iCollapseFac,
	int bDebugFlag
	);

/* osemrecon.c */

void vOsemRecon(
	Parms_t *psParms,
	View_t *psViews,
	Srf_t *psSrf,
	float *pfPrjImage,
	float *pfScatterEstimate,
	float *pfAtnMap,
	float *pfReconImage,
	void (*pIterationCallbackFunc)(int IterationNumber, float *pfCurrentEst)
	);

/* pad2d.c */

float *pfPad2d(
	float *pfInPix,
	int iPadX,
	int iPadY,
	int iXdim,
	int iYdim,
	float *pfOutPix
	);

float *pfUnPad2d(
	float *pfInPix,
	int iPadX,
	int iPadY,
	int iOutXdim,
	int iOutYdim,
	float *pfOutPix
	);

float *pfAvgRows(
	float *pfInPix,
	int iXdim,
	int iNRows,
	float *pfOutPix
	);

float *pfCopyRows(
	float *pfInPix,
	int iXdim,
	int iNRows,
	float *pfOutPix
	);

/* padswap.c */

void PadSwap(
	float *pfIn,
	int inXdim,
	int inYdim,
	float *pfOut,
	int outXdim,
	int outYdim
	);

void vCopy2dArray(
	float *pfInPix,
	int nx,
	int ny,
	float *pfOutPix,
	int outnx,
	int outny
	);

/* prjbck.c */

void InitPrj(
	int iNumSlices,
	int iNumPixels,
	int iNumRotPixs,
	int iModel
	);

void Project(
	int iNumSlices,
	int iNumPixels,
	int iNumRotPixs,
	int iPrjModel,
	float fPixelWidth,
	int iAngle,
	float fAngle,
	float fCFCR,
	DrfTab_t *psDrfTab,
	Srf_t *psSrf,
	float *atnmatrix,
	float *reconmatrix,
	float *calcprj,
	float fPrimaryFac,
	float fAtnFac,
	int iAxialPadLength,
	int iAxialAvgLength,
	int bAtnMapRotated,
	int bActImgRotated,
	int bFastRotate
	);

void BackProject(
	int iNumSlices,
	int iNumPixels,
	int iNumRotPixs,
	int iBckModel,
	float fPixelWidth,
	int iAngle,
	float fAngle,
	float fCFCR,
	DrfTab_t *psDrfTab,
	Srf_t *psSrf,
	float *atnmatrix,
	float *ratioprj,
	float *imagematrix,
	float fPrimaryFac,
	float fAtnFac,
	int iAxialPadLength,
	int iAxialAvgLength,
	int bAtnMapRotated,
	int bFastRotate
	);

float *pfgConvMatrix;

int igConvMatrixSize;

void FreePrj(void);

/* rot3dpar.c */

void rot3dpar(
	float fAngle,
	float *OldMatrix,
	float *NewMatrix,
	int PB_FLAG,
	int NumSlices,
	int NumPixels,
	int NumRotPixs
	);

#endif
