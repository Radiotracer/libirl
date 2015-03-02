/*
	irl.h

	This file declares the public API for the MIP libirl

	$Id: irl.h 113 2010-11-11 17:43:22Z frey $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#ifndef IRL_H
#define IRL_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef TRUE
#define TRUE	1
#endif

#ifndef FALSE
#define FALSE	0
#endif

	/* Masks that can be used to test iBckModel and iPrjModel to see which
		 effects are modeled */

#define MODEL_ATN 1
#define MODEL_DRF 2
#define MODEL_SRF 4

	/* general parameters that are common to all recon and
		 projection programs. They are divided into those that are input
		 by the user and those that are computed
		 user supplied parameters (via parameter file or command line
	*/

typedef struct {
	int NumPixels;         /* # of image pixels per row/col, NDIM(U)
														This must be the same as the number of
														projection bins in the transaxial direction */
	int NumSlices;
	int NumViews;			  /* number of angles in the reconstructed image */
	float BinWidth;        /* projection bin width in units of cm */
	/* this will be the same size as the pixels
		 in the reconstructed image */
	/* parameters for attenuation map */
	float fAtnScaleFac;		  /* factor to use to scale the attenuation map
														 to units of attenuation per pixel at the
														 primary photon energy*/
	int iNumSrfIterations; /* number of iterations in which to include scatter
														modeling*/
	float fScatEstFac;				  /* factor to be used to multiply ScatEst */
	
	int NumAngPerSubset;       /* Number of angles in a subset. Must be an
																integer factor*/
	int NumIterations;
	char *pchNormImageBase; /*Name of the directory and file basename into which
														The normalization images will be written */
	
	int SrfCollapseFac; /* Factor by which to collapse the reconstructed image
												 when doing scatter compensation. Factors of 1 (no
												 collapse) 2 or 4 are supported. Note that for factors
												 other than 1, the number of pixels and slices must
												 be divisible by this. This significanty speeds
												 up scatter modeling. The collapse factor should be
												 chosen so the pixel size after collapse is not
												 bigger than ~1.3 cm.
											*/
	
	/* parameters required if drf modeling using a grf table computed on-the-fly
		 will be used. (i.e., bModelDrf is true and pchDrfTabFname is NULL) */
	float fHoleDiam; /* diameter collimator hole having same area as hole
											of collimator for which GRF compensation is desired*/
	float fHoleLen;
	float fBackToDet;
	float fIntrinsicFWHM;
} IrlParms_t;
	
typedef struct {
	int bModelDrf; 			/*true if the DRF or Grf will be modeled */
	int bReconIsInitEst; 	/*nonzero if the reconstructed image contains an
													initial estimate. If so, it must have no negative
													values (they will be set to 0 otherwise).
													If false, then the initial estimate will be
													computed.
												*/
	int bUseGrfInBck;			/* if nonzero, then the GRF is computed on using
													 the parameters in IrlParms and used during
													 the backprojection. This can be used to speed
													 up reconstruction without a significant
													 sacrifice of image quality when DRF modeling
													 is used
												*/
	int bFFTConvolve;			/* If true, then FFT-based convolution is used
													 when modeling the DRF. Generally, when the
													 GRF is modeled, the kernels are so small that
													 spatial convolution is faster. However, for
													 the full DRF, FFT convolution is faster.
												*/
	int bUseContourSupport; /* defines the region where the initial estimate
														 is nonzero based on the orbit (i.e., everything
														 inside the orbit is nonzero */
	float fAtnMapThresh;		/* If positive, then everything above this value
														 in the attenuation map is assumed to be part of
														 the reconstruction. */
	int iMsgLevel; /* a number from 1 to 9 describing the level of verbosity
										for writing status information to the message file. 
										Smaller numbers result in less information. The default
										level is 4. Levels greater than 6 may result in a large
										amount of information. */
	int iAxialPadLength; /* Pads the reconstruction image, attenuation map,
													and projection data with iAxialPadLength pixels
													on both axial ends. The padding is done by averaging
													iAxialPadAverage pixels in the axial direction.
													This must be >=0. If iAxialPadLength == 0, then no
													padding is done.
											 */
	int iAxialAvgLength; /* Number of pixels in axial direction to average
													when performing axial padding. if iAxialPadLength
													is > 0, then this must be >= 1.
											 */
} Options_t;
	
	/* the following structure must be filled for each projection view */
	typedef struct {
		float Angle;            /* Projection angle relative to the x-axis in 
															 radians. Note that when reconstructed image
															 slices are viewed with pixel 0,0 in the upper
															 left hand corner, then the rotation direction
															 must be negative.
														*/
		float Left;					/* location in projections of start of
													 measured projection bins. This can be used
													 to take into account the dead area and reduce
													 truncation artifacts that result when this
													 area is treated as zero rather than unmeasured */
		float Right;				/* location in modified projections where end of
													 measured projection bin projects.*/
		float CFCR;             /* distance from the CoR to collimator face in (cm)
														 */
	} PrjView_t;

	int IrlOsem
	(
	 IrlParms_t *psIrlParms,	 /* reconstruction parameters, see below */
	 Options_t *psOptions,	 /*structure containing reconstruction options */
	 PrjView_t *psViews, /*array of View_t structures, one for each projection
												 angle (see below) in the same order as the actual
												 images
											 */
	 char *pchDrfTabFname, /*filename containing the drf table. If NULL, and
													 drf modeling is requested, then the grf is
													 computed on the fly using the parameters contained
													 in the IrlParms structure (fHoleDiam, fHoleLen,
													 fBackToDet, and fIntrinsicFwhm). Note that the
													 DRF table file must have the same pixel size as
													 the reconstructed image.
												 */
									
	 char *pchScatKrnlFname, /*filename containing the scatter kernel file.
														 if NULL, then scatter is not modeled.
														 The pixel. The scatter kernel should be for the
														 appropriate isotope and must have the same pixel
														 size as the reconstructed image after factoring
														 in the SrfCollapseFac. For example, if the
														 pixel size in the  reconstructed image is 0.6 cm,
														 and a collapse factor of 2 is used, then the
														 pixel size in the scatter kernel must be 1.2 cm 
													 */
	 void (*pIterationCallbackFunc)(int IterationNumber, float *pfCurrentEst),
	 /*callback function that is called at the end of every iteration.
		 The callback function has two arguments: the iteration number
		 and a pointer to the Current Estimate image. If this is Null,
		 then no callback will be made */
	 float *pfScatterEstimate, /* Estimate of scatter that will be added to
																the reconstructed image estimate after
																the computed projection. Must be the same
																size and pixel order as the
																projection image. If NULL, then no scatter
																estimate is used.
																Notes: 1. The estimate is scaled by fScatEstFac
																before being added.
																2. This is indendent of scatter modeling
																and can be used in addition or in
																lieu of it.
														 */
	 float *pfAtnMap,				/* pointer to attenuation map image. This
														 must be the same size as the reconstructed
														 image (NumPixels*NumPixels*NumSlices).
														 The pixel ordering is (fastest to slowest)
														 x, y, z (slice). After scaling by AtnScaleFac, 
														 the units should be in per
														 pixel at the energy of the primary photons.
													*/
	 float *pfPrjImage,		 	/* pointer to projection image. This must have
														 the size NumPixels*NumPixels*NumSlices. The
														 pixel ordering is (fastest to slowest):
														 transaxial bin, slice, view
													*/
	 float *pfReconImage,			/* pointer to reconstructed image initial
															 estimate (if bReconIsInitEst is true) and
															 memory where the reconstructed image will
															 be stored.
														*/
										  
	 char *pchLogFname,				/* file name where warning, information, and
															 fatal error messages are written. If NULL,
															 this will not be created */
	 char *pchMessageFname		/* file name where status messages are printed. 
															 if Null, these are not printed */
	 );

	int IrlGenprj(
								IrlParms_t *psIrlParms, /* reconstruction parameters, see below */
								Options_t *psOptions, /*structure containing reconstruction options */
								PrjView_t *psPrjViews, /*array of View_t structures, one for each projection
																				 angle (see below) in the same order as the actual
																				 images */
								char *pchDrfTabFname, /*filename containing the drf table. If NULL, and
																				drf modeling is requested, then the grf is
																				computed on the fly using the parameters contained
																				in the IrlParms structure (fHoleDiam, fHoleLen,
																				fBackToDet, and fIntrinsicFwhm). Note that the
																				DRF table file must have the same pixel size as
																				the reconstructed image.
																			*/

								char *pchScatKrnlFname, /*filename containing the scatter kernel file.
																					if NULL, then scatter is not modeled.
																					The pixel. The scatter kernel should be for the
																					appropriate isotope and must have the same pixel
																					size as the reconstructed image after factoring
																					in the SrfCollapseFac. For example, if the
																					pixel size in the  reconstructed image is 0.6 cm,
																					and a collapse factor of 2 is used, then the
																					pixel size in the scatter kernel must be 1.2 cm
																				*/
								void (*pPrjViewCallbackFunc)(int PrjViewNumber, float *pfCurrentPrjView),
								/*callback function that is called at the end of each
									projection view. The callback function has two arguments:
									the view number and a pointer to the Current projection
									view image. If this is Null, then no callback will be made
								*/
								float *pfScatterEstimate, /* Estimate of scatter that will be added to
																						 the the computed projection. Must be the same
																						 size and pixel order as the
																						 projection image. If NULL, then no scatter
																						 estimate is used.
																						 Notes: 1. The estimate is scaled by fScatEstFac
																						 before being added.
																						 2. This is indendent of scatter modeling
																						 and can be used in addition or in
																						 lieu of it.
																						 3. Currently only supports same bin size as
																						 pixel size. The interpolation from different
																						 bin size is expected to be done by user
																						 before calling irlGenprj
																					*/
								float *pfActImage,         /* pointer to image of activity distribution
																						*/
								float *pfAtnMap,           /* pointer to attenuation map image. This
																							must be the same size as the reconstructed
																							image (NumPixels*NumPixels*NumSlices).
																							The pixel ordering is (fastest to slowest)
																							x, y, z (slice).
																					 */
								float *pfPrjImage,         /* pointer to projection image. This must have
																							the size NumPixels*NumPixels*NumSlices. The
																							pixel ordering is (fastest to slowest):
																							transaxial bin, slice, view. Currently only
																							supports using of bin-size equals pixel-size.
																							For different bin and pixel size, the interpolation
																							is expected to be done by user after calling irlGenprj.
																							It is assumed that the activity image represents the total
																							number of decays times the sensitivity of the collimator.
																							As a result, each projection review represents
																							1/NumAngles times the number of decays.
																					 */
								float fPrimaryFac,			// Primary factor, use 0 when generate scatter estimation only
								char *pchLogFname,         /* file name where warning, information, and
																							fatal error messages are written. If NULL,
																							this will not be created */
								char *pchMessageFname      /* file name where status messages are printed.
																							if Null, these are not printed */

								);


	int ComputeOrbit(float *pfPixels, /* 2D image used to compute orbit */
		IrlParms_t *psParms,
		float fMargin, /* Gap between edge of patient and camera face
											that is used in computing the CFCR for each view
											Ignored if pfMargins is non-null
											For circular orbits, this is treated as the CFCR. */
		float *pfMargins, /* If non-null, then this is a vector of lenght NumViews
												 which contains the value of fMargin for each view.
												 If null, then it is ignored and the value fMargin is used
												 for all views */
		float fThresh, 	/* Threshold specifying used in computing edge of object
											 If positive, then it is the percent of maximum
											 If negative, then it is an absolute threshold */
	
		int OrbitType, /* 0: circular orbit. fMargin is treated as the CFCR for
											each view is placed in the CFCR element for each view
											1: noncircular orbit with camera at each view moved 
											independently so it is fMargin from the edge of the
											patient as defined by the pixel values in fImage */
		PrjView_t *psViews /* Array of view structures. The CFCR is computed
													for each view and put in the CFCR element for the
													view. These must already have the projection angle
													placed in the view structure */
									 ); /* returns 0 if there is an error in orbit computation, 
												 1 if there was a memory allocation error, and 2 if there was
												 an illegal orbit type */
	
	char *pchIrlErrorString(void); /* returns error message for last fatal IRL error */

#ifdef __cplusplus
}
#endif

#endif /* IRL_H */
