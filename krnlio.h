/*
	krnlio.h

	$Id: krnlio.h 30 2005-06-17 15:30:03Z mjs $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

/* structure for storing kernels */
typedef struct {
		char *pchTitle; /* string containing information about the kernel */
		char *pchIsotope; /* string containing isotope for which kernel was
									generated */
		int iNumXYoffs, iNumDists; /* xy, and z dimensions, respectively */
		float fKrnlMu0; /* atten coef for kernel */
		float fPixSize; /* pixel size (in cm) for kernel */
		float fHighEnergy, fLowEnergy; /* energy window for which kernel
													 was computed (in keV) */
		float fPriFac, /* primary photon image scaled by this before adding
								to scatter source */
				fScatFac, /* effective source scaled by this before adding
								 to primary photon image */
				fCntrPixFac; /* if nonzero, add primary photon image scaled
									by this value to image after convolution
									with scatter kernel. This is to allow improved
									coarse grid modeling.  Usually zero.*/
		float *pfKrnl, *pfKrnlMu; /* NumXYoffs by NumXYoffs by NumDists images
											 containing kernel and average atten
											 coef kernel */
} KrnlData_t;
/* prototypes*/
void vWriteKrnlFile(KrnlData_t *psKrnl, char *pchKrnlFname);
KrnlData_t *psReadKrnlFile(char *pchKrnlFname, int bReadArrays);
KrnlData_t *psReadKrnlIndirect(char *pchParFile, int iParLine);
KrnlData_t *psFreeKrnlData(KrnlData_t *psKrnl,int bFreeArrays);
