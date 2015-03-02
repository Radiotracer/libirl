#ifndef DRFTABIO_H
#define DRFTABIO_H

/* drftabio.c */

/* structure containing DRF for a single distance */
typedef struct {
	int MaxOffsX;           /* Psf has circular symmetry */
	int MaxOffsY;
	int bFFTdone;
	float *Psf;             /* function of x, y, offsets */
} Drf_t;

/* structure used to store DRF table data. Used both for tables stored in
	files and those generated on the fly by gendrf. 
*/
typedef struct DrfTab_s {
	int iNumDistances;
	float fStartDistance;
	float fDrfSpacing;
	float fPixSize;
	int bIncKrnl;
	int bDrfFromFile;
	struct DrfTab_s *psBckDrfTab; /* drf table to be used for backprojection */
	float fGap, fHoleDiam, fColThick, fIntrinsicFWHM, fMaxFracErr;
	Drf_t *psDrfs; /* pointer to DrfTable data*/
	int bFFTconvolve; /* True if drf convolution is to be done using FFT */
}  DrfTab_t;

/* standard routines used to read and write DRF tables */
void vWriteDrfTab(
	DrfTab_t *psOutDrfTab,
	char *pchOutFname,
	char *pchTitle,
	char *pchColName
	);

DrfTab_t sReadDrfTab(
	char *pchInFname,
	char **ppchInTitle,
	char **ppchInColName
	);
#endif
