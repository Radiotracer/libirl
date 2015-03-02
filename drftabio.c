/*
	drftabio.c
	
	$Id: drftabio.c 109 2010-04-22 16:15:01Z frey $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mip/miputil.h>
#include <mip/errdefs.h>
#include <mip/printmsg.h>
#include "drftabio.h"

#include "irlprivate.h"

#define DRF_TITLESIZE 80
#define COLSIZE 20
#define DRFTAB_VERSION (int)1
#define SIGNATURE "DRFT"
#define SIGSIZE 4

/* The routines vWriteDrfTab and sReadDrfTab read and write drf tables.
	The format of the files is as follows below:
		signature: char[4]:4 bytes=DRFT 
		version: int: binary value equal to 1 in this version. Note that this
			can also be used to see if the byte order is correct.
		title: char[DRF_TITLESIZE]:null terminated string containing info about
		title: char[COLSIZE]: null terminated string containing collimator name
		float fPixSize: size of pixels in a plane in cm.
		float fDrfSpacing: distance between planes in the table in cm
		float fStartDistance: starting distance for table in cm
		int iNumDistances: number of distances for which there is a DRF
		int bIncrKrnl: nonzero if the DRF is to be used as an incremental
			kernel, zero otherwise
	Following this there are iNumDstances items that contain the Drf for
 	each distance.
		int iXoffs, iYoffs: number of offsets from the center in the X and Y 
			directions. Note that the  size of the table will then be
			2*iXoffs-1 by 2*iYoffs-1.
		there are then (2*iXoffs-1)*(2*iYoffs)-1 floating point values. These
		are stored with X varying fastest. Note that, conventionally, X is the
		transverse direction (perp to axis of rotation) and Y is the axial
		direction.
*/
static vCheckIoError(FILE *fp, char *pchRoutine, char *pchFile, char *pchMsg)
{
	if (ferror(fp)){
		vErrorHandler(ECLASS_FATAL, ETYPE_IO, pchRoutine, 
			"IO Error for file %s %s", pchFile, pchMsg);
	}
}

void vWriteDrfTab(DrfTab_t *psOutDrfTab, char *pchOutFname, char *pchTitle, 
	char *pchColName)
/* writes a DRF Table in the format described above to pchOutFname. The input
	is a pointer to a DrfTab structure (defined in irl.h)
*/
{
	FILE *fpOut;
	char pchOutTitle[DRF_TITLESIZE];
	char pchOutColName[COLSIZE];
	int iLen, i;
	int iXoffs, iYoffs;
	const int iDrfTabVersion=DRFTAB_VERSION;
	char pchErrBuf[256];
	Drf_t *pDrfs;

	vPrintMsg(1,"WriteDrfTab: %s\n",pchOutFname);
	fpOut = fopen(pchOutFname,"wb");
	if (fpOut == NULL){
		vErrorHandler(ECLASS_FATAL, ETYPE_IO, "WriteDrfTab",
			"Error opening DRF Table file: %s", pchOutFname);
	}
	fwrite(SIGNATURE, sizeof(char), SIGSIZE, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", ": writing signature");
	fwrite(&iDrfTabVersion, sizeof(int),1,fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", ": writing version");

	iLen = strlen(pchTitle);
	if (iLen > DRF_TITLESIZE){
		vErrorHandler(ECLASS_WARN, ETYPE_INFO, "vWriteDrfTab",
			"DRF Table Title is too long (%d characters): truncating to %d",
			iLen-1, DRF_TITLESIZE-1);
	}
	strncpy(pchOutTitle, pchTitle, DRF_TITLESIZE);
	/* fill the rest of the string with titles. Start with byte iLen since
		iLen also includes the null terminator*/
	for(i=iLen; i < DRF_TITLESIZE; ++i)
		pchOutTitle[i] = '\0';
	fwrite(pchOutTitle, sizeof(char),DRF_TITLESIZE,fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", ": writing title");

	iLen = strlen(pchColName);
	if (iLen > COLSIZE){
		vErrorHandler(ECLASS_WARN, ETYPE_INFO, "vWriteDrfTab",
			"Collimator Name is too long (%d characters): truncating to %d",
			iLen-1, COLSIZE-1);
	}
	strncpy(pchOutColName, pchColName, COLSIZE);
	/* fill the rest of the string with titles. Start with byte iLen since
		iLen also includes the null terminator*/
	for(i=iLen; i < COLSIZE; ++i)
		pchOutColName[i] = '\0';
	fwrite(pchOutColName, sizeof(char),COLSIZE,fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", ": writing col name");

	fwrite(&(psOutDrfTab->fPixSize), sizeof(float),1, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", ": writing pixsize");
	fwrite(&(psOutDrfTab->fDrfSpacing), sizeof(float),1, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", ": writing DrfSpacing");
	fwrite(&(psOutDrfTab->fStartDistance), sizeof(float),1, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", ": writing StartDistance");

	fwrite(&(psOutDrfTab->iNumDistances), sizeof(int),1, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", ": writing NumDistances");

	fwrite(&(psOutDrfTab->bIncKrnl), sizeof(int),1, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", ": writing bIncKrnl");

	pDrfs = psOutDrfTab->psDrfs; /* this is the array of Drf structures*/
	for(i=0; i<psOutDrfTab->iNumDistances; ++i){
		iXoffs = pDrfs[i].MaxOffsX;
		iYoffs = pDrfs[i].MaxOffsY;
#ifdef DEBUG
		fprintf(stderr,"dist=%d, xoffs=%d, yoffs=%d, sum=%.4g\n",
			i, iXoffs, iYoffs, sum_float(pDrfs[i].Psf, (2*iXoffs-1)*(2*iYoffs-1)));
#endif

		fwrite(&iXoffs, sizeof(int), 1, fpOut);
		sprintf(pchErrBuf,": writing iXoffs, dist %d", i);
		vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", pchErrBuf);

		fwrite(&iYoffs, sizeof(int), 1, fpOut);
		sprintf(pchErrBuf,": writing iYoffs, dist %d", i);
		vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", pchErrBuf);

		fwrite(pDrfs[i].Psf, sizeof(float), (2*iXoffs-1)*(2*iYoffs-1), fpOut);
		sprintf(pchErrBuf,": writing PSF, dist %d", i);
		vCheckIoError(fpOut, pchOutFname, "WriteDrfTab", pchErrBuf);
	}
	/*
	fclose(fpOut);
	*/
}

DrfTab_t sReadDrfTab(char *pchInFname,  char **ppchInTitle, 
	char **ppchInColName)
/* reads a Drf table from the file pchInFname and returns a pointer to the 
	DrfTab_t structure that contains the data. Also returns pointers to
	the title and collimator name strings.
*/
{
	FILE *fpIn;
	char pchSig[SIGSIZE+1];
	char pchInTitle[DRF_TITLESIZE];
	char pchInColName[COLSIZE];
	int i;
	int iXoffs, iYoffs;
	int iPsfSize; /* total number of elements in Psf*/
	int bSwapped;/* set to true if the file is in the wrong byte order*/
	int iDrfTabVersion;
	char pchErrBuf[256];
	Drf_t *pDrfs;
	DrfTab_t sInDrfTab;

	/* initialized structure items in case reads below fail*/
	sInDrfTab.fPixSize=sInDrfTab.fDrfSpacing=sInDrfTab.fStartDistance=0.0;
	sInDrfTab.iNumDistances=0;
	sInDrfTab.psDrfs=NULL;
	sInDrfTab.bDrfFromFile=TRUE;

	fpIn = fopen(pchInFname,"rb");
	if (fpIn == NULL){
		vErrorHandler(ECLASS_FATAL, ETYPE_IO, "ReadDrfTab",
			"Error opening DRF Table file: %s", pchInFname);
	}
	fread(pchSig,sizeof(char),SIGSIZE,fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadDrfTab", ": reading signature");
	pchSig[SIGSIZE] = '\0';
	if (strncmp(SIGNATURE,pchSig,SIGSIZE)){
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "ReadDrfTab",
			"Signature in Drf File %s is %s, should be %s", pchInFname,
			pchSig, SIGNATURE);
	}

	fread(&iDrfTabVersion, sizeof(int),1,fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadDrfTab", ": reading version");
	if (iDrfTabVersion == DRFTAB_VERSION){
		bSwapped=FALSE;
		vPrintMsg(6,"Input DrfTab doesn't need byte swapping\n");
	}else{
		/* try swapping the bytes in the version to see if the file is
			from a dift. byte order */
		vSwapWords((void *)&iDrfTabVersion,1);
		if (iDrfTabVersion == DRFTAB_VERSION){
			bSwapped=TRUE;
			vPrintMsg(6,"swapping bytes in DrfTab file %s\n",pchInFname);
		}else{
			/* even after swapping the version doesn't match*/
			vSwapWords((void *)&iDrfTabVersion,1);/*swap back before reporting err*/
			vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "ReadDrfTab",
				"Version in Drf File %s is  %x hex, should be %x", pchInFname,
				iDrfTabVersion, DRFTAB_VERSION);
		}
	}
	pchInTitle[0]='\0';
	fread(pchInTitle, sizeof(char), DRF_TITLESIZE, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadDrfTab", ": reading title");
	*ppchInTitle = pchIrlStrdup(pchInTitle);

	pchInColName[0]='\0';
	fread(pchInColName, sizeof(char), COLSIZE, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadDrfTab", ": reading colname");
	*ppchInColName = pchIrlStrdup(pchInColName);
				

	fread(&(sInDrfTab.fPixSize), sizeof(float),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadDrfTab", ": reading PixSize");
	if (bSwapped) vSwapWords((void *)&(sInDrfTab.fPixSize),1);

	fread(&(sInDrfTab.fDrfSpacing), sizeof(float),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadDrfTab", ": reading DrfSpacing");
	if (bSwapped)vSwapWords((void *)&(sInDrfTab.fDrfSpacing),1);

	fread(&(sInDrfTab.fStartDistance), sizeof(float),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadDrfTab", ": reading StartDistance");
	if (bSwapped)vSwapWords((void *)&(sInDrfTab.fStartDistance),1);

	fread(&(sInDrfTab.iNumDistances), sizeof(int),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadDrfTab", ": reading NumDistances");
	if (bSwapped)vSwapWords((void *)&(sInDrfTab.iNumDistances),1);

	fread(&(sInDrfTab.bIncKrnl), sizeof(int),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadDrfTab", ": reading bIncKrnl");
	if (bSwapped)vSwapWords((void *)&(sInDrfTab.bIncKrnl),1);

	vPrintMsg(6,"ReadDrfTab:\n   Drf=%s\n   Col=%s\n",
		pchInTitle, pchInColName);
	vPrintMsg(6, 
		"   Pixsize=%.3g, DrfSpacing=%.3g, fStartDist=%.2g, NumDist=%d\n",
		sInDrfTab.fPixSize, sInDrfTab.fDrfSpacing, sInDrfTab.fStartDistance,
		sInDrfTab.iNumDistances);

	pDrfs = (Drf_t *)pvIrlMalloc(sInDrfTab.iNumDistances*sizeof(Drf_t),
						"ReadDrfTab:Drfs");
	for(i=0; i<sInDrfTab.iNumDistances; ++i){
		fread(&iXoffs, sizeof(int), 1, fpIn);
		if (bSwapped)vSwapWords((void *)&iXoffs,1);
		fread(&iYoffs, sizeof(int), 1, fpIn);
		if (bSwapped)vSwapWords((void *)&iYoffs,1);

		pDrfs[i].MaxOffsX = iXoffs;
		pDrfs[i].MaxOffsY = iYoffs;
		pDrfs[i].bFFTdone=0;
		iPsfSize=(2*iXoffs-1)*(2*iYoffs-1); 
		pDrfs[i].Psf = vector(iPsfSize, "ReadDrfTab:Drfs[i].Psf");
		fread(pDrfs[i].Psf, sizeof(float), iPsfSize, fpIn);
		sprintf(pchErrBuf,": reading drf %d",i);
		vCheckIoError(fpIn, pchInFname, "ReadDrfTab", pchErrBuf);

#ifdef DEBUG
		fprintf(stderr,"dist %d: read %d bytes, sum=%.3g\n", i,
			iPsfSize, sum_float(pDrfs[i].Psf, iPsfSize));
#endif
		if (bSwapped)vSwapWords(pDrfs[i].Psf,iPsfSize);
	}
	sInDrfTab.psDrfs=pDrfs;

	fclose(fpIn);
	return sInDrfTab;
}
