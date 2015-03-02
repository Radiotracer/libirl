/*
	krnlio.c

	$Id: krnlio.c 58 2006-01-30 21:06:01Z mjs $
 
 	Copyright 2002-2005 The Johns Hopkins University. ALL RIGHTS RESERVED.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mip/miputil.h>
#include <mip/errdefs.h>
#include <mip/printmsg.h>

#include "krnlio.h"

#define KRNL_TITLESIZE 80
#define ISOTOPE_SIZE 20
#define KRNLFILE_VERSION (int)1
#define SIGNATURE "KRNL"
#define SIGSIZE 4
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* The routines vReadKrnls and vWriteKrnls
	The format of the files is as follows below:
		signature: char[4]:4 bytes=KRNL
		version: int: binary value equal to 1 in this version. Note that this
			can also be used to see if the byte order is correct.
		title: char[KRNL_TITLESIZE]:null terminated string containing info about
		title: char[ISOTOPE_SIZE]: null terminated string containing collimator name
		float fKrnlMu0: attenuation coef for the energy window
		float fHighEnergy, fLowEnergy: energy window
		float fPixSize: size of pixels in the kernel
		float fPriFac: factor to scale primary photons by when adding to
							effective scatter source. This is used especially
							in the case in non-photopeak energy windows where
							it would a small number (often zero).
		float fScatFac: number to scale effective source by when adding to
							primary photon image. Usually set to 1.
		float fCntrPixFac: After doing the convolutions to form the effective
							scatter source but before recombining the, the
							primary source is added in after being multiplied
							by this factor.  This can be used when
							a coarse scatter grid is used to allow "convolving"
							the central pixel in the kernel by a different factor.
							Currently this is not used and thus is usually set
							to zero. Note that if this is used then this value
							must be subtracted from the central pixel of the
							scatter kernel when it is used.
		int iNumXYoffs: number of offsets parallel to collimato
		int iNumDists: number of offsets perpendicular to 

	Following this there are two 3d floating point kernels with sizes
		iNumXYoffs*iNumXYoffs*iNumDists. They are stored with X, Y, and
		   Z (distance parallel to detector) varying from fastest to slowest, 
			respectively.
	Data types: (defined in krnlio.h)
		KrnlData_t: stores kernels and information about them
	Subroutines:
		void vWriteKrnlFile(KrnlData_t *psKrnl,pchKrnlFname): Writes kernel and 
			associated data from psKrnl to file whose name is in pchKrnlFname
			
		KrnlData_t *psReadKrnlFile(char *pchKrnlFname,int bReadArrays):
			Reads kernel and 
			associated data from file whos name is in pchKrnlFname and returns 
			an pointer to a KrnlData_t structure. If bReadArrays is true the
			kernel and kernel mu arrays are read in and stored in the
			structure.
			

		vFreeKrnlData(KrnlData_t *psKrnl,int bFreeArrays);
			Frees the structure psKrnl and memory used by items it points. If
			bFreeArrays is TRUE, then the memory used by the arrays pfKrnl
			and pfKrnlMu is also freed.
*/

static void vCheckIoError(FILE *fp, char *pchRoutine, char *pchFile, char *pchMsg)
{
	if (ferror(fp)){
		vErrorHandler(ECLASS_FATAL, ETYPE_IO, pchRoutine, 
			"IO Error for file %s %s", pchFile, pchMsg);
	}
}

void vWriteKrnlFile(KrnlData_t *psKrnl, char *pchOutFname)
/* writes a Kernel File in the format described above to pchOutFname. The input
	is a pointer to a KrnlData structure (defined in irl.h)
*/
{
	FILE *fpOut;
	char pchOutTitle[KRNL_TITLESIZE];
	char pchOutIsotope[ISOTOPE_SIZE];
	char *pchTitle=psKrnl->pchTitle;
	char *pchIsotope=psKrnl->pchIsotope;
	int iLen, i,iKrnlSize;
	const int iKrnlFileVersion=KRNLFILE_VERSION;

	vPrintMsg(1,"WriteKrnlFile: %s\n",pchOutFname);
	fpOut = fopen(pchOutFname,"wb");
	if (fpOut == NULL){
		vErrorHandler(ECLASS_FATAL, ETYPE_IO, "WriteKrnlFile",
			"Error opening Krnl file: %s", pchOutFname);
	}
	fwrite(SIGNATURE, sizeof(char), SIGSIZE, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing signature");
	fwrite(&iKrnlFileVersion, sizeof(int),1,fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing version");

	iLen = strlen(pchTitle);
	if (iLen >= KRNL_TITLESIZE){
		vErrorHandler(ECLASS_WARN, ETYPE_INFO, "vWriteKrnlFile",
			"Kernel Title is too long (%d characters): truncating to %d",
			iLen-1, KRNL_TITLESIZE-1);
		iLen = KRNL_TITLESIZE - 1;
	}
	strncpy(pchOutTitle, pchTitle, KRNL_TITLESIZE);
	/* fill the rest of the string with titles. Start with byte iLen since
		iLen does not include the null terminator*/
	for(i=iLen; i < KRNL_TITLESIZE; ++i)
		pchOutTitle[i] = '\0';
	fwrite(pchOutTitle, sizeof(char),KRNL_TITLESIZE,fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnlTab", ": writing title");

	iLen = strlen(pchIsotope);
	if (iLen >= ISOTOPE_SIZE){
		vErrorHandler(ECLASS_WARN, ETYPE_INFO, "vWriteKrnlFile",
			"Isotope string is too long (%d characters): truncating to %d",
			iLen-1, ISOTOPE_SIZE-1);
		iLen = ISOTOPE_SIZE - 1;
	}
	strncpy(pchOutIsotope, psKrnl->pchIsotope, ISOTOPE_SIZE);
	/* fill the rest of the string with titles. Start with byte iLen since
		iLen does not include the null terminator*/
	for(i=iLen; i < ISOTOPE_SIZE; ++i)
		pchOutIsotope[i] = '\0';
	fwrite(pchOutIsotope, sizeof(char),ISOTOPE_SIZE,fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing isotope string");

	fwrite(&(psKrnl->fKrnlMu0), sizeof(float),1, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing KrnlMu0");

	fwrite(&(psKrnl->fHighEnergy), sizeof(float),1, fpOut);
	fwrite(&(psKrnl->fLowEnergy), sizeof(float),1, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing energy window");

	fwrite(&(psKrnl->fPixSize), sizeof(float),1, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing pixsize");

	fwrite(&(psKrnl->fPriFac), sizeof(float),1, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing primary factor");

	fwrite(&(psKrnl->fScatFac), sizeof(float),1, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing scatter factor");

	fwrite(&(psKrnl->fCntrPixFac), sizeof(float),1,fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing center pixel factor");

	fwrite(&(psKrnl->iNumXYoffs), sizeof(int),1, fpOut);
	fwrite(&(psKrnl->iNumDists), sizeof(int),1, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing kernel size");

	/* now write the two kernels */
	iKrnlSize=psKrnl->iNumXYoffs * psKrnl->iNumXYoffs * psKrnl->iNumDists;
	fwrite(psKrnl->pfKrnl, sizeof(float), iKrnlSize, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing kernel");

	fwrite(psKrnl->pfKrnlMu, sizeof(float), iKrnlSize, fpOut);
	vCheckIoError(fpOut, pchOutFname, "WriteKrnl", ": writing kernel mu");

	fclose(fpOut);
}

static void vPrintKrnlInfo(KrnlData_t *psKrnl)
{
	vPrintMsg(6,"ReadKrnl:\n   Krnl=%s\n   Isotope=%s\n",
		psKrnl->pchTitle, psKrnl->pchIsotope);
	vPrintMsg(6, 
		"   Pixsize=%.3g, KrnlMu0=%.3g, LowEnergy=%.1f, HighEnergy=%.1f\n",
		psKrnl->fPixSize, psKrnl->fKrnlMu0, psKrnl->fLowEnergy,
		psKrnl->fHighEnergy);
	vPrintMsg(6, 
		"   PriFac=%.3g, ScatFac=%.3g, CntrPixFac=%.3g\n",
		psKrnl->fPriFac, psKrnl->fScatFac, psKrnl->fCntrPixFac);
	vPrintMsg(6, 
		"   NumXYoffs=%d, NumDists=%d\n",
		psKrnl->iNumXYoffs, psKrnl->iNumDists);
}

KrnlData_t *psReadKrnlFile(char *pchInFname, int bReadArrays)
/* reads the kernels and associated information from the file pchInFname and 
   returns a pointer to the KrnlData_t structure that contains the data
*/
{
	FILE *fpIn;
	char pchSig[SIGSIZE+1];
	char pchInTitle[KRNL_TITLESIZE+1];
	char pchInIsotope[ISOTOPE_SIZE+1];
	int bSwapped;/* set to true if the file is in the wrong byte order*/
	int iKrnlVersion;
	int iKrnlSize;
	KrnlData_t *psKrnl=NULL;

	psKrnl = pvIrlMalloc(sizeof(KrnlData_t),"ReadKrnl:KrnlData");
	/* initialized structure items in case reads below fail*/
	psKrnl->fPixSize = 0.0;
	psKrnl->pchTitle = psKrnl->pchIsotope=NULL;
	psKrnl->iNumDists = psKrnl->iNumXYoffs=0;
	psKrnl->fKrnlMu0=psKrnl->fHighEnergy=psKrnl->fLowEnergy=0.0;
	psKrnl->fKrnlMu0 = psKrnl->fPriFac = psKrnl->fScatFac = 0.0;
	psKrnl->fCntrPixFac = 0.0;
	psKrnl->pfKrnl=NULL;
	psKrnl->pfKrnlMu=NULL;

	fpIn = fopen(pchInFname,"rb");
	if (fpIn == NULL){
		vErrorHandler(ECLASS_FATAL, ETYPE_IO, "ReadKrnl",
			"Error opening Krnl Table file: %s", pchInFname);
	}
	fread(pchSig,sizeof(char),SIGSIZE,fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading signature");
	pchSig[SIGSIZE] = '\0';
	if (strncmp(pchSig,SIGNATURE,SIGSIZE)){
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "ReadKrnl",
			"Signature in Krnl File %s is %s, should be %s", pchInFname,
			pchSig, SIGNATURE);
	}

	fread(&iKrnlVersion, sizeof(int),1,fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading version");
	if (iKrnlVersion == KRNLFILE_VERSION){
		bSwapped=FALSE;
		vPrintMsg(6,"Input Krnl doesn't need byte swapping\n");
	}else{
		/* try swapping the bytes in the version to see if the file is
			from a dift. byte order */
		vSwapWords((void *)&iKrnlVersion,1);
		if (iKrnlVersion == KRNLFILE_VERSION){
			bSwapped=TRUE;
			vPrintMsg(6,"swapping bytes in Krnl file %s\n",pchInFname);
		}else{
			/* even after swapping the version doesn't match*/
			vSwapWords((void *)&iKrnlVersion,1);/*swap back before reporting err*/
			vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "ReadKrnl",
				"Version in Krnl File %s is  %x hex, should be %x", pchInFname,
				iKrnlVersion, KRNLFILE_VERSION);
		}
	}
	pchInTitle[KRNL_TITLESIZE]='\0';
	fread(pchInTitle, sizeof(char), KRNL_TITLESIZE, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading title");
	psKrnl->pchTitle = pchIrlStrdup(pchInTitle);

	pchInIsotope[ISOTOPE_SIZE]='\0';
	fread(pchInIsotope, sizeof(char), ISOTOPE_SIZE, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading isotope name");
	psKrnl->pchIsotope = pchIrlStrdup(pchInIsotope);
				
	fread(&(psKrnl->fKrnlMu0), sizeof(float),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading KrnlMu0");
	if (bSwapped) vSwapWords((void *)&(psKrnl->fKrnlMu0),1);

	fread(&(psKrnl->fHighEnergy), sizeof(float),1, fpIn);
	fread(&(psKrnl->fLowEnergy), sizeof(float),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading Energy Window");
	if (bSwapped) vSwapWords((void *)&(psKrnl->fHighEnergy),1);
	if (bSwapped) vSwapWords((void *)&(psKrnl->fLowEnergy),1);

	fread(&(psKrnl->fPixSize), sizeof(float),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading PixSize");
	if (bSwapped) vSwapWords((void *)&(psKrnl->fPixSize),1);

	fread(&(psKrnl->fPriFac), sizeof(float),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading PriFac");
	if (bSwapped)vSwapWords((void *)&(psKrnl->fPriFac),1);

	fread(&(psKrnl->fScatFac), sizeof(float),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading ScatFac");
	if (bSwapped)vSwapWords((void *)&(psKrnl->fScatFac),1);

	fread(&(psKrnl->fCntrPixFac), sizeof(float),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading CntrPixFac");
	if (bSwapped)vSwapWords((void *)&(psKrnl->fCntrPixFac),1);

	fread(&(psKrnl->iNumXYoffs), sizeof(int),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading NumXYoffs");
	if (bSwapped)vSwapWords((void *)&(psKrnl->iNumDists),1);

	fread(&(psKrnl->iNumDists), sizeof(int),1, fpIn);
	vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading NumDistances");
	if (bSwapped)vSwapWords((void *)&(psKrnl->iNumDists),1);

	if (psKrnl->iNumDists <= 0 || psKrnl->iNumXYoffs < 0){
			vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "ReadKrnl",
					"sizes of kernel must be > 0 not numdists=%d, numXYoffs=%d",
					psKrnl->iNumDists, psKrnl->iNumXYoffs);
	}


	if (bReadArrays){
		iKrnlSize = psKrnl->iNumDists * psKrnl->iNumXYoffs * psKrnl->iNumXYoffs;
		psKrnl->pfKrnl = (float *)pvIrlMalloc(iKrnlSize*sizeof(float), 
				"ReadKrnl:Krnl");
		psKrnl->pfKrnlMu = (float *)pvIrlMalloc(iKrnlSize*sizeof(float), 
							"ReadKrnl:KrnlMu");

		fread(psKrnl->pfKrnl,sizeof(float),iKrnlSize,fpIn);
		vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading Krnl");
		fread(psKrnl->pfKrnlMu,sizeof(float),iKrnlSize,fpIn);
		vCheckIoError(fpIn, pchInFname, "ReadKrnl", ": reading KrnlMu");
		if (bSwapped){
			vSwapWords((void *)psKrnl->pfKrnl, iKrnlSize);
			vSwapWords((void *)psKrnl->pfKrnlMu, iKrnlSize);
		}
		vPrintMsg(6,"    sum of krnl=%.4g, sum of krnlmu=%.4g\n",
				sum_float(psKrnl->pfKrnl,iKrnlSize),
				sum_float(psKrnl->pfKrnlMu,iKrnlSize));
			
	}else{
		psKrnl->pfKrnl = NULL;
		psKrnl->pfKrnlMu=NULL;
	}
	fclose(fpIn);
	return psKrnl;
}

#ifdef ALLOW_INDIRECT_SRF_KRNL
KrnlData_t *psReadKrnlIndirect(char *pchParFile, int iParLine)
{
	FILE *fp;
	int iline;
	char c;
	float fPriFac, fScatFac, fCntrPixFac, fKrnlMu0, fKrnlPixSize;
	float *pfKrnl, *pfKrnlMu;
	char pchKrnlFname[256], pchMuFname[256];
	int ixdim, iydim, izdim, iNumXYoffs, iNumDists;
	KrnlData_t *psKrnl;

	fp = fopen(pchParFile,"r");
	if (fp == NULL)
		vErrorHandler(ECLASS_FATAL,ETYPE_IO, 
				"SetSrfParms","error opening srf parameterfile");
	for(iline=1; iline<iParLine; ++iline)
	  while((c = getc(fp)) != '\n' && c != EOF);
	if (fscanf(fp,"%256s %256s %f %f %f %f %f",
			pchKrnlFname, pchMuFname, &fKrnlPixSize, 
			&fKrnlMu0, &fPriFac, &fScatFac, &fCntrPixFac) != 7)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "ReadKrnlIndirect",
			"Error parsing krnl parfile line %d in file: %s", 
			iParLine, pchParFile);
	vPrintMsg(6,
		"krnl=%s mu=%s pixsize=%.3f mu0=%.4f CntrPixFac=%.6g scfac=%.6g\n",
		pchKrnlFname, pchMuFname, fKrnlPixSize,fKrnlMu0,fCntrPixFac,fScatFac);
	fclose(fp);

	pfKrnl = readimage3d(pchKrnlFname, &ixdim, &iydim, &izdim);
	if (ixdim != iydim){
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "ReadKrnlIndirect",
			"x and y dimensions of kernel not equal");
	}
	iNumDists = izdim;
	iNumXYoffs = ixdim;
	pfKrnlMu = readimage3d(pchMuFname, &ixdim, &iydim, &izdim);
	if (ixdim != iNumXYoffs || iydim != iNumXYoffs || izdim != iNumDists){
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "ReadKrnlIndirect",
			"kernel and kernelmu must be same size\n");
	}
	psKrnl = pvIrlMalloc(sizeof(KrnlData_t),"ReadKrnlIndirect:psKrnl");
	psKrnl->pchTitle = pvIrlMalloc(KRNL_TITLESIZE+1, 
				"ReadKrnlIndirect:pchTitle");
#ifdef NOSNPRINTF
	sprintf(psKrnl->pchTitle, 
#else
	snprintf(psKrnl->pchTitle, KRNL_TITLESIZE, 
#endif
			"krnl from %s, line %d (krnl=%s, mu=%s)",
			pchParFile, iParLine,pchKrnlFname,pchMuFname);
	psKrnl->pchIsotope = pchIrlStrdup("unknown");

	psKrnl->iNumXYoffs = iNumXYoffs;
	psKrnl->iNumDists =  iNumDists;
	psKrnl->fPixSize = fKrnlPixSize;
	psKrnl->fHighEnergy = -1.0;
	psKrnl->fLowEnergy = -1.0;
	psKrnl->fPriFac = fPriFac;
	psKrnl->fScatFac = fScatFac;
	psKrnl->fCntrPixFac = fCntrPixFac;
	psKrnl->fKrnlMu0=fKrnlMu0;
	psKrnl->pfKrnl=pfKrnl;
	psKrnl->pfKrnlMu=pfKrnlMu;
	vPrintKrnlInfo(psKrnl);
	return(psKrnl);
}
#endif

KrnlData_t *psFreeKrnlData(KrnlData_t *psKrnl,int bFreeArrays)
/* Frees kernel title and isotope strings and Krnl structure. If
 * bFreeArrays is true then it also frees the pfKrnl and pfKrnlMu
 * memory.
 */
{
	IrlFree(psKrnl->pchTitle);
	IrlFree(psKrnl->pchIsotope);
	if (bFreeArrays){
		IrlFree(psKrnl->pfKrnl);
		IrlFree(psKrnl->pfKrnlMu);
	}
	IrlFree(psKrnl);
	return (NULL);
}


