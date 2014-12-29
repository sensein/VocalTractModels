/*
*********************************************************
File: vtsynth.cpp
Author: Satrajit S. Ghosh
Email: satra@cns.bu.edu
Department of Cognitive and Neural Systems
Boston University
Date: 10.12.2000


Description:
The purpose of this program is to synthesize the utterance 
by interpolating through area functions.

Modifications:

Bugs:
*********************************************************
*/

#include	<stdio.h>
#include	<math.h>
#include	<stdlib.h>
//#include	<malloc.h>
#include	"matrix.h"

#include	"mex.h"

#include "vtconf.h"
#include "vcv_lib.h"

#define	NOISE_SOURCE_LEVEL 2.0e-8;
#define	mxGetSize(m)	(mxGetN(m) * mxGetM(m))

static double	*pnss;
static double	*pTAg0;
static double	*pTAgP;
static double	*pTF0;
static double	*pTAF1;
static double	*pTAF2;
static double	*pdx;

static mxArray	*oSIG = 0;

static double	*pSIG;

/* =====================================================================*/
/*	Function to determine if arguments are valid. 			*/
/*	Also fill in some pointers we will need later.			*/

static  int CheckArguments(int nlhs, mxArray *plhs[], 
						   int nrhs, const mxArray *prhs[])
{     
	if (nrhs<7 || nrhs > 8 || nlhs < 1 || nlhs > 1){
		printf("Incorrect calling syntax:\n [SIG] = ");
		printf("vtsynth(nss,TAg0,TAgP,TF0,TAF1,TAF2,dx)\n");
		printf("nlhs is %d, nrhs is %d.\n", nlhs, nrhs);
		return 1;
	}
	if (nrhs<8)
	    smpfrq = 10000.0f;
	else
	    smpfrq = *mxGetPr(prhs[7]);    
    simfrq = 3*smpfrq;
    
	/*------------ Check input is not empty ------------------------------ */	
	if (mxGetSize(prhs[0]) != 1)
		return 1;
	if (mxGetM(prhs[1]) != 3)
		return 1;
	if (mxGetM(prhs[2]) != 3)
		return 1;
	if (mxGetM(prhs[3]) != 3)
		return 1;
	if (mxGetM(prhs[4]) != 3)
		return 1;

	pnss = mxGetPr(prhs[0]);
	pTAg0 = mxGetPr(prhs[1]);
	pTAgP = mxGetPr(prhs[2]);
	pTF0 = mxGetPr(prhs[3]);
	pTAF1 = mxGetPr(prhs[4]);
	pTAF2 = mxGetPr(prhs[5]);
	pdx = mxGetPr(prhs[6]);

	/*-------------------- Create the output matrix ----------------------- */
	return 0;
}  

/* =======================================================================*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])  
{ 	
	DWORD	sampCount, xPoint;
	float	F0;
	int		T0;
	interpoMem	TAg0, TAgP, TF0;
	interpoMemVTAF	TAF;

	int		i,j;
	float	Ag0, AgP;
	float AgMem = 0;
	float amin;
	float fsig, *sig;
	int allocsz=0;

	float noiseSourceLevel = 0.5; /* 1.5; */
	
	if (CheckArguments(nlhs, plhs, nrhs, prhs) ){
		mexErrMsgTxt("vtsynth argument checking failed.");
		return ;
	}

	
	nbp = 9;
	/* simulation options */
	vocal_tract = TIME_VARYING;
	nasal_tract = OFF;
	wall = YIELDING;
	rad_boundary = RL_CIRCUIT;
	wall_radiation = OFF;
	noise_source = ON;		/* or OFF */
	noiseAmp = noiseSourceLevel*NOISE_SOURCE_LEVEL;


	nvt = (int)(*pnss);
	TAF.dx = (float)(*pdx);

	TAg0.NControlPoints	= mxGetN(prhs[1]);
	TAgP.NControlPoints	= mxGetN(prhs[2]);
	TF0.NControlPoints	= mxGetN(prhs[3]);
	TAF.NControlPoints	= mxGetN(prhs[4]);

	/* Create storage space for each structure */
	TAg0.time = (float *) calloc(TAg0.NControlPoints , sizeof(float));
	TAg0.value = (float *) calloc(TAg0.NControlPoints , sizeof(float));
	TAg0.type = (int *) calloc(TAg0.NControlPoints , sizeof(int));

	TAgP.time = (float *) calloc(TAgP.NControlPoints , sizeof(float));
	TAgP.value = (float *) calloc(TAgP.NControlPoints , sizeof(float));
	TAgP.type = (int *) calloc(TAgP.NControlPoints , sizeof(int));

	TF0.time = (float *) calloc(TF0.NControlPoints , sizeof(float));
	TF0.value = (float *) calloc(TF0.NControlPoints , sizeof(float));
	TF0.type = (int *) calloc(TF0.NControlPoints , sizeof(int));

	TAF.time = (float *) calloc(TAF.NControlPoints , sizeof(float));
	TAF.type = (int *) calloc(TAF.NControlPoints , sizeof(int));
	TAF.nloc = (int *) calloc(TAF.NControlPoints , sizeof(int));
	TAF.af = (float **) calloc(TAF.NControlPoints,sizeof(float *));
	for(i=0;i<TAF.NControlPoints;i++)
		TAF.af[i] = (float *) calloc(nvt,sizeof(float));

	sig = (float *)calloc(10000,sizeof(float));
	allocsz = 10000;

	/* Copy values over */
	for(i=0;i<TAg0.NControlPoints;i++)
	{
		TAg0.time[i]	= pTAg0[3*i];
		TAg0.value[i]	= pTAg0[3*i+1];
		TAg0.type[i]	= pTAg0[3*i+2];
	}
	for(i=0;i<TAgP.NControlPoints;i++)
	{
		TAgP.time[i]	= pTAgP[3*i];
		TAgP.value[i]	= pTAgP[3*i+1];
		TAgP.type[i]	= pTAgP[3*i+2];
	}
	for(i=0;i<TF0.NControlPoints;i++)
	{
		TF0.time[i]		= pTF0[3*i];
		TF0.value[i]	= pTF0[3*i+1];
		TF0.type[i]		= pTF0[3*i+2];
	}

	for(nAcc=0,i=0;i<TAF.NControlPoints;i++)
	{
		TAF.time[i] = pTAF1[3*i];
		TAF.type[i] = pTAF1[3*i+1];
		TAF.nloc[i] = pTAF1[3*i+2];

		/*if (TAF.nloc[i] >= 0)
		for(amin = 999., j=0;j<nvt;j++)
		{
		TAF.af[i][j] = pTAF2[nvt*i+j];
		if (amin >= TAF.af[i][j])
		{
		amin=TAF.af[i][j]; 
		nAcc=j;
		}
		}
		*/	
		for(j=0;j<nvt;j++)
			TAF.af[i][j] = pTAF2[nvt*i+j];

	}
	//printf("nAcc[%d]\n",nAcc);

	/* Initialize more */
	TAg0.endFlag = OFF;
	TAg0.cpCount = 0;
	TAgP.endFlag = OFF;
	TAgP.cpCount = 0;
	TF0.endFlag = OFF;
	TF0.cpCount = 0;
	TAF.endFlag = OFF;
	TAF.cpCount = 0;
	sampCount = 0;
	xPoint = 0;

	/* Calculate */
	do
	{
		/* compute glottal area by interpolations */
		Ag0 = interpolX(sampCount, &TAg0);
		AgP = interpolX(sampCount, &TAgP);
		F0  = interpolX(sampCount, &TF0);
		if(TF0.Ttype == RESET) 
			xPoint = TF0.nextInterPoint;
		if(xPoint == sampCount)
		{ 
			T0 = 0;
			if(F0 > 0.)T0 = (int)(smpfrq/F0 + .5);	/* refresh pitch */
			xPoint += T0;
		}
		Ag = FantGlottalArea(AgP, &T0) + Ag0;
		Ag = (Ag + 0.8*AgMem)/1.8; AgMem = Ag;		/* smoothing */

		/* compute VT area-function by interpolation */
		interpolVTAF(sampCount, &TAF);

		if(sampCount == 0)vtt_ini();	/* initialize simulater */
		fsig = vtt_sim();

		if (sampCount == allocsz)
		{
			sig = (float *)realloc(sig,(allocsz+10000)*sizeof(float));
			allocsz += 10000;
		}
		sig[sampCount] = fsig;

		sampCount++;

		if(TAg0.endFlag == ON && TAgP.endFlag == ON)
			break;
	}while(TAg0.endFlag == OFF || TAgP.endFlag == OFF);

	/* termination procedurs */
	vtt_term();
	free(afvt);

	/* Output */
	oSIG = mxCreateDoubleMatrix(sampCount, 1, mxREAL);
	pSIG = mxGetPr(oSIG);

	for(i=0;i<sampCount;i++)
	{
		pSIG[i] = sig[i];
	}

	/* Assign output pointers */
	plhs[0] = oSIG; 

	// do cleanup
	free(TAg0.time);
	free(TAg0.value);
	free(TAg0.type);

	free(TAgP.time);
	free(TAgP.value);
	free(TAgP.type);

	free(TF0.time);
	free(TF0.value);
	free(TF0.type);

	free(TAF.time);
	free(TAF.type);
	free(TAF.nloc);
	for(i=0;i<TAF.NControlPoints;i++)
		free(TAF.af[i]);
	free(TAF.af);
	free(sig);
}
