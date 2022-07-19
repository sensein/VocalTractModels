/*
*********************************************************
File: AMgetdata.c
Author: Satrajit S. Ghosh
Email: satra@cns.bu.edu
Department of Cognitive and Neural Systems
Boston University
Date: 9.27.1999

Original code from Dr. Shinji Maeda's vtcalcs program

Description:
This purpose of this function is to take LAM parameters 
(AMpar) and convert them to a area functions of a 
tube (afvt) and calculate the transfer function (tf) of 
the tube and the Formants (Fmt), their bandwidths (Bw)
and their Amplitudes (Aw) and the data required for 
plotting the AM model.

Modifications:
Transfer function's x-axis data, freq (tff) is added as an output.
Transfer functoin's y-axis data's name (tf) changed to tfm.

Bugs:

*********************************************************
*/

#include	<stdio.h>
#include	<math.h>
#include	"mex.h"
#include	"plot_lib.h"
#include	"lam_lib.h"
#include	"vtconfig.h"
#include	"vsyn_lib.h"
#include	"mydir.h"

#include "vtcalcs.h"

#define	mxGetSize(m)	(mxGetN(m) * mxGetM(m))

static double	*pTCcfg;
static double	*pPCcfg;
static double	*pAMcfg;

static mxArray	*oAf = 0;
static mxArray	*oTfm = 0;
static mxArray	*oTff = 0;
static mxArray	*oFmt = 0;
static mxArray	*oBw = 0;
static mxArray	*oAmp = 0;
static mxArray	*oPdat1 = 0;
static mxArray	*oPdat2 = 0;

static double	*pAf;
static double	*pTfm;
static double	*pTff;
static double	*pFmt;
static double	*pBw;
static double	*pAmp;
static double	*pPdat1;
static double	*pPdat2;

static int size = 0;

#define NP	 29					/* = np        */
extern float2D	ivt[NP];		/* VT inside contours             */
extern float2D	evt[NP];		/* VT exterior contours           */
extern float	lip_w, lip_h;			/* lip-tube width and height      */

/* =====================================================================*/
/*	Function to determine if arguments are valid. 			*/
/*	Also fill in some pointers we will need later.			*/

static  int CheckArguments(int nlhs, mxArray *plhs[], 
				int nrhs, const mxArray *prhs[])
{     
	if (nrhs<3 || nrhs > 3 || nlhs < 8 || nlhs > 8){
		printf("Incorrect calling syntax:\n [AreaF,TransF,Fmts,BW,Ampl,pdat1,pdat2] = ");
		printf("AMgetdata(tractcfg,phycons,p3cfg)\n");
		printf("nlhs is %d, nrhs is %d.\n", nlhs, nrhs);
		return 1;
	}

/*------------ Check input is not empty ------------------------------ */	
	if (mxGetSize(prhs[0]) != 4)
		return 1;
	if (mxGetSize(prhs[1]) != 5)
		return 1;
	if (mxGetSize(prhs[2]) != 8)
		return 1;

	pTCcfg = mxGetPr(prhs[0]);
	pPCcfg = mxGetPr(prhs[1]);
	pAMcfg = mxGetPr(prhs[2]);
	
/*-------------------- Create the output matrix ----------------------- */
	return 0;
}  

/* =======================================================================*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])  
{ 	
	register int i;
	float x;
	int max_len = 6;
	int p;
	int ns0 = 29;
	static area_function *af0;

	if (CheckArguments(nlhs, plhs, nrhs, prhs) ){
		mexErrMsgTxt("AMgetdata argument checking failed.");
		return ;
	}

	/* read nasal tract area function from the file */
	//printf("Reading nasal tract file.\n");
	read_af(NTAFpath , &nna, &afnt );
	//printf("Finished Reading nasal tract file.\n");

	/* Initialization */
	//read_rad();


	// update all constants and variables
	if((int)(pTCcfg[0]) == RL_CIRCUIT)	rad_boundary = RL_CIRCUIT;
	else	rad_boundary = SHORT_CIRCUIT;
	//printf("Rad_boundary[%d]\n",rad_boundary);
	if((int)(pTCcfg[1]) == YIELDING)	wall = YIELDING;
	else	wall = RIGID;
	//printf("wall[%d]\n",wall);
	if((int)(pTCcfg[2]) == ON)	nasal_tract = ON;
	else	nasal_tract = OFF;
	//printf("nasal_tract[%d]\n",nasal_tract);
	if((int)(pTCcfg[3]) == CLOSE)	glt_boundary = CLOSE;
	else	glt_boundary = OPEN;
	//printf("glt_boundary[%d]\n",glt_boundary);

	ro = (float)(pPCcfg[0]);
	//printf("Air density[%f]\n",ro);
	c = (float)(pPCcfg[1]);
	//printf("Sound velocity[%f]\n",c);
	wall_resi = (float)(pPCcfg[2]);
	//printf("wall_resi[%f]\n",wall_resi);
	wall_mass = (float)(pPCcfg[3]);
	//printf("wall_mass[%f]\n",wall_mass);
	wall_comp = (float)(pPCcfg[4]);
	//printf("wall_comp[%f]\n",wall_comp);

	for(i=0;i<7;i++)
	{
		AMpar[i] = (float)(pAMcfg[i]);
		//printf("AMpar[%d][%f]\n",i,AMpar[i]);
	}
	anc = (float)(pAMcfg[7]);
	//printf("nasal area[%f]\n",anc);


	/* Initialization */
	read_model_spec();
	af0 = (area_function *) calloc( ns0, sizeof(area_function) );
	nph = 9;
	nbu = 8;
	nss = nbu + nph;
	afvt  = (area_function *) calloc( nss, sizeof(area_function) );
	convert_scale();
	semi_polar();

	
	oAf = mxCreateDoubleMatrix(2,nss, mxREAL);
	pAf = mxGetPr(oAf);

	/* Compute VT profile and area function, and plot them */
	lam( AMpar );				/* profile	*/
	
	sagittal_to_area( &ns0, af0 );		/* area function */
	for(i=0;i<ns0;i++)
	{
		//printf("AreaNS0[%d][%f]\n",i,af0[i]);
	}
	appro_area_function( ns0, af0, nss, afvt);
	for(i=0;i<nss;i++)
	{
		//printf("AreaAF[%d][%f]\n",i,afvt[i]);
		pAf[2*i] = afvt[i].A;
		pAf[2*i+1] = afvt[i].x;
	}
	if( nasal_tract == ON )
	{  
		anc = (float) min( anc, afvt[nph].A );
		afvt[nph].A -= anc;
		pAf[2*nph] = afvt[nph].A;
		pAf[2*nph+1] = afvt[nph].x;
	}

	calplot_tf_FBA(nfrmmax, frm, bw, amp, &nfrms, tfmag, tffreq, &ntf);

	//printf("Finished calculations.\n");

	oPdat1 = mxCreateDoubleMatrix(4,NP, mxREAL);
	pPdat1 = mxGetPr(oPdat1);
	for(i=0;i<NP;i++)
	{
		pPdat1[4*i] = ivt[i].x;
		pPdat1[4*i+1] = ivt[i].y;
		pPdat1[4*i+2] = evt[i].x;
		pPdat1[4*i+3] = evt[i].y;
	}
	oPdat2 = mxCreateDoubleMatrix(1,2,mxREAL);
	pPdat2 = mxGetPr(oPdat2);
	pPdat2[0] = lip_w;
	pPdat2[1] = lip_h;


	oTfm = mxCreateDoubleMatrix(ntf, 1, mxREAL);
	pTfm = mxGetPr(oTfm);
	oTff = mxCreateDoubleMatrix(ntf, 1, mxREAL);
	pTff = mxGetPr(oTff);
	for(i=0;i<ntf;i++)
	{
		pTfm[i] = tfmag[i];
		pTff[i] = tffreq[i];
	}

	oFmt = mxCreateDoubleMatrix(nfrms, 1, mxREAL);
	pFmt = mxGetPr(oFmt);
	oBw = mxCreateDoubleMatrix(nfrms, 1, mxREAL);
	pBw = mxGetPr(oBw);
	oAmp = mxCreateDoubleMatrix(nfrms, 1, mxREAL);
	pAmp = mxGetPr(oAmp);
	for(i=0;i<nfrms;i++)
	{
		pFmt[i] = frm[i];
		pBw[i] = bw[i];
		pAmp[i] = amp[i];
	}

	/* Assign output pointers */
	plhs[0] = oAf; 
	plhs[1] = oTfm;
	plhs[2] = oTff;
	plhs[3] = oFmt; 
	plhs[4] = oBw;
	plhs[5] = oAmp; 
	plhs[6] = oPdat1; 
	plhs[7] = oPdat2;

	// do cleanup
	//free( rad_re );
	//free( rad_im );
	free( afnt );
	free( af0 );
	free( afvt );
}
