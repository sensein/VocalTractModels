/*
*********************************************************
File: AS2F.c
Author: Satrajit S. Ghosh
Email: satra@cns.bu.edu
Department of Cognitive and Neural Systems
Boston University
Date: 9.27.1999

Original code from Dr. Shinji Maeda's vtcalcs program

Description:
This purpose of this function is to take saggital area
functions (af0) and convert them to a area functions of a 
tube (afvt) and calculate the transfer function (tf) of 
the tube and the Formants (Fmt), their bandwidths (Bw)
and their Amplitudes (Aw).

Modifications:

Bugs:
Saggital to area does not work at present. Currently, it 
takes area functions which have different lengths and 
converts them to area functions of the same length.
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
static double	*pAFcfg;
static double	*pAFcfgmisc;

static mxArray	*oAf = 0;
static mxArray	*oTf = 0;
static mxArray	*oFmt = 0;
static mxArray	*oBw = 0;
static mxArray	*oAmp = 0;

static double	*pAf;
static double	*pTf;
static double	*pFmt;
static double	*pBw;
static double	*pAmp;

static int size = 0;

/* =====================================================================*/
/*	Function to determine if arguments are valid. 			*/
/*	Also fill in some pointers we will need later.			*/

static  int CheckArguments(int nlhs, mxArray *plhs[], 
				int nrhs, const mxArray *prhs[])
{     
	if (nrhs<4 || nrhs > 4 || nlhs < 5 || nlhs > 5){
		printf("Incorrect calling syntax:\n [AreaF,TransF,Fmts,BW,Ampl] = ");
		printf("AS2F(tractcfg,phycons,safcfg,anc)\n");
		printf("nlhs is %d, nrhs is %d.\n", nlhs, nrhs);
		return 1;
	}

/*------------ Check input is not empty ------------------------------ */	
	if (mxGetSize(prhs[0]) != 4)
		return 1;
	if (mxGetSize(prhs[1]) != 5)
		return 1;
	size = mxGetSize(prhs[2])/2;
	if (mxGetSize(prhs[3]) != 2)
		return 1;

	pTCcfg = mxGetPr(prhs[0]);
	pPCcfg = mxGetPr(prhs[1]);
	pAFcfg = mxGetPr(prhs[2]);
	pAFcfgmisc = mxGetPr(prhs[3]);
	
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
	float length = 0.0f;
	static area_function *af0;

	if (CheckArguments(nlhs, plhs, nrhs, prhs) ){
		mexErrMsgTxt("AS2F argument checking failed.");
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

	nss = (int)(pAFcfgmisc[0]);
	//printf("no. of sections[%d]\n",nss);
	anc = (float)(pAFcfgmisc[1]);
	//printf("nasal area[%f]\n",anc);


	/* Initialization */
	af0 = (area_function *) calloc( size, sizeof(area_function) );
	ns0 = size;
	//printf("Number of sections provided[%d]",ns0);
	length = 0.0f;
	for(i=0;i<ns0;i++)
	{
		af0[i].A = (float)(pAFcfg[2*i]);
		af0[i].x = (float)(pAFcfg[2*i+1]);
		length += af0[i].x;
	}
	//printf("Tube length[%.1f]",length);

	x = length/(float) nss;
	nph = (int)(9.0/x);	/* nasal brantch point is 9 cm above glottis */
	nbu = nss - nph;

	afvt  = (area_function *) calloc( nss, sizeof(area_function) );
	
	oAf = mxCreateDoubleMatrix(2,nss, mxREAL);
	pAf = mxGetPr(oAf);

	/* Compute VT profile and area function, and plot them */
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

	calplot_tf_FBA(nfrmmax, frm, bw, amp, &nfrms, tfunc, &ntf);

	//printf("Finished calculations.\n");


	oTf = mxCreateDoubleMatrix(ntf, 1, mxREAL);
	pTf = mxGetPr(oTf);
	for(i=0;i<ntf;i++)
	{
		pTf[i] = tfunc[i];
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
	plhs[1] = oTf; 
	plhs[2] = oFmt; 
	plhs[3] = oBw; 
	plhs[4] = oAmp; 

	// do cleanup
	//free( rad_re );
	//free( rad_im );
	free( afnt );
	free( af0 );
	free( afvt );
}
