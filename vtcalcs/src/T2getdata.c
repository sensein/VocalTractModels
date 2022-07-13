/*
*********************************************************
File: T2getdata.c
Author: Satrajit S. Ghosh
Email: satra@cns.bu.edu
Department of Cognitive and Neural Systems
Boston University
Date: 9.27.1999

Original code from Dr. Shinji Maeda's vtcalcs program

Description:
This purpose of this function is to take two tube model
parameters UTpar and convert them to a area functions of a 
tube (afvt) and calculate the transfer function (tf) of 
the tube and the Formants (Fmt), their bandwidths (Bw)
and their Amplitudes (Aw).

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
//#include	"lam_lib.h"
#include	"vtconfig.h"
#include	"vsyn_lib.h"
#include	"mydir.h"

#include "vtcalcs.h"

#define	mxGetSize(m)	(mxGetN(m) * mxGetM(m))

static double	*pTCcfg;
static double	*pPCcfg;
static double	*pT2cfg;

static mxArray	*oAf = 0;
static mxArray	*oTfm = 0;
static mxArray	*oTff = 0;
static mxArray	*oFmt = 0;
static mxArray	*oBw = 0;
static mxArray	*oAmp = 0;

static double	*pAf;
static double	*pTfm;
static double	*pTff;
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
	if (nrhs<3 || nrhs > 3 || nlhs < 6 || nlhs > 6){
		printf("Incorrect calling syntax:\n [AreaF,TransF,Fmts,BW,Ampl] = ");
		printf("T2getdata(tractcfg,phycons,t2cfg)\n");
		printf("nlhs is %d, nrhs is %d.\n", nlhs, nrhs);
		return 1;
	}

/*------------ Check input is not empty ------------------------------ */	
	if (mxGetSize(prhs[0]) != 4)
		return 1;
	if (mxGetSize(prhs[1]) != 5)
		return 1;
	if (mxGetSize(prhs[2]) != 5)
		return 1;

	pTCcfg = mxGetPr(prhs[0]);
	pPCcfg = mxGetPr(prhs[1]);
	pT2cfg = mxGetPr(prhs[2]);
	
/*-------------------- Create the output matrix ----------------------- */
	return 0;
}  

/* =======================================================================*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])  
{ 	
	register int i;
	float x;
	int max_len = 5;
	int	max_vtlen = 40;
	int n1,n2;
	float dx1,dx2,ratio;

	if (CheckArguments(nlhs, plhs, nrhs, prhs) ){
		mexErrMsgTxt("T2getdata argument checking failed.");
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

	T2par.A1 = (float)(pT2cfg[0]);
	//printf("A1[%f]\n",T2par.A1);
	T2par.x1 = (float)(pT2cfg[1]);
	//printf("x1[%f]\n",T2par.x1);
	T2par.A2 = (float)(pT2cfg[2]);
	//printf("A2[%f]\n",T2par.A2);
	T2par.x2 = (float)(pT2cfg[3]);
	//printf("x2[%f]\n",T2par.x2);
	anc = (float)(pT2cfg[4]);
	//printf("nasal area[%f]\n",anc);

	/* acclocate memory for area function */
	afvt = (area_function *) calloc( max_vtlen, sizeof(area_function) );

	n1 = (int)(ceil( T2par.x1 ));
	n2 = (int)(ceil( T2par.x2 ));
	nss = n1 + n2;
	if( nss > max_vtlen )
	{
		ratio = (float) nss / (float) max_vtlen;
		n1 = (int)(n1*ratio);
		n2 = (int)(n2*ratio);
		nss = n1 + n2;
	}
	
	oAf = mxCreateDoubleMatrix(2,nss, mxREAL);
	pAf = mxGetPr(oAf);

	dx1 = T2par.x1/n1;
	dx2 = T2par.x2/n2;
	nph = 9;
	nbu = nss - nph;
	
	for(i=0; i<n1; i++)
	{
		afvt[i].x = dx1;
		afvt[i].A = T2par.A1;
		pAf[2*i] = afvt[i].A;
		pAf[2*i+1] = afvt[i].x;
	}

	for(i=n1; i<nss; i++)
	{
		afvt[i].x = dx2;
		afvt[i].A = T2par.A2;
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

	// do cleanup
	//free( rad_re );
	//free( rad_im );
	free( afnt );
	free( afvt );
}
