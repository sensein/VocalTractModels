// Save this file as the name of the function you want to implement with it
// eg myfun.c
// compile the file from Matlab as:
// mex myfun.c
/*
 *	Template file for mex.
 */

#include	<stdio.h>
#include	<math.h>
#include	"mex.h"

#define	mxGetSize(m)	(mxGetN(m) * mxGetM(m))


// input variable pointers
static double	*pi1;
static double	*pi2;
static double	*pi3;

// output variables
static mxArray	*o1 = 0;
static mxArray	*o2 = 0;
static mxArray	*o3 = 0;

// output variable pointers
static double	*po1;
static double	*po2;
static double	*po3;

static int size = 0;

/* =====================================================================*/
/*	Function to determine if arguments are valid. 			*/
/*	Also fill in some pointers we will need later.			*/

static  int CheckArguments(int nlhs, mxArray *plhs[], 
				int nrhs, const mxArray *prhs[])
{  
	if (nrhs<3 || nrhs > 3 || nlhs < 5 || nlhs > 5){
		printf("Incorrect calling syntax:\n [o1,o2,o3] = ");
		printf("mexfun(i1,i2,i3)\n");
		printf("nlhs is %d, nrhs is %d.\n", nlhs, nrhs);
		return 1;
	}

/*------------ Check input is not empty ------------------------------ */	
	size = mxGetSize(prhs[0]);

// Here you get pointers to your inputs
	pi1 = mxGetPr(prhs[0]);
	pi2 = mxGetPr(prhs[0]);
	pi3 = mxGetPr(prhs[0]);
	
/*-------------------- Create the output matrix ----------------------- */
	return 0;
}  

/* =======================================================================*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])  
{ 	
	register int i,j;
	int nrows = 2, ncols = 2;

	if (CheckArguments(nlhs, plhs, nrhs, prhs) ){
		mexErrMsgTxt("mexfun argument checking failed.");
		return ;
	}

	// These are real types you can complex types too
	o1 = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
	po1 = mxGetPr(o1);
	o2= mxCreateDoubleMatrix(nrows, ncols, mxREAL);
	po2= mxGetPr(o2);
	o3= mxCreateDoubleMatrix(nrows, ncols, mxREAL);
	po3= mxGetPr(o3;

	// Matlab uses column major
	for (j=0;j<ncols;j++)
	{
		for(i=0;i<nrows;i++)
		{
			po1[j*nrows+i] = rand();
		}
	}

	/* Assign output pointers */
	plhs[0] = o1; 
	plhs[1] = o2; 
	plhs[2] = o3; 

}
