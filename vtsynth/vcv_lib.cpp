/************************************************************************
*	File: VCV_lib.c														*
*	Note: This file contains interface functions for synthesis of VCV's.*							*
************************************************************************/

//#include "stdafx.h"

#include	<stdio.h>
#include	<stdlib.h>
//#include <conio.h>
//#include <malloc.h>
#include	<math.h>
#include	<string.h>

#include	"iowav.h"
#include	"vcv_lib.h"

#define	OFF	0
#define	ON	1

// This is a hack to make it compile on linux
// especially since these functions are not used
// in the mex version
#define strcmpi(x,y) (0)
#define _getch()  1
#define strlwr(x) 1
 
char	*stopConsonants="dfghjkpqstvxz";
const	float	oqF = 0.36,	/* opening & closing quotients for Fant's */
cqF = 0.26;

const	float	wp  = 2.0,	/* time warping coefficients for Maeda's. */
oqM = 0.5,
cqM = 0.2;
static	int		yIndexLen = 4;

/******************(externally defined global variables)******************/

extern	float	smpfrq;
extern	float	Ag;
extern	int		nvt;
extern	area_function	*afvt;
extern	int		noiseSourceLoc;

extern	char	TPfileDirectory[80];

/****************************(functions)*********************************/

/*****
*	Function: readVTAFalloc
*	Note: Read area function data (for vocal trac and nasal tract )
*		from a text file. Memory for AF wil be allocated in this
*		function, which must be freed in the caling program.
*		When the file is found, it returns 1, else 0
*****/
int	readVTAFalloc(
				  char	*AFpath,			/* complete file specification */
				  char*	AFtitle,			/* a text */
				  int		*nss,				/* # of sections */
				  area_function	**af )		/* area function */
{
	FILE	*in;
	char	*FileType;
	area_function	*p;
	float	a, dx;

	if((in = fopen(AFpath, "rt"))==NULL) return(0);
	fscanf(in, "%s", AFtitle );
	fscanf(in, "%d", nss);
	*af = (area_function *) calloc(*nss, sizeof(area_function));

	FileType = strtok(AFpath, ".");
	FileType = strtok(NULL, ".");

	if(!strcmpi(FileType, "ARE"))
	{ fscanf(in, "%f", &dx);		/* read dx */
	for(p = *af; p<(*af)+(*nss); p++)
	{ fscanf(in, "%f", &a);
	p->A = a;
	p->x = dx;
	}
	}
	else if(!strcmpi(FileType, "AX"))
	{ for(p = *af; p<(*af)+*nss; p++)	/* read (A, x) */
	{ fscanf(in, "%1*s%f%1*s %f%1*s", &a, &dx);
	p->A = a;
	p->x = dx;
	}
	}
	fclose (in);
	return(1);
}


/*****
*	Function: readVTAF
*	Note: Read area function data (for vocal trac and nasal tract )
*		from a text file, which contains
*		1) a text less than 80 characters
*		2) number of VT sections
*		3) area function with one of the following 4 formats:
*		If typ = .ARE, dx is constant, and followed by A's from Glt to Lip.
*		If     = .XY,  dx's followed by y's; A is derived from y.
*		If     = .XA,  dx's folloed by A's.
*		If     = .AX,  (A, dx)'s, which are paired in parentheses.
*		4) If the file is for a stop consonant, reads noise source location
*		relative to the constriction section. Consonant or not is identified
*		by the first letter of the file name.
*****/
void	readVTAF(
				 char	*AFpath,			/* complete file specification	*/
				 char	*AFtitle,			/* title od AF file		*/
				 int		nvt,				/* # of sections		*/
				 area_function	*af)		/* area function		*/
{
	FILE	*in;
	char	*FileType;
	int		nss, i;
	float	a, dx, y0, y, y_ave;
	char	cCode[2];

	if((in = fopen(AFpath, "rt"))==NULL)
	{ printf("AF path, %s, not found\n", AFpath);
	_getch();
	exit(-1);
	}
	fscanf(in, "%s", AFtitle );	/* read a title text */
	fscanf(in, "%d", &nss);		/* read # of sections */
	if(nss != nvt)
	{ puts("Unexpected number of VT sections");
	_getch();
	exit(-1);
	}
	FileType = strtok(AFpath, ".");
	FileType = strtok(NULL, ".");

	if(!strcmpi(FileType, "ARE"))
	{ fscanf(in, "%f", &dx);		/* read dx */
	for(i=0; i<nvt; i++)			/* read area and copy dx */
	{ fscanf(in, "%f", &af[i].A);
	af[i].x = dx;
	}
	}

	else if(!strcmpi(FileType, "XA"))
	{ for(i=0; i<nvt; i++) fscanf(in, "%f", &af[i].x);	/* x's */
	for(i=0; i<nvt; i++) fscanf(in, "%f", &af[i].A);	/* A's */
	}

	else if(!strcmpi(FileType, "XY"))
	{ for(i=0; i<nvt; i++) fscanf(in, "%f", &af[i].x);	/* x's */
	fscanf(in, "%f", &y0);		/* the initial y (glottis) */
	for(i=0; i<nvt; i++)
	{ fscanf(in, "%f", &y);
	y_ave = (y0 + y)/2;
	y0 = y;
	af[i].A = 2.0*pow(y_ave, 1.5);	/* sagittal-to-area */
	}
	}
	else if(!strcmpi(FileType, "AX"))
	{ for(i=0; i<nvt; i++)			/* read (A, x) */
	{ fscanf(in, "%1*s%f%1*s %f%1*s", &af[i].A, &af[i].x);
	}
	}

	else
	{ puts("File type for VT AF not found: Hit any key to quit");
	_getch();
	exit(-1);
	}

	/* Noise source location */

	cCode[0] = *(strrchr(AFpath, '\\')+1);
	cCode[1] = '\0';
	strlwr(cCode);
	if(strpbrk(cCode, stopConsonants)) fscanf(in, "%d", &noiseSourceLoc);
	fclose(in);
}

/*****
*	Function: FantGlottalArea
*	Note: To calculate oscilating glottal area (Ag cm2) at sample
*		point n during a fundamental period (t0 samples).
*		The glottal pulse shape is specified by Fant's model.
*		The calculated area is returned by value to the calling
*		program.
*****/

float	FantGlottalArea(
						float	Ap,			/* peak glottal area (cm2) */
						int		*T0)		/* the fundamental period in samples	*/
						/* The non-zero value signals a new		*/
						/* glottal cycle, then T0 will be zeroed */
{
	static	DWORD	n;
	static	float	amp;
	static	int	t1, t2, t3;
	static	float	a, b, A;
	double My_PI = 3.1415;

	float	t;

	if(*T0 > 0)		/* set a new glottal cycle */
	{ A  = Ap;
	t1 = oqF * (*T0);
	t2 = cqF * (*T0);
	t3 = t1 + t2;
	a  = My_PI/t1;
	b  = 1./(1. - cos(a*t2));
	*T0 = 0;
	n   = 0;
	}

	if(n < t1) amp = 0.5*A*(1.0 - cos(a*n));		/* opening */
	if(n >= t1 && n < t3)							/* closing */
	{ t = n - t1;
	amp = A*(1. - b + b*cos(a*t));
	}
	if(n >= t3) amp = 0.0;							/* closed */
	n++;
	return(amp);
}

/*****
*	Function : MaedaGlottalArea  (NOT WORKING RIGHT)
*	Note :	To calculate oscilating glottal area (Ag cm2) at sample
*		point n during a fundamental period (t0 samples).
*		The glottal pulse shape is specified by Maeda's model.
*		The calculated area is returned by value to the calling
*		program.
*****/

float	MaedaGlottalArea(
						 float	Ap,		/* peak glottal area (cm2)		*/
						 int	*T0)		/* the fundamental period in samples	*/
						 /* The non-zero value signals a new	*/
						 /* glottal cycle, then T0 will be zeroed*/
{
	static	DWORD	n;
	static	float	amp;
	static	int	t1, t2, t3;
	static	float	a, b, A;

	float	t;

	if(*T0 > 0)		/* set a new glottal cycle */
	{ A = 0.5*Ap;
	t1 = oqM * (*T0);
	t2 = cqM * (*T0);
	t3 = t1 + t2;
	a  = 3.141593/t1;
	b  = (float) (t1 - t2)/pow(t2, wp);
	*T0 = 0;
	n   = 0;
	}

	if(n < t1) amp = A*(1.0 - cos(a*n));		/* opening */
	if(n >= t1 && n < t3)				/* closing */
	{ t = n - t1;
	amp = A*(1. + cos(a*(t+pow(t,wp))));
	}
	if(n >= t3) amp = 0.0;				/* closed */
	n++;
	return(amp);
}

/*****
*	function: interpolX
*	note: Interolates a temporal pattern of a single variable.
*		The pattern is specified by a list of triplets; time(in ms),
*		value to be interpolated, and type of interpolation
*		(LIN, COS, EXP, etc) in a text file, for example,
*
*			3							// # of items
*			0.		100.	***			// always t=0. and type=***
*			54.		250.	COS
*			100.	155.	LIN
*
*		The input, file path and flag, etc. is passed by a pointer
*		to structure "interpoMem".  This function returns interpolated
*		value at the "sampCount" sample point, with the rate
*		specified by "smpfrq" in Hz.
*****/

float	interpolX(
				  DWORD		sampCount,		/* 0, 1, ...... */
				  interpoMem*	m)
{
	char	typ[4];
	float	time, value;
	double	dt;

	if(sampCount == 0)								/* Initialization */
	{
		m->cpCount = 0;
		if(m->cpCount >= m->NControlPoints)			/* no more data */
		{
			m->endFlag = ON;
			return(0.);
		}
		//fscanf(m->in, "%f %f %s", &time, &(m->b), typ);
		time = m->time[m->cpCount];
		m->b = m->value[m->cpCount];
		m->cpCount++;
		m->nextInterPoint = (smpfrq/1000.)*time;
	}

	if(sampCount == m->nextInterPoint)
	{
		if(m->cpCount == m->NControlPoints)			/* no more data */
		{
			m->endFlag = ON;
			return(m->b);
		}
		//fscanf(m->in, "%f %f %s", &time, &value, typ);
		time = m->time[m->cpCount];
		value = m->value[m->cpCount];
		m->Ttype = (transiType)(m->type[m->cpCount]);

		m->cpCount++;
		m->interPoint = m->nextInterPoint;
		m->nextInterPoint = (smpfrq/1000.)*time;

		//if(!strcmpi(typ, "SET"))m->Ttype = RESET;
		//if(!strcmpi(typ, "LIN"))m->Ttype = LINEAR;
		//if(!strcmpi(typ, "COS"))m->Ttype = COSINE;

		/* coefs. for interpolation */
		dt = (double)m->nextInterPoint - (double)m->interPoint;
		switch(m->Ttype)
		{ case RESET :
		m->b = value;
		break; 
		case LINEAR :
			m->a = (value - m->b)/dt;
			m->b = value;
			break;
		case COSINE :
			m->a = (value - m->b)/2;
			m->b = value;
			m->w = 3.141593/dt;
			break;
		default : break;
		}
	}

	/* Interpolation */
	dt = (double)sampCount - (double)m->nextInterPoint;
	switch(m->Ttype)
	{ case RESET  : return(0.);
	case LINEAR : return(m->a*dt + m->b);
	case COSINE : return(m->a*(cos(m->w*dt)-1) + m->b);
	default  : break;
	}
	return(0.);
}


/*****
*	function: interpolAx
*	note:
*		Interolates a temporal variation of a section (A, x) of VT area
*		function.
*		The pattern is specified by a list of quadplets; time(in ms),
*		value to be interpolated, and type of interpolation
*		(LIN, COS, EXP, etc) in a text file, for example,
*
*			3									// # of items
*			0.		2.3		5.5		***			// always t=0.
*			54.		4.7		5.0		COS
*			100.	2.3		5.5		LIN
*
*		The input, file path and flag, etc. is passed by a pointer
*		to structure "interpoMemTube".  This function returns
*		interpolated value, {A, x}, at the "sampCount" sample point,
*		with the rate specified by "smpfrq" in Hz.
*****/

area_function	interpolAx(
						   DWORD			sampCount,		/* 0, 1, ...... */
						   interpoMemTube*	m)
{
	char	typ[4];
	float	time;
	double	dt;
	area_function	af, afNull = {0., 0.};

	if(sampCount == 0)								/* Initialization */
	{ if((m->in = fopen(m->path, "rt")) == NULL)
	{ printf("%s not found", m->path);
	_getch();
	exit(-1);
	}
	fscanf(m->in, "%d", &(m->NControlPoints));
	m->cpCount = 0;
	if(m->cpCount >= m->NControlPoints)			/* no more data */
	{ m->endFlag = ON;
	fclose(m->in);
	return(afNull);
	}
	fscanf(m->in, "%f %f %f %s", &time, &(m->bA), &(m->bx), typ);
	m->cpCount++;
	m->nextInterPoint = (smpfrq/1000.)*time;
	}

	if(sampCount == m->nextInterPoint)
	{ if(m->cpCount == m->NControlPoints)			/* no more data */
	{ m->endFlag = ON;
	fclose(m->in);
	af.A = m->bA; af.x = m->bx;
	return(af);
	}
	fscanf(m->in, "%f %f %f %s", &time, &(af.A), &(af.x), typ);
	m->cpCount++;
	m->interPoint = m->nextInterPoint;
	m->nextInterPoint = (smpfrq/1000.)*time;
	if(!strcmpi(typ, "SET"))m->Ttype = RESET;
	if(!strcmpi(typ, "LIN"))m->Ttype = LINEAR;
	if(!strcmpi(typ, "COS"))m->Ttype = COSINE;

	/* coefs. for interpolation */
	dt = (double)m->nextInterPoint - (double)m->interPoint;
	switch(m->Ttype)
	{ case RESET :
	m->bA = af.A;
	m->bx = af.x;
	break;
	case LINEAR :
		m->aA = (af.A - m->bA)/dt;
		m->bA = af.A;
		m->ax = (af.x - m->bx)/dt;
		m->bx = af.x;
		break;
	case COSINE :
		m->aA = (af.A - m->bA)/2.;
		m->bA = af.A;
		m->w  = 3.141593/dt;
		m->ax = (af.x - m->bx)/2.;
		m->bx = af.x;
		break;
	default : break;
	}
	}

	/* Interpolation */
	dt = (double)sampCount - (double)m->nextInterPoint;
	switch(m->Ttype)
	{ case RESET  :
	return(afNull);
	case LINEAR :
		af.A = m->aA*dt + m->bA;
		af.x = m->ax*dt + m->bx;
		return(af);
	case COSINE :
		af.A = m->aA*(cos(m->w*dt)-1) + m->bA;
		af.x = m->ax*(cos(m->w*dt)-1) + m->bx;
		return(af);
	default  : break;
	}
	return(afNull);
}

/****
function: coscos
note: returns the value of cos(x) with a peak shape modification
****/
double coscos(int x, int N)
{
	double pi= 3.141593;
	float alpha=0.2;
	double r, c;
	int S;

	S=alpha*N;

	r=cos((pi/N)*x);
	if(x<=S)
	{
		c=0.01*(1.-cos((2.*pi/S)*x))/2.;
		r=r+c;
	}
	if(x>N-S)
	{
		c=0.01*(1.-cos((2.*pi/S)*(x-N+S)))/2.;
		r=r-c;
	}
	return(r);
}


/*****
*	function: interpolVTAF
*	note:
*	Interpolates a temporal variation of the VT area funvtion.
*	The variation is defined by a list of AF's with
*	the specifications of non-uniform sample points (in ms)
*	and the transition type, {***, LIN, COS, EXP}.
*	These infomation is described in a text file, for example,
*
*			3							// # of AF's
*			2							// # of sections
*										// line space
*			0.							// time (always 0)
*			(2.0, 8.0)	(6.0, 8.0)		// (A, x) <note>
*			***							// transi. type
*										// line space
*			150.
*			(6.0, 8.0)	(2.0, 8.0)
*			COS
*
*	<note> The (A, x) description can be replaced by the
*	name of standard area-function file with the type *.AX,
*	*.ARE, or *.XA.  The number of sections of AF in every
*	file must be equal to that specified in the second line.
*
*	The input, file path and flag, etc. is passed by a pointer
*	to structure "interpoMemAF".  This function returns
*	interpolated value, {A, x}, in the global variable "afvt"
*	with the number of sections, "nvt", at the "sampCount"
*	sample point, with the rate specified by "smpfrq" in Hz.
*****/

void	interpolVTAF(
					 DWORD			sampCount,		/* 0, 1, ...... */
					 interpoMemVTAF*	m)
{
	char	typ[4], text[80], filename[14];
	float	time;
	area_function	af[NSSMAX];
	int		i;
	char	ch;
	static	char	AFfilepath[80];
	double	dt;

	if(sampCount == 0)							/* Initialization */
	{ 
		m->cpCount = 0;
		if(m->cpCount >= m->NControlPoints)		/* no more data */
		{
			m->endFlag = ON;
			return;
		}

		afvt = (area_function *)calloc(nvt, sizeof(area_function));

		time = m->time[m->cpCount];
		m->Ttype = (transiType)(m->type[m->cpCount]);

		//fscanf(m->in, "%f", &time);
		//ch = fgetc(m->in);					/* the next line */
		//ch = fgetc(m->in);					/* get the first character */
		//ungetc(ch, m->in);					/* restore the character */

		if (m->nloc[m->cpCount]>=0)
			noiseSourceLoc = m->nloc[m->cpCount];
		//	printf("nloc[%d]\n",noiseSourceLoc);

		for(i=0; i<nvt; i++)
		{
			m->bA[i] = m->af[m->cpCount][i];
			m->bx[i] = m->dx;
		}

		//fscanf(m->in, "%s", typ);
		m->cpCount++;
		m->nextInterPoint = (smpfrq/1000.)*time;
	}

	if(sampCount == m->nextInterPoint)
	{
		if(m->cpCount == m->NControlPoints)			/* no more data */
		{
			m->endFlag = ON;
			return;
		}

		//fscanf(m->in, "%f", &time);
		time = m->time[m->cpCount];
		m->Ttype = (transiType)(m->type[m->cpCount]);

		//readVTAF(AFfilepath, text, nvt, af);
		//fscanf(m->in, "%s", typ);

		if (m->nloc[m->cpCount]>=0)
			noiseSourceLoc = m->nloc[m->cpCount];
		//	printf("nloc[%d]\n",noiseSourceLoc);

		for(i=0; i<nvt; i++)
		{
			af[i].A = m->af[m->cpCount][i];
			af[i].x = m->dx;
		}

		m->cpCount++;
		m->interPoint = m->nextInterPoint;
		m->nextInterPoint = (smpfrq/1000.)*time;

		//if(!strcmpi(typ, "SET"))m->Ttype = RESET;
		//if(!strcmpi(typ, "LIN"))m->Ttype = LINEAR;
		//if(!strcmpi(typ, "COS"))m->Ttype = COSINE;

		/* coefs. for interpolation */
		dt = (double)m->nextInterPoint - (double)m->interPoint;
		switch(m->Ttype)
		{
		case RESET :
			for(i=0; i<nvt; i++)
			{
				m->bA[i] = af[i].A;
				m->bx[i] = af[i].x;
			}
			break;
		case LINEAR :
			for(i=0; i<nvt; i++)
			{
				m->aA[i] = (af[i].A - m->bA[i])/dt;
				m->bA[i] = af[i].A;
				m->ax[i] = (af[i].x - m->bx[i])/dt;
				m->bx[i] = af[i].x;
			}
			break;
		case COSINE :
			for(i=0; i<nvt; i++)
			{
				m->aA[i] = (af[i].A - m->bA[i])/2.;
				m->bA[i] = af[i].A;
				m->ax[i] = (af[i].x - m->bx[i])/2.;
				m->bx[i] = af[i].x;
			}
			m->w = 3.141593/dt;
			m->N = dt;
			break;
		default : break;
		}
	}

	/* Interpolation */
	dt = (double)sampCount - (double)m->nextInterPoint;
	switch(m->Ttype)
	{
	case RESET  :
		return;
	case LINEAR :
		for(i=0; i<nvt; i++)
		{
			afvt[i].A = m->aA[i]*dt + m->bA[i];
			afvt[i].x = m->ax[i]*dt + m->bx[i];
		}
		return;
	case COSINE :
		for(i=0; i<nvt; i++)
		{
			/*afvt[i].A = m->aA[i]*(cos(m->w*dt)-1) + m->bA[i];*/
			/*afvt[i].x = m->ax[i]*(cos(m->w*dt)-1) + m->bx[i];*/
			afvt[i].A = m->aA[i]*(coscos(dt, m->N)-1) + m->bA[i];
			afvt[i].x = m->ax[i]*(coscos(dt, m->N)-1) + m->bx[i];
		}
		return;
	default  : break;
	}
}

/*****
*	function: FourTubeAF
*	note:
*	Computes a standard VT area function (with a fixed number
*	of sections) from the 4 tube model.  The results are in
*	"afvt", which are global variables used in VT  simulation.
*****/
void	FourTubeAF(
				   area_function	T1,
				   area_function	T2,
				   area_function	T3,
				   area_function	T4)
{
	area_function	af[4];
	float	dx;
	float	x0, x, z0, z, v;
	int		i, j;

	/* copy area-function specified by the 4 tubes */
	af[0].A = T1.A;	af[0].x = T1.x;
	af[1].A = T2.A; af[1].x = T2.x;
	af[2].A = T3.A;	af[2].x = T3.x;
	af[3].A = T4.A; af[3].x = T4.x;

	/* compute dx */
	for(x=0., i=0; i<4; i++) x = x + af[i].x;
	dx = x/nvt;
	for(i=0; i<nvt; i++) afvt[i].x = dx;

	/* compute areas */
	for(x=0, i=0; i<nvt; i++)
	{ x0 = x;						/* x lower limit */
	x += dx;						/* x upper limit */
	for(z=0, j=0; j<4; j++)
	{ z += af[j].x;
	if(x0 <= z) break;
	}

	if(x <= z || j == 3)
		afvt[i].A = af[j].A;		/* also x is in section j */
	else
	{ v = (z - x0)*af[j].A;
	for(++j; j<4; j++)
	{ z0 = z;
	z += af[j].x;
	if(x <= z) break;
	else v += af[j].x*af[j].A;
	}
	v += (x-z0)*af[j].A;
	afvt[i].A = v/dx;
	}
	}
}
