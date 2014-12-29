/************************************************************************
*									*
*	File :	vt_time.c						*
*	Note :	Functions for a time-domain simulatio of the vocal	*
*		tract with the sub-glottal system and a nasal side	*
*		branch.							*
*		Voiced sounds are synthesized by time-varying the	*
*		glottal area, Ag (cm2).					*
*		Fricative and aspirated sounds are generated		*
*		automatically.  The value of nonlinear resistances is	*
*		assumed to be proportional to DC flow Udc. The resistans*
*		and Udc are calculated using second-order equation,	*
*		based on a simple low frequency model of the vocal ract.*
*		The amplitude of random noise sources, at the glottis	*
*		and at the super-glottal constriction, then are		*
*		modulated by the cubic law, Udc**3/A**2.5 (Stevevs).	*
*		(Note: Flanagan proposes that the amplitude of noise	*
*		sources is proportional to the square of Re, which is a	*
*		Reynolds number, i.e., Udc*Udc/Ag and Udc*Udc/Ac,	*
*		respectively.)						*
************************************************************************/

//#include "stdafx.h"

#include	<stdlib.h>
//#include <malloc.h>
#include	<math.h>

#include	"vtsimul.h"
#ifndef __max
#define __max(a,b)	((a<b)?a:b)
#endif

/************(Global variables which can be passed to main)**************/

extern	int	nAcc;		/* the main constriction section	*/

double	Rg_v,		/* viscous resistance at glottis	*/
Rg_k,		/* Bernoulli resistence at glottis	*/
Rc_v,		/* viscous resistance at constriction	*/
Rc_k,		/* Bernoulli resistance at constriction	*/
Udc,		/* a DC air flow in cm3/s		*/
Ug,		/* flow at glottis in cm3/s		*/
Uac;		/* flow at supraglottal constriction	*/

/*********************(structure array definitions)***********************/

typedef struct { double	Rs,	/* series viscous (flow) resistance	*/
Ls,	/* series inductance (acoustic mass)	*/
els,	/* voltage source associated with Ls	*/
Ns,	/* dipole noise pressure sources	*/
Ca,	/* parallel capacitance (compliance)	*/
ica,	/* current source associated with Ca	*/
Ud,	/* parallel flow source du to dA/dt	*/
Rw,	/* wall mechanical resistance		*/
Lw,	/* wall mass (inductance)		*/
elw,	/* voltage source associated with Lw	*/
Cw,	/* wall compiance			*/
ecw,	/* voltage source associated with Cw	*/
Gw;	/* total wall conductance, 1/(Rw+Lw+Cw)	*/
}  td_acoustic_elements;

typedef	struct { double	s,	/* forces (interlaced voltage-current	*/
/* sources)				*/
w,	/* matrix coefficients			*/
x,	/* variables (interlaced U and P's)	*/
S,	/* s after elimination procedure	*/
W;	/* w after elimination procedure	*/
}  td_linear_equation;

/******( global constants and variables for functions in this file )*****/

static	int	deci;	     /* decimation rate = simfrq/smpfrq */
static	double	dt_sim;
static	double	K_Bernoulli, K_viscous_glotte, K_viscous_const;
static	double	La, Ca, Grad, Srad, Rw, Lw, Cw, Kr;
static	double	Ag2, Ac2, rv_sum, rk_sum;    /* flow calculation*/
static	double	Rt, Lt, Ct;		     /* tracheal walls  */
static	double	Rlungs;		/* (No lung cavity simulated) */
static	int	Ng;
static	constriction	Ac_pharynx, Ac_mouth;
constriction	Ac;	/* detected super_glottal constriction	*/

/* vocal tract : {tracha} + glottis + pharynx + {mouth} */

/** fixed tracheal tube **/
static	int			ntr1;
static	td_acoustic_elements	*actr;

/** glotis (a single section) **/
static	int			nglt0 = 1, nglt1=2;
static	area_function		afglt[1], dafglt[1];
static	td_acoustic_elements	acglt[2];

/** vocal tract : pharynx + {mouth} **/
static	int			nvt0, nvt1;
static	area_function		*af, *daf;
static	td_acoustic_elements	*ac;

/** entire vocal tract : {trachea} +glottis + pharyx + {mouth} */
static	int			n0, n2, n3, n4, n5;
td_linear_equation	*eq;

/** mouth (bucal) tube **/
static	int			nbu0, nbu1, nbu3, nbu4, nbu5, nbu6;
static	area_function		*afbu, *dafbu;
static	td_acoustic_elements	*acbu;
td_linear_equation	*eqbu;

/** lip openng **/
static	double			Grad_lip, Lrad_lip, iLrad_lip;
static	double			U0_lip = 0, U1_lip;

/** total wall flow **/
double			U0_wall = 0, U1_wall;

/* nasal tract */
/** (fixed part) **/
static	int			nna0, nna1;
static	area_function		*afna;
static	td_acoustic_elements	*acna;

/** (nasal coupling section) **/
static	int			nnc1;
static	area_function		*afnc, *dafnc;
static	td_acoustic_elements	*acnc;

/** (entire nasal tract) **/
static	int			nnt3, nnt4, nnt5, nnt6;
static	td_linear_equation	*eqnt;
/** nostrils **/
static	double			Grad_nos, Lrad_nos, iLrad_nos;
static	double			U0_nos, U1_nos;

/* arrays for noise sources, scale factors and filter memories */
static	float	noiseScale;
static	float	inMemGlt[4][3], outMemGlt[4][3], lowMemGlt;
static	float	inMemAc[4][3], outMemAc[4][3], lowMemAc;

/****************( Global variables declared in "main" )*****************/

extern	float	simfrq;		/* simulation frequency in Hz		*/
extern	float	smpfrq;		/* sampling freq. (simfrq = n*smpfrq)	*/

extern	float	Psub;		/* subglottal air pressure in cmH2O	*/
extern	float	lung_vol;	/* lung volume in cm3			*/

extern	int	ntr;		/* # of tracheal sections		*/
extern	area_function	aftr[];	/* area function from lungs to glottis	*/

extern	float	Ag;		/* grottis area in cm**2		*/
extern	float	xg;		/* glottis thickness in cm (=0.3 cm)	*/
extern	float	lg;		/* fold length in cm			*/

extern	int	nvt;		/* # of vocal tract sections		*/
extern	area_function	*afvt;	/* AF from glottis to lips (or to nbp)	*/

extern	int	nbp;	/* (nasal branch point) = # of paryngeal sections*/
extern	int	nnc;		/* # of coupling sections		*/
extern	int	nnt;		/* # of fixed nasal tract sections	*/
extern	area_function	*afnt;	/* fixed AF from left to right(nostrils)*/
extern	float	sinus_vol;	/* volume of the nasal side cavity (cm3)*/
extern	float	sinus_apt;	/* aperture of the sinus canal (cm2)	*/
extern	float	sinus_len;	/* length of the sinus canal (cm)	*/
extern	int	sinus_pos;	/* sinus position in sections from nbp	*/

/*( simulation options )*/

extern	int	subGLTsystem;		/* = OFF or ON			*/
extern	int	nasal_tract;		/* = OFF or ON			*/
extern	int	sinus_cavity;		/* = OFF or ON			*/
extern	int	wall;			/* = YIELDING or RIGID		*/
extern	int	rad_boundary;		/* = RL_CIRCUIT or SHORT_CIRCUIT*/
extern	int	wall_radiation;		/* = Off or ON			*/
extern	int	vocal_tract;		/* = STATIONARY or TIME_VARYING	*/
extern	int	dynamic_term;		/* = OFF or ON			*/

extern	int	noise_source;		/* = OFF or ON			*/
extern	int	noiseSourceLoc;		/* = # of sections downstream	*/
/*   from Ac.N			*/
extern	float	noiseAmp;		/* = bound for noise source:
								for glottal noise, -16dB	*/

/*( an extra heat loss factor in the nasal tract )*/

extern	float	extra_loss_factor;	/* = normally 50 % */

/*( physical constants )*/

extern	float	ro;		/* = 1.14e-3; air density, gm/cm**3	*/
extern	float	c;		/* = 3.5e+4; sound velocity, cm/s	*/
extern	float	eta;		/* = 1.4; adiabatic constant		*/
extern	float	cp;		/* = 0.24; specific heat, cal/gm.degree	*/
extern	float	lamda;		/* = 5.5e-5; heat conduction, cal/cm.sec.degree */
extern	float	mu;		/* = 1.86e-4; viscosity coef, dyne.sec/cm**2 */
extern	float	wall_resi;	/* = 1600.; wall resistance, gm/s/cm2	*/
extern	float	wall_mass;	/* = 1.5; wall mass per area, gm/cm2	*/
extern	float	wall_comp;	/* = 3.0e+5; wall compliance		*/
extern	float	trachea_resi;	/*   for walls of the tracheal tube	*/
extern	float	trachea_mass;
extern	float	trachea_comp;
extern	float	H2O_bar;
extern	float	Kc;		/* a coeff. for Bernoulli resistance	*/

/***************************( Local functions )***************************/

/*****
*	Function: clear_sources
*	Note	: Clear source terms of acoustic reactance elements.
*****/
void	clear_sources (
					   int	ns1,			/* # of sections + 1 */
					   td_acoustic_elements  e[] )	/* acoustic elements */
{
	int	i;

	for(i=0; i<ns1; i++)
	{  e[i].els = 0;
	e[i].Ns  = 0;
	e[i].ica = 0;
	e[i].elw = 0;
	e[i].ecw = 0;
	e[i].Ud  = 0;
	}
}

/*****
*	Function: clear_PU
*	Note	: Clear pressures and volume vlocities inside VT.
*****/
void	clear_PU (
				  int	ns5,			/* =2*ns+3, ou ns = # of sections */
				  td_linear_equation  eq[])	/* matrix elements */
{
	int	i;

	for(i=0; i<ns5; i++) eq[i].x = 0;
}

/*****
*	Function: nonzero_t
*	Note	: limit x to a small nonzero positive value.
*****/
float	nonzero_t( float x )
{

	if(x >= AREA_MIN)return(x);
	else             return(AREA_MIN);
}

/*****
*	Function: copy_initial_af_t
*	Note :	Copy an initial input vocal tract area function, af,
*		having nss ( = nph + nbu) sections onto af and, if
*		nasal_tract == 0N, afbu with with check zeros.
****/

void	copy_initial_af_t(void)
{
	int	i, j;

	afglt[0].A = nonzero_t(Ag);			/* glottis */
	afglt[0].x = nonzero_t(xg);

	for(i=0; i<nvt0; i++)			/* pharynx + {mouth} */
	{  af[i].A = nonzero_t(afvt[i].A);
	af[i].x = nonzero_t(afvt[i].x);
	}

	if(nasal_tract == ON)
	{  for(j=0; j<nbu0; i++, j++)			/* mouth */
	{  afbu[j].A = nonzero_t(afvt[i].A);
	afbu[j].x = nonzero_t(afvt[i].x);
	}
	for(i=0; i<nnc; i++)
	{  afnc[i].A = nonzero_t(afnt[i].A);		/* nasal coupling */
	afnc[i].x = nonzero_t(afnt[i].x);
	}
	for(i=0; i<nna0; i++)			/* fixed sections */
	{  afna[i].A = nonzero_t(afnt[i+nnc].A);
	afna[i].x = nonzero_t(afnt[i+nnc].x);
	}
	}
}

/*****
*	Function: dax
*	Note :	Using the current vocal tract area function, afvt[], and the
*		previous one, af[] and afbu[], compute increment or
*		decriment for section area and length.  The calculated
*		values are used to refresh area function during the deci
*		simulation cycles by the linear interpolation.
*
*		It calculates Ud also.
*****/

void	dax (void)
{
	int	i, j;
	float	fdeci;

	fdeci = (float)deci;
	/* glottis */
	dafglt[0].A = (nonzero_t(Ag) - afglt[0].A)/fdeci;
	dafglt[0].x = 0.;				/* fixed length */
	if(dynamic_term == ON)
		acglt[1].Ud = smpfrq*dafglt[0].A*afglt[0].x;

	if(vocal_tract == TIME_VARYING)
	{
		/* pharynx + {mouth} */
		for(i=0; i<nvt0; i++)
		{  daf[i].A = (nonzero_t(afvt[i].A) - af[i].A)/fdeci;
		daf[i].x = (nonzero_t(afvt[i].x) - af[i].x)/fdeci;
		if(dynamic_term == ON)
			ac[i+1].Ud = smpfrq
			*(daf[i].A*af[i].x+(af[i].A+daf[i].A)*daf[i].x);
		}

		/* mouth and nasal coupling section */
		if(nasal_tract == ON)
		{  for(j=0; j<nbu0; i++, j++)
		{  dafbu[j].A = (nonzero_t(afvt[i].A) - afbu[j].A)/fdeci;
		dafbu[j].x = (nonzero_t(afvt[i].x) - afbu[j].x)/fdeci;
		if(dynamic_term == ON)
			acbu[j+1].Ud = smpfrq
			*(dafbu[j].A*afbu[j].x+(afbu[j].A+dafbu[j].A)*dafbu[j].x);
		}
		for(i=0; i<nnc; i++)
		{  dafnc[i].A = (nonzero_t(afnt[i].A) - afnc[i].A)/fdeci;
		dafnc[i].x = 0.;			/* fixed length */
		if(dynamic_term == ON)
			acnc[i+1].Ud = smpfrq*dafnc[i].A*afnc[i].x;
		}
		}
	}
}

/****
*	function : find_constriction
*	Note :	Find the constriction along the vocal tract that has the
*		minimum non-zero cross_sectional area and is farest from
*		the glottis.  The detected constriction dimensions are used
*		in DC resistance and noise source calculation.
*
*		The minimum senction in the glottis and in the nasal tract
*		is not included.
****/
void	find_constriction(void)
{
	int i;

	for(i=0, Ac_pharynx.A=100.; i<nvt0; i++)
		if(af[i].A > AREA_MIN && af[i].A <= Ac_pharynx.A)
		{  Ac_pharynx.A = af[i].A;
	Ac_pharynx.x = af[i].x;
	Ac_pharynx.N = i;
		}

		if(nasal_tract == OFF)		/* pharynx+mouth cavities */
		{  Ac.A = Ac_pharynx.A;
		Ac.x = Ac_pharynx.x;
		Ac.N = Ac_pharynx.N;
		Ac.Loc = PHARYNX;
		}

		if(nasal_tract == ON)
		{
			for(i=0, Ac_mouth.A=100.; i<nbu0; i++) /* minimum area in mouth */
				if(afbu[i].A > AREA_MIN && afbu[i].A <= Ac_mouth.A)
				{  Ac_mouth.A = afbu[i].A;
			Ac_mouth.x = afbu[i].x;
			Ac_mouth.N = i;
				}

				if(Ac_mouth.A <= Ac_pharynx.A)
				{  Ac.A = Ac_mouth.A;
				Ac.x = Ac_mouth.x;
				Ac.N = Ac_mouth.N;
				Ac.Loc = MOUTH;
				}
				else
				{  Ac.A = Ac_pharynx.A;
				Ac.x = Ac_pharynx.x;
				Ac.N = Ac_pharynx.N;
				Ac.Loc = PHARYNX;
				}
		}
}

/****
*	function : refresh_af
*	note :	Refresh area function by linear interporation.  It is
*		used with "dax".  It seraches the constriction location
*		and refresh Ac and Nc also.
****/
void	refresh_af(void)
{
	int	i;

	for(i=0; i<nvt0; i++)			/* pharynx + {mouth} */
	{  af[i].A = nonzero_t(af[i].A + daf[i].A);
	af[i].x = nonzero_t(af[i].x + daf[i].x);
	}

	if(nasal_tract == ON)
	{  for(i=0; i<nbu0; i++)		/* mouth */
	{ afbu[i].A = nonzero_t(afbu[i].A + dafbu[i].A);
	afbu[i].x = nonzero_t(afbu[i].x + dafbu[i].x);
	}
	for(i=0; i<nnc; i++)			/* nasal-coupling */
	{  afnc[i].A = nonzero_t(afnc[i].A + dafnc[i].A);
	afnc[i].x = nonzero_t(afnc[i].x + dafnc[i].x);
	}
	}
	find_constriction();		/* find new Ac */
}

/*****
*	Function: acous_elements_t
*	Note  :	compute acoustic elements of the tansmission line for
*		the time_domain simulation.
*		e[0].{Rs, Ls, ..} corresponds to the left_half of the
*		initial section, and e[ns].{} to the right-half
*		of the last section.
*****/
void	acous_elements_t (
						  int			tract_typ,	/* TRACHEA, VOCAL_TRACT etc.*/
						  int			ns,		/* # of sections	*/
						  area_function		AF[],
						  td_acoustic_elements	e[])
{
	double	r0, L0, r1, L1;
	float	xda, ax;
	int	i, j;

	/* acoustic elements */

	r0 = 0.;
	L0 = 0.;
	for(i=0, j=1; i<ns; i++, j++)
	{  xda = AF[i].x / AF[i].A;
	if(tract_typ == GLOTTIS)
	{  r1 = K_viscous_glotte/pow(AF[i].A, 3.);
	L1 = 0;		/** acous. mass causes an explosion!!?? **/
	}
	else
	{  r1 = K_viscous_const*xda/AF[i].A;
	L1 = La*xda;
	}
	e[i].Rs = r0 + r1; r0 = r1;		/* series elements	*/
	e[i].Ls = L0 + L1; L0 = L1;
	e[j].Ca = Ca*AF[i].A*AF[i].x;	/* parallel elements	*/
	}
	e[ns].Rs = r1;		       /* right arm of the last section */
	e[ns].Ls = L1;

	/* add an extra loss to nasal sections */
	if(tract_typ == NOSE)
	{  e[0].Rs += 0.05;
	for(i=1; i<ns; i++) e[i].Rs += 0.1;
	e[ns].Rs += 0.05;
	}

	/* walls */
	for(i=0, j=1; i<ns; i++, j++)
	{  if(tract_typ == TRACHEA)
	{  ax       = AF[i].x * sqrt(AF[i].A);
	e[j].Rw = Rt/ax;
	e[j].Lw = Lt/ax;
	e[j].Cw = Ct/ax;
	e[j].Gw = 1.0/(e[j].Rw + e[j].Lw + e[j].Cw);
	}
	else
	{  if(wall == YIELDING)
	{  ax       = AF[i].x * sqrt(AF[i].A);
	e[j].Rw = Rw/ax;
	e[j].Lw = Lw/ax;
	e[j].Cw = Cw/ax;
	e[j].Gw = 1.0/(e[j].Rw + e[j].Lw + e[j].Cw);
	}
	else e[j].Gw = 0.;
	}
	}
}

/*****
*	Function: NTacous_elements_t
*	Note  :	compute acoustic elements of the tansmission line for
*		the time_domain simulation of the fixed part of the nasal
*		tract that can have a side cavity (Helmholtz resonator)
*		representing the nasal sinus.
*****/

void	NTacous_elements_t(
						   int			ns,		/* # of sections	*/
						   area_function		AF[],
						   td_acoustic_elements	e[])
{
	double	r0, L0, r1, L1;
	float	xda, ax;
	int	i, j, Psinus;

	/* sinus position in the fixed part ofthe NT */

	Psinus = __max(sinus_pos - nnc, 0);

	/* acoustic elements */

	r0 = 0.;
	L0 = 0.;
	for(i=0, j=1; i<ns; i++, j++)
	{  if(sinus_cavity == ON && j == Psinus)
	{  r1 = 0;
	L1 = 0;
	e[i].Rs = r0;  r0 = r1;
	e[i].Ls = L0;  L0 = L1;
	e[j].Ca = 0.;
	}
	else
	{  xda = AF[i].x / AF[i].A;
	r1 = K_viscous_const*xda/AF[i].A;
	L1 = La*xda;
	e[i].Rs = r0 + r1; r0 = r1;	/* series elements	*/
	e[i].Ls = L0 + L1; L0 = L1;
	e[j].Ca = Ca*AF[i].A*AF[i].x;	/* parallel elements	*/
	}
	}
	e[ns].Rs = r1;		       /* right arm of the last section */
	e[ns].Ls = L1;

	/* add an extra loss to nasal sections */
	e[0].Rs += 0.05;
	for(i=1; i<ns; i++) e[i].Rs += 0.1;
	e[ns].Rs += 0.05;

	/* walls */
	for(i=0, j=1; i<ns; i++, j++)
	{
		if(wall == YIELDING)
		{  ax       = AF[i].x * sqrt(AF[i].A);
		e[j].Rw = Rw/ax;
		e[j].Lw = Lw/ax;
		e[j].Cw = Cw/ax;
		e[j].Gw = 1.0/(e[j].Rw + e[j].Lw + e[j].Cw);
		}
		else e[j].Gw = 0.;
	}
	if(sinus_cavity == ON)		/* elements for side cavity */
	{  e[Psinus].Cw = Ca*sinus_vol;
	e[Psinus].Lw = La*sinus_len/sinus_apt;
	e[Psinus].Rw = K_viscous_const*sinus_len/pow(sinus_apt, 1.5);
	e[Psinus].Gw = 1.0/(e[Psinus].Rw + e[Psinus].Lw + e[Psinus].Cw);
	}
}

/*****
*	function : refresh_acoustic_elements
*	note :	Refresh values of acoustic elements for given vocal tract
*		area function including the nasal coupling section.
****/
void	refresh_acoustic_elements(void)
{
	float	Alip;

	/* vocal tract : parynx + {mouth} */

	acous_elements_t(VOCAL_TRACT, nvt0, af, ac);
	if(nasal_tract == ON)
	{  acous_elements_t(VOCAL_TRACT, nbu0, afbu, acbu);	/* mouth */
	acous_elements_t(NOSE, nnc, afnc, acnc);		/* n-coup.*/
	}
	/* radiation load at lips */

	if( rad_boundary == RL_CIRCUIT )
	{  if(nasal_tract == OFF) Alip = af[nvt0-1].A;
	else			  Alip = afbu[nbu0-1].A;
	Grad_lip = Grad*Alip;		/* radiation conductance*/
	Lrad_lip = Srad*sqrt(Alip);		/* radiation suceptance	*/
	}
	else
	{  Grad_lip = 5.;			/* short circuit	*/
	Lrad_lip = 0.;
	}
	if(nasal_tract == OFF) eq[n4].w = Grad_lip + Lrad_lip;
	else		      eqbu[nbu4].w = Grad_lip + Lrad_lip;
}


/*****
*	Function: eq_elements_t
*	Note  :	matrix coefficients of the linear equation, for a tube
*		having ns sections.  Its important function is to connect
*		tubes to form a single non-branched tube.
*****/
void	eq_elements_t (
					   int			tract_typ,
					   int			*k1,		/* index of right arm */
					   int			ns,		/* # of sections */
					   td_acoustic_elements	e[],
					   td_linear_equation	q[])
{
	int	i, j;

	if(*k1 == 1) q[*k1].w  =  e[0].Rs + e[0].Ls;	/* right arm */
	else	     q[*k1].w += (e[0].Rs + e[0].Ls);

	for(i=1, j=*k1; i<=ns; i++)
	{  q[++j].w = e[i].Ca;				/* parallel element */
	if(wall==YIELDING || tract_typ==TRACHEA) q[j].w += e[i].Gw;
	q[++j].w = e[i].Rs + e[i].Ls;		/* series element */
	}
	*k1 = j;
}

/*****
*	Function: force_constants
*	Note	: Refresh current/voltage source of the reactive elements
*		  and specify force constants (eq[i].s) in the linear
*		  equation, s = wx.
*****/
void	force_constants (
						 int			tract_type,
						 int			*k1,
						 int			ns,		/* # of sections */
						 td_acoustic_elements	e[],
						 td_linear_equation	q[] )
{
	int	i, j;
	double	Uw;

	/* Refresh current and voltage sources */

	/* right arm */
	e[0].els  = 2.0*e[0].Ls*q[*k1].x - e[0].els;

	for(i=1, j=*k1; i<=ns; i++)
	{  e[i].ica = 2.0*e[i].Ca*q[++j].x - e[i].ica; /* acoustic C */
	if( wall == YIELDING || tract_type == TRACHEA )
	{  Uw = e[i].Gw * (q[j].x - e[i].ecw + e[i].elw);
	e[i].elw  = 2.0*e[i].Lw*Uw - e[i].elw;
	e[i].ecw += 2.0*e[i].Cw*Uw;
	if(tract_type == VOCAL_TRACT) U1_wall += Uw;
	}
	e[i].els = 2.0*e[i].Ls*q[++j].x - e[i].els; /* acoustic L */
	}

	/* Copy force terms */

	if(*k1 == 1) q[*k1].s  =  e[0].els + e[0].Ns;	/* right arm */
	else	     q[*k1].s += (e[0].els + e[0].Ns);

	for(i=1, j=*k1; i<=ns; i++)
	{  q[++j].s = e[i].ica + e[i].Ud;
	if(wall == YIELDING || tract_type == TRACHEA)
		q[j].s += e[i].Gw*(e[i].ecw - e[i].elw);
	q[++j].s = e[i].els + e[i].Ns;
	}
	*k1 = j;
	/* zero noise source (NOT WORKING)*/
	/**for(i=0; i<=ns; i++)e[i].Ns = 0;**/
}

/******
*	Function: forward_elimination_t
*	Note	: a forward elimination proceduere to solve
*		  a linear algebric equation, s = Wx.
******/
void	forward_elimination_t(
							  int			n_fin,		/* final index */
							  td_linear_equation	q[])

{
	int	i;

	q[1].S = q[0].S + q[1].s;
	q[1].W = q[1].w;

	for(i=2; i<=n_fin; i++)
	{  q[i].S = q[i-1].S + q[i-1].W*q[i].s;
	q[i].W = q[i-2].W + q[i-1].W*q[i].w;
	}
}

/******
*	Function: backward_substitution_t
*	Note	: a backward substitution proceduere to solve
*		  a linear algebric equation, s = wx.
******/
void	backward_substitution_t(
								int			n_ini,		/* initial index */
								td_linear_equation	q[] )
{
	int	i;

	for(i=n_ini; i>0; i--)
		q[i].x = (q[i].S - q[i-1].W*q[i+1].x)/q[i].W;
}

/******
*	Function: backward_elimination_t
*	Note	: a backward elimination proceduere to solve
*		  a linear algebric equation, s = Wx.
******/
void	backward_elimination_t(
							   int			n_ini,		/* initial index	 */
							   td_linear_equation	q[])

{
	int	i;

	q[n_ini].S = q[n_ini].s;
	q[n_ini].W = q[n_ini].w;

	for(i=n_ini-1; i>0; i--)
	{  q[i].S = -q[i+1].S + q[i+1].W*q[i].s;
	q[i].W =  q[i+2].W + q[i+1].W*q[i].w;
	}
}

/******
*	Function: forward_substitution_t
*	Note	: a backward or forward substitution proceduere to solve
*		  a linear algebric equation, s = wx.
******/
void	forward_substitution_t(
							   int			n_fin,		/* final index */
							   td_linear_equation	q[] )
{
	int	i;

	for(i=1; i<=n_fin; i++)
		q[i].x = (q[i].S + q[i+1].W*q[i-1].x)/q[i].W;
}

/***********(signal processing functions, decimater and filters)*********/

/*****
*	Functions : decimation
*	Note	: The following two functions, decim_init and decim, are
*		  to reduce sampling rate by a factor "deci"
*		  (= simfrq/smpfrq), using an FIR for the interpolation.
*		  The filter length is fixed to a odd samples ( e.g.,101)
*		  and its coefficients are calculated with "decim_init".
*		  The output sample is returned by value.
*****/

int	count_decim;
int	q_decim = 51, p_decim = 101;	/* q = (p-1)/2 + 1 */
float	h_decim[51],  v_decim[101];

int	decim_init( void )
{
	float	cutoff, hd;
	int	i, q1;
	float	temp, pi = 3.141593;

	count_decim = 0;
	q1 = q_decim - 1;
	temp = 2.0*pi/(p_decim-1);
	for( i=0; i<p_decim; i++) v_decim[i] = 0;

	cutoff = 0.9*pi/deci;			/* cutoff frequency */
	for( i=0; i<q1; i++)
	{  hd   = sin(cutoff*(i-q1))/(pi*(i-q1));
	h_decim[i] = hd*( 0.54 - 0.46*cos(temp*i) );
	}
	h_decim[q1] = 0.5*cutoff/pi;

	/* return constant delay in output samples */
	return( (int) ((float)q_decim/(float)deci +0.5) );
}

float	decim(
			  int   out_flag,	/* = 0 for storing x, = 1 for filtering */
			  float x   )	/* input sample with the rate of simfrq Hz */
{
	float	sum;
	int	i, j, k;

	/* Store input sample in the filter memory */

	if( count_decim == p_decim ) count_decim = 0;
	v_decim[count_decim] = x;

	/* Filtering and output y */

	if( out_flag == 1 )
	{  j = count_decim;
	k = count_decim - 1;
	sum = 0;
	for( i=0; i<q_decim; i++)
	{  --j; if(j == -1)      j = p_decim - 1;
	++k; if(k == p_decim) k = 0;
	sum = sum + h_decim[i]*(v_decim[j] + v_decim[k]);
	}
	}
	count_decim++;
	return( sum );
}


/****
*	function: filter2
*	note:	Second order filter (direct form I). The size of input &
*		output buffer is 3.
****/
void	filter2(
				float	*x,		/* input buffer, x[0] is input	*/
				float	*y,		/* output buffer, y[0] is output*/
				float	*a,		/* coefficints of denominator, a[0]=1	*/
				float	*b)		/* coefficints of numerator	*/
{
	y[2] = y[1]; y[1] = y[0];
	y[0]=b[0]*x[0]+b[1]*x[1]+b[2]*x[2] - a[1]*y[1]-a[2]*y[2];
	x[2] = x[1]; x[1] = x[0];
}


/****
*	function: frictionNoise3
*	note:	generation of filtered noise to mimic turbulent noise
*		inside the vocal tract. Computer generated random numbers
*		undergo two filters: a highpass to cutoff low freaquency
*		components below 200 Hz and a lowpass filter (-6dB/oct)
*		to shape noise at high frequencies.
*
*		The low pass filter is a simple first order integrater
*		with alpha = 0.8
*
*		The 3rd order highpass filter is used for the low cut.
*		The cutoff frequency is 200 Hz when the sampling frequency
*		is 10 kHz.
****/
float	frictionNoise3(
					   int	Entry,		/* if =0, initialization */
					   float	mem1[4][3],	/* fixed work array (input memories) */
					   float	mem2[4][3],     /* fixed work array (output memories)*/
					   float	*mem)		/* fixed memory for lowcut filter */
{
	static	float	a[2][3] ={{1.0000, -1.8977, 0.9069},
	{1.0000, -0.9067, 0.0000}};
	static	float	b[2][3] ={{0.9546, -1.9093, 0.9546},
	{0.9974, -0.9974, 0.0000}};
	static	float	max, alpha=0.8;
	float	s;
	int	i, j;

	if(Entry == 0)
	{ 
		//randomize();
		max = 2./(float)RAND_MAX;
		for(j=0; j<3; j++)
		{
			for(i=0; i<2; i++)		/* lowcut filter */
			{ if(i==0) mem1[0][0] = max*(float)rand() - 1.;
			else     mem1[i][0] = mem2[i-1][0];
			filter2(mem1[i], mem2[i], a[i], b[i]);
			}
			s = mem2[1][0];
			s += alpha*(*mem); *mem = s; 	/* -6dB/oct */
		}
	}
	for(i=0; i<2; i++)
	{ if(i==0) mem1[0][0] = max*(float)rand() - 1.;
	else     mem1[i][0] = mem2[i-1][0];
	filter2(mem1[i], mem2[i], a[i], b[i]);
	}
	s = mem2[1][0];
	s += alpha*(*mem); *mem = s;	/* -6dB/oct */

	return(s);
}


/****
*	function: frictionNoise7
*	note:	generation of filtered noise to mimic turbulent noise
*		iside the vocal tract. Computer generated random numbers
*		undergo two filters: a highpass to cutoff low freaquency
*		components below 500 Hz and a lowpass filter (-6dB/oct)
*		to shape noise at high frequencies.
*
*		The low pass filter is a simple first order integrater
*		with alpha = 0.8
*
*		The 7-th order highpass filter is used for the low cut.
*		The cutoff frequency is 0.02 (normalized by sampling
*		frequency.) Actaul cutoff frequency varyies depending on
*		the sampling frequency as below:
*			20 kHz: 400 Hz
*			16 kHz: 320 Hz
*			10 kHz: 200 Hz
****/
float	frictionNoise7(
					   int	Entry,		/* if =0, initialization */
					   float	inBuf[4][3],	/* fixed work array (input memories) */
					   float	outBuf[4][3],	/* fixed work array (output memories)*/
					   float	*mem)		/* fixed memory for lowcut filter */
{
	static	float	a[4][3] ={{1.000000, -0.836060, 0.000000},
	{1.000000, -1.697388, 0.724731},
	{1.000000, -1.772522, 0.801086},
	{1.000000, -1.893677, 0.924164}};

	static	float	b[4][3] ={{0.459015, -0.918030, 0.459015},
	{1.711121, -1.711121, 0.000000},
	{0.893433, -1.786865, 0.893433},
	{1.042725, -2.085449, 1.042725}};
	static	float	max, alpha=0.8;
	float	s;
	int	i, j;

	if(Entry == 0)
	{ //randomize();
		max = 2./(float)RAND_MAX;
		for(j=0; j<7; j++)
		{
			for(i=0; i<4; i++)		/* lowcut filter */
			{ if(i==0) inBuf[0][0] = max*(float)rand() - 1.;
			else     inBuf[i][0] = outBuf[i-1][0];
			filter2(inBuf[i], outBuf[i], a[i], b[i]);
			}
			s = outBuf[3][0];
			s += alpha*(*mem); *mem = s; 	/* -6dB/oct */
		}
	}
	for(i=0; i<4; i++)
	{ if(i==0) inBuf[0][0] = max*(float)rand() - 1.;
	else     inBuf[i][0] = outBuf[i-1][0];
	filter2(inBuf[i], outBuf[i], a[i], b[i]);
	}
	s = outBuf[3][0];
	s += alpha*(*mem); *mem = s;	/* -6dB/oct */

	return(s);
}


/******************( Functions called from a main )***********************/

/*****
*	Function : vtt_ini
*	Note :	Initialize the vocal-tract state. It returns the constant
*		delay due to the decimation filter.
*****/

int	vtt_ini (void)
{
	int	i, j, cnst_delay;
	float	Alip;
	float	pi = 3.141593;

	/*** Initialize decimation filter ***/

	deci = simfrq/smpfrq;
	dt_sim = 1./simfrq;
	cnst_delay = decim_init();

	/*** Initialize random number generation ***/

	if(noise_source == ON)
	{  frictionNoise3(0, inMemGlt, outMemGlt, &lowMemGlt);
	frictionNoise3(0, inMemAc, outMemAc, &lowMemAc);
	}

	/*** Coefficients for computing acoustic-aerodynamic elements ***/

	/* flow registance */

	K_Bernoulli = Kc*ro;		/* Bernoulli kinetic resistance */
	/* The Kc value depends on the cross-section shape: =1.42 for
	the glottis (rectangular) and more close to 1. for a supraglottal
	constriction.  For the simplicity sake, the single value is
	used for the two cases. R = K_Bernoulli*Udc/(A*A) for full length
	*/

	K_viscous_glotte = 12.*lg*lg*xg*mu/2.;
	/* The K formula assumes a rectangular duct, and is applicable
	only at the glottis.  R = K_viscous_glotte/(A*A*A) for a half
	length */

	K_viscous_const = (8.*pi*mu)/2.0;	/* viscous resistance  */
	/* The K formula assumes a circular duct, and is applicable
	only at the superglottal constriction (Stevens, 1971).
	R = K_viscous_const*dx/(A*A) for a half length.  Note that
	viscous resistance is not same as "wall friction loss" used
	in the frequency domain calculation, which is a function of
	wall surface and frequency */

	/* acoustic elements */
	La = (2.0/dt_sim)*(ro/2.0);	/* acoustic mass (La)		*/
	Ca = (2.0/dt_sim)/(ro*c*c);	/* acoustic stiffness (1/Ca)	*/

	/* walls */
	Rw = wall_resi/(2.0*sqrt(pi));			/* VT walls */
	Lw = (2.0/dt_sim)*wall_mass/(2.0*sqrt(pi));
	Cw = (dt_sim/2.0)*wall_comp/(2.0*sqrt(pi));

	if(subGLTsystem == ON)
	{  Rt = trachea_resi/(2.0*sqrt(pi));		/* VT walls */
	Lt = (2.0/dt_sim)*trachea_mass/(2.0*sqrt(pi));
	Ct = (dt_sim/2.0)*trachea_comp/(2.0*sqrt(pi));
	}

	/* radiation impedance; 1/G_rad and 1/S_rad in parallel */
	Grad = (9.0*pi*pi)/(128.0*ro*c);	  /* conductance (G_rad) */
	Srad = (dt_sim/2.0)*(3.0*pi*sqrt(pi))/(8.0*ro);/* suceptance  (S_rad) */

	/* radiated sound pressure at 1 m */
	Kr = ro*simfrq/(2.0*pi*100.0);

	/*** address and memory allocation ***/

	/* tracheal tube */

	Ng = 0;
	if(subGLTsystem == ON)
	{  Ng    = ntr;
	ntr1  = ntr + 1;
	actr  = (td_acoustic_elements *) calloc(ntr1, sizeof(td_acoustic_elements));
	}

	/* main tract : parynx + {mouth} */

	if(nasal_tract == ON) nvt0 = nbp;	/* pharynx */
	else		      nvt0 = nvt;	/* parynx + mouth */
	nvt1 = nvt0 + 1;

	af  = (area_function *) calloc(nvt0, sizeof(area_function));
	daf = (area_function *) calloc(nvt0, sizeof(area_function));
	ac  = (td_acoustic_elements *) calloc(nvt1, sizeof(td_acoustic_elements));

	/* entire tract for equation elements : {trachea}+glottis+pyarynx+{mouth} */

	if(nasal_tract == ON) n0 = nbp;		/* pharynx */
	else		      n0 = nvt;		/* parynx + mouth */
	n0++;					/* +1 for glottal section */
	if(subGLTsystem == ON) n0 += ntr;	/* + tracheal sections */
	n2 = 2*n0;
	n3 = n2 + 1;
	n4 = n2 + 2;
	n5 = n2 + 3;
	eq  = (td_linear_equation *) calloc(n5, sizeof(td_linear_equation));

	if(nasal_tract == ON)
	{
		/* mouth cavity */
		nbu0 = nvt - nbp;			/* # of bucal sections */
		nbu1 = nbu0 + 1;
		nbu3 = 2*nbu0 + 1;
		nbu4 = 2*nbu0 + 2;
		nbu5 = 2*nbu0 + 3;
		nbu6 = 2*nbu0 + 4;
		afbu  = (area_function *) calloc(nbu0, sizeof(area_function));
		dafbu = (area_function *) calloc(nbu0, sizeof(area_function));
		acbu = (td_acoustic_elements *) calloc(nbu1, sizeof(td_acoustic_elements));
		eqbu = (td_linear_equation *) calloc(nbu6, sizeof(td_linear_equation));

		/* nasal coupling sections */
		nnc1 = nnc + 1;
		afnc   = (area_function *) calloc(nnc, sizeof(area_function));
		dafnc  = (area_function *) calloc(nnc, sizeof(area_function));
		acnc = (td_acoustic_elements *) calloc(nnc1, sizeof(td_acoustic_elements));

		/* fixed nasal tract */
		nna0 = nnt - nnc;
		nna1 = nna0 + 1;
		afna = (area_function *) calloc(nna0, sizeof(area_function));
		acna = (td_acoustic_elements *) calloc(nna1, sizeof(td_acoustic_elements));

		/* entire nasal tract */
		nnt3 = 2*nnt + 1;
		nnt4 = 2*nnt + 2;
		nnt5 = 2*nnt + 3;
		nnt6 = 2*nnt + 4;
		eqnt = (td_linear_equation *) calloc(nnt6, sizeof(td_linear_equation));
	}

	/***  Clear memories ***/

	/* current/voltage sources associated with reactances */

	if(subGLTsystem == ON) clear_sources(ntr1, actr);
	clear_sources(nglt1, acglt);
	clear_sources(nvt1, ac);
	iLrad_lip = 0.;

	if(nasal_tract == ON)
	{  clear_sources(nbu1, acbu);
	clear_sources(nnc1, acnc);
	clear_sources(nna1, acna);
	iLrad_nos = 0.;
	}

	/* clear volume velocities and central pressures */

	clear_PU(n5, eq);
	U1_lip = 0;
	if(nasal_tract == ON)
	{  clear_PU(nbu6, eqbu);
	clear_PU(nnt6, eqnt);
	U1_nos = 0;
	}

	/* set constant values for S and W */

	eq[0].S = H2O_bar*Psub;		/* lung air-pressure */
	eq[0].W = 1.0;
	if(nasal_tract == ON)
	{  eqbu[nbu5].S = 0.;
	eqbu[nbu5].W = 1.0;
	eqnt[nnt5].S = 0;
	eqnt[nnt5].W = 1.0;
	}

	/*** Copy initial area function of time-varying tubes ***/

	copy_initial_af_t();

	/*** find the constriction forward and farest from the glottis ***/

	find_constriction();

	/**** Fixed acoustic elements ****/

	/* lungs */
	Rlungs = 1.5/sqrt(600);			/* ?? */

	/* trachea */
	if(subGLTsystem == ON) acous_elements_t(TRACHEA, ntr, aftr, actr);

	/* fixed part of the nasal tract */
	if(nasal_tract == ON)
	{  NTacous_elements_t(nna0, afna, acna);
	/* radiation load at the nostrils */
	if(rad_boundary == RL_CIRCUIT)
	{  Grad_nos = Grad*afnt[nnt-1].A;
	Lrad_nos = Srad*sqrt(afnt[nnt-1].A);
	}
	else
	{  Grad_nos = 5000.;				/* short circuit */
	Lrad_nos = 0.;
	}
	eqnt[nnt4].w = Grad_nos + Lrad_nos;
	}

	/* stationary case */

	if(vocal_tract == STATIONARY) refresh_acoustic_elements();

	return( cnst_delay );
}

/*****
*	Function: vtt_sim
*	Note	: time-domain simulation of the vocal tract. Returns
*		  a single speech sample as value with the rate of
*		  smpfrq (Hz).
*****/

float	vtt_sim( void )
{
	int	i, j, time, r_arm;
	double	p, q, sound, sound_decim;

	/*** compute da and dx with a new area function, and Ud=d(A*x)/dt ***/

	dax();

	/*** Simulate deci (=simfrq/smpfrq) cycles with interpolation of a and x ***/

	for(time=0; time<deci; time++)
	{
		/*** refresh area function by a linear interpolation ***/

		afglt[0].A = nonzero_t(afglt[0].A + dafglt[0].A);	/* glottis */
		if(vocal_tract == TIME_VARYING) refresh_af();

		/*** time-varying acoustic elements ***/

		acous_elements_t(GLOTTIS, nglt0, afglt, acglt);
		if(vocal_tract == TIME_VARYING) refresh_acoustic_elements();

		/*** Bernoulli (kinetic) resistance at glottis and at constriction ***/

		/* compute a DC flow considering flow resistances at the 2 constrictions */

		Rg_v    = acglt[0].Rs + acglt[1].Rs;
		Ag2     = 1./(afglt[0].A*afglt[0].A);
		Ac2     = 1./(Ac.A*Ac.A);
		Rc_v    = 2.*K_viscous_const*Ac.x*Ac2;
		rv_sum  = Rg_v + Rc_v;
		rk_sum  = K_Bernoulli*(Ag2 + Ac2);
		Udc = (-rv_sum + sqrt(rv_sum*rv_sum + 4.*rk_sum*H2O_bar*Psub))
			/(2.*rk_sum);

		Rg_k = K_Bernoulli*Udc*Ag2;
		Rc_k = K_Bernoulli*Udc*Ac2;

		acglt[1].Rs += Rg_k;		/* add Bernoulli res. to glottis */
		if(Ac.Loc == PHARYNX) ac[Ac.N+1].Rs = Rc_v + Rc_k;
		if(Ac.Loc == MOUTH)   acbu[Ac.N+1].Rs = Rc_v + Rc_k;

		/* Noise sources at the exit of glottis and at the constriction */
		if(noise_source == ON)
		{
			if(Ac.Loc == PHARYNX)
			{ j = Ac.N+noiseSourceLoc+1;   /* at n sections downstream */
			if(j > nvt0)j = Ac.N + 1;    /* at the exit */
			ac[j].Ns = noiseAmp
				*frictionNoise3(1, inMemAc, outMemAc, &lowMemAc)
				*(Udc*Udc*Udc/pow(Ac.A, 2.5));
			}
			if(Ac.Loc == MOUTH)
			{ j = Ac.N+noiseSourceLoc+1;   /* at n sections downstream */
			if(j > nbu0)j = Ac.N + 1;    /* at the exit */
			acbu[j].Ns = noiseAmp
				*frictionNoise3(1, inMemAc, outMemAc, &lowMemAc)
				*(Udc*Udc*Udc/pow(Ac.A, 2.5));
			}
			/* -10dB = 0.315, -16dB = 0.158, -20dB = 0.1 */
			acglt[1].Ns = 0.1*noiseAmp
				*frictionNoise3(1, inMemGlt, outMemGlt, &lowMemGlt)
				*(Udc*Udc*Udc/pow(afglt[0].A, 2.5));
		}

		/*** w : matrix elements ***/

		/* vocal tract : {trachea} + glottis + pharynx + {mouth} */
		r_arm = 1;
		if(subGLTsystem == ON)
		{  eq_elements_t(TRACHEA, &r_arm, ntr, actr, eq);
		eq[1].w += Rlungs;
		}
		eq_elements_t(GLOTTIS, &r_arm, nglt0, acglt, eq);
		if(subGLTsystem == OFF) eq[1].w += Rlungs;
		eq_elements_t(VOCAL_TRACT, &r_arm, nvt0, ac, eq);

		if(nasal_tract == ON)
		{
			/* bucal and nasal tract */
			r_arm = 1;
			eq_elements_t(VOCAL_TRACT, &r_arm, nbu0, acbu, eqbu);
			r_arm = 1;
			eq_elements_t(NOSE, &r_arm, nnc, acnc, eqnt);
			eq_elements_t(NOSE, &r_arm, nna0, acna, eqnt);
		}

		/*** Refresh force constants ***/

		U1_wall = 0;		/* initialization for total wall "flow" */

		/* vocal tract : {trachea} + glottis + pharynx + {mouth} */
		r_arm = 1;
		if(subGLTsystem == ON)
			force_constants(TRACHEA, &r_arm, ntr, actr, eq);
		force_constants(GLOTTIS, &r_arm, nglt0, acglt, eq);
		force_constants(VOCAL_TRACT, &r_arm, nvt0, ac, eq);

		if(nasal_tract == OFF && rad_boundary == RL_CIRCUIT )
		{  iLrad_lip = 2.0*Lrad_lip*eq[n4].x + iLrad_lip;
		eq[n4].s = -iLrad_lip;		/* rad. admitance */
		}

		if( nasal_tract == ON )
		{  r_arm = 1;
		force_constants(VOCAL_TRACT, &r_arm, nbu0, acbu, eqbu);
		if( rad_boundary == RL_CIRCUIT )
		{  iLrad_lip = 2.0*Lrad_lip*eqbu[nbu4].x + iLrad_lip;
		eqbu[nbu4].s = -iLrad_lip;		/* rad. admitance */
		}

		r_arm = 1;
		force_constants(NOSE, &r_arm, nnc, acnc, eqnt);
		force_constants(NOSE, &r_arm, nna0, acna, eqnt);
		if( rad_boundary == RL_CIRCUIT )
		{  iLrad_nos = 2.0*Lrad_nos*eqnt[nnt4].x + iLrad_nos;
		eqnt[nnt4].s = -iLrad_nos;		/* rad. admitance */
		}
		}

		/*** solve s = Wx ***/

		if(nasal_tract == OFF)
		{  forward_elimination_t(n4, eq);
		eq[n4].x = eq[n4].S/eq[n4].W;
		backward_substitution_t(n3, eq);
		}
		else
		{  forward_elimination_t(n3, eq);
		backward_elimination_t(nbu4, eqbu);
		backward_elimination_t(nnt4, eqnt);

		p = eq[n3].S/eq[n3].W-eqbu[1].S/eqbu[1].W-eqnt[1].S/eqnt[1].W;
		q = eq[n2].W/eq[n3].W+eqbu[2].W/eqbu[1].W+eqnt[2].W/eqnt[1].W;
		eq[n4].x = eqbu[0].x = eqnt[0].x = p/q;

		backward_substitution_t(n3, eq);
		forward_substitution_t(nbu4, eqbu);
		forward_substitution_t(nnt4, eqnt);
		}
		/* decimation of the radiated sound */

		U0_lip = U1_lip;
		if(nasal_tract == OFF) U1_lip = -eq[n3].x;
		else			  U1_lip = -eqbu[nbu3].x;
		sound  = U1_lip - U0_lip;

		if(wall_radiation == ON)
		{ sound += (U1_wall - U0_wall);     /* radiation from walls */
		U0_wall = U1_wall;
		}

		if(nasal_tract == ON)
		{  U0_nos = U1_nos;		     /* add radiation from nose */
		U1_nos = -eqnt[nnt3].x;
		sound += (U1_nos - U0_nos);
		}

		if( time == deci - 1 ) sound_decim = decim( 1, Kr*sound );
		else                 	        decim( 0, Kr*sound );
	}

	/* sample Ug and Uac at fixed points, 1st and nAcc-th section */
	j = 1;
	if(subGLTsystem==ON) j += ntr;
	j = 2*j + 1;
	Ug = eq[j].x;

	if(nasal_tract==OFF)
	{ j = nAcc + 1;
	if(subGLTsystem==ON) j += ntr;
	j = 2*j + 1;
	Uac = eq[j].x;
	}
	else
	{ if(nAcc <= nbp)
	{ j = nAcc + 1;
	if(subGLTsystem==ON) j+=ntr;
	j = 2*j + 1;
	Uac = eq[j].x;
	}
	else
	{ j = 2*nAcc + 1;
	Uac = eqbu[j].x;
	}
	}

	return(sound_decim);
}

/*****
*	Function : vtt_term
*	Note :	free memories
****/

void	vtt_term ( void )
{
	if(subGLTsystem == ON) free(actr);

	free(af);
	free(daf);
	free(ac);

	free(eq);

	if(nasal_tract == ON)
	{  free(afbu);
	free(dafbu);
	free(acbu);
	free(eqbu);

	free(acna);
	free(afnc);
	free(dafnc);
	free(acnc);

	free(eqnt);
	}
}