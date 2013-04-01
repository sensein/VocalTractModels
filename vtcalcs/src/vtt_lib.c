


/***************************************************************************
*                                                                          *
*	File :	vtt_lib.c                                                  *
*	Note :	Functions for a time-domain simulatio of the vocal tract   *
*		with a nasal side branch.				   *
*                                                                          *
***************************************************************************/

#include	<stdlib.h>
//#include	<alloc.h>
#include	<math.h>

#include	<plot_lib.h>
#include	<vtsimul.h>

/*********************(stractue array definitions)***********************/

typedef struct { float	Rs,	/* series (flow) resistance		*/
			Ls,	/* series inductance (acoustic mass)	*/
			els,	/* voltage source associated with Ls	*/
			Ns,	/* dipole noise pressure sources	*/
			Ca,	/* parallel capacitance (compliance)	*/
			ica,	/* current source associated with Ca	*/
			Ud,	/* parzllel flow source dur to dA/dt	*/
			Rw,	/* wall mechanical resistance		*/
			Lw,	/* wall mass (inductance)		*/
			elw,	/* voltage source associated with Lw	*/
			Cw,	/* wall compiance			*/
			ecw,	/* voltage source associated with Cw	*/
			Gw;	/* total wall conductance, 1/(Rw+Lw+Cw)	*/
		}  td_acoustic_elements;

typedef	struct { float	s,	/* forces (interlaced voltage-current	*/
				/* sources)				*/
			w,	/* matrix coefficients			*/
			x,	/* variables (interlaced U and P's)	*/
			S,	/* s after elimination procedure	*/
			W;	/* w after elimination procedure	*/
		}  td_linear_equation;

/*******( global constants and variables for the following functions)*****/

	static	int	deci;	/* decimation rate = simfrq/smpfrq */
	static	float	H2O_bar= 980.39f;
	static	float	dt_sim;
	static	float	Rk, Rv, La, Ca, Grad, Srad, Rw, Lw, Cw, Kr;

/* pharyngeal tube */
	static	int			nph2, nph3, nph4;
	static	area_function		*afph, *dph;
	static	td_acoustic_elements	*acph;
	static	td_linear_equation	*eqph;
/* bucal tube	*/
	static	int			nbu2, nbu3, nbu4;
	static	area_function		*afbu, *dbu;
	static	td_acoustic_elements	*acbu;
	static	td_linear_equation	*eqbu;
	static	float			Grad_lips, Lrad_lips, irad_lips;
	static	float			U0_lips, U1_lips;
/* nasal tract */
	static	area_function		*dna;		     /* zeros */
	static	int			nna2, nna3, nna4;
	static	area_function		afnc[1], dnc[1];     /* NT inlet */
	static	td_acoustic_elements	*acna;
	static	float			Rs_na, Ls_na;
	static	td_linear_equation	*eqna;
	static	float			Grad_nose, Lrad_nose, irad_nose;
	static	float			U0_nose, U1_nose;

/****************( Global variables declared in "main" )*****************/

extern	float	Psub;		/* subglottal air pressure in cmH2O	*/
extern	float	Ag;		/* grottis area in cm**2		*/
extern	float	xg;		/* glottis thickness in cm (=0.3 cm)	*/

extern	float	simfrq;		/* simulation frequency in Hz		*/
extern	float	smpfrq;		/* sampling freq. (simfrq = n*smpfrq)	*/

extern	int	nbu;		/* # of sections in the bucal tube	*/
extern	int	nph;		/* # of sections in the phryngeal tube	*/
extern	int	source_loc;	/* source location by section number	*/
extern	area_function *afvt;	/* AF from the glottis to the lips	*/
extern	float	anc;		/* nasal coupling area in cm2		*/
extern	int	nna;		/* # of sections in the nasal tract	*/
extern	area_function *afnt;	/* AF from nostrils to coupling point	*/

	/*( simulation options )*/

extern	int	nasal_tract;		/* = OFF or ON			*/
extern	int	wall;			/* = YIELDING or RIGID		*/
extern	int	rad_boundary;		/* = RL_CIRCUIT or SHORT_CIRCUIT*/
extern	int	vocal_tract;		/* = STATIONARY or TIME_VARYING	*/
extern	int	dynamic_term;		/* = OFF or OFF			*/

	/*( an extra heat loss factor in the nasal tract )*/

extern	float	extra_loss_factor;	/* = 50 % */

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


/***************************( Local functions )***************************/

/*****
*	Function: nonzero_t
*	Note	: limit x to a small nonzero positive value.
*****/
float	nonzero_t( float x )
{
	float	limit = 0.0001f;

	if( x >= limit ) return( x );
	else             return( limit );
}

/*****
*	Function: copy_initial_af_t
*	Note :	Copy an initial input vocal tract area function, af,
*		having nss ( = nph + nbu) sections onto afph1 and afbu1
*		with with check zeros. Note that the address for the bucal
*		tube will be reversed, afbu1[0] corresponds to the lip
*		end and afbu1[nbu-1] to the nasal branch point.
****/

void	copy_initial_af_t ( void )
{
	int	i, j;

	for(i=0; i<nph; i++)
	{  afph[i].A = nonzero_t(afvt[i].A);
	   afph[i].x = nonzero_t(afvt[i].x);
	}
	for(i=0, j=nph+nbu-1; i<nbu; i++, j--)
	{  afbu[i].A = nonzero_t(afvt[j].A);
	   afbu[i].x = nonzero_t(afvt[j].x);
	}
	afnc[0].A = nonzero_t( anc );
	afnc[0].x = nonzero_t( afnt[nna-1].x );
}

/*****
*	Function: dax
*	Note :	Using a new vocal tract area function, afvt[], and the
*		previous one, afph[], and afbu[], compute increment or
*		decriment for section area and length.  The calculated
*		values are used to reflesh area function during the deci
*		simulation cycles by the linear interpolation.
*****/

void	dax ( void )
{
	int	i, j;

	for(i=0; i<nph; i++)
	{  dph[i].A = (afvt[i].A - afph[i].A)/(float) deci;
	   dph[i].x = (afvt[i].x - afph[i].x)/(float) deci;
	}
	for(i=0, j=nph+nbu-1; i<nbu; i++, j--)
	{  dbu[i].A = (afvt[j].A - afbu[i].A)/(float) deci;
	   dbu[i].x = (afvt[j].x - afbu[i].x)/(float) deci;
	}
	dnc[0].A = (anc - afnc[0].A)/(float) deci;
	dnc[0].x = 0.0;				/* constant length */
}

/*****
*	Function: Ud
*	Note :	Compute the parallele current source due to a change in
*		the volume of each section. Since the linear interpolation
*		of the area function is used during "deci" simulation
*		cycles, Ud is constant during the same cycles.
*****/

void	Ud ( void )
{
	int	i, j;

	for(i=0; i<nph; i++)
	   acph[i].Ud = smpfrq*(afvt[i].A*afvt[i].x - afph[i].A*afph[i].x);

	for(i=0, j=nph+nbu-1; i<nbu; i++, j--)
	   acbu[i].Ud = smpfrq*(afvt[j].A*afvt[j].x - afbu[i].A*afbu[i].x);

	acna[nna].Ud = smpfrq*(anc*afnt[nna-1].x - afnc[0].A*afnc[0].x);
}

/*****
*	Function: acou_mtrx
*	Note  :	compute acoustic elements of the tansmission line and
*		matrix coefficients of the linear equation, for a tube
*		having ns sections.
*****/
void	acou_mtrx (
	int			ns,		/* # of sections */
	area_function		af[],
	area_function		daf[],		/* increment or decrimant */
	td_acoustic_elements	ac[],
	td_linear_equation	eq[],
	float			r0,		/* arm of previous section */
	float			L0   )
{
	float	r1, L1, xda, ax;
	int	i, j;

/* compute the current area function by a linear interpolation */

	for(i=0; i<ns; i++)
	{  af[i].A += daf[i].A;
	   af[i].x += daf[i].x;
	}

/* acoustic elements */

	for(i=0, j=1; i<ns; i++, j++)
	{  xda = af[i].x / af[i].A;
	   r1  = Rv*xda/af[i].A;
	   L1  = La*xda;
	   ac[i].Rs = r0 + r1; r0 = r1;		/* series elements */
	   ac[i].Ls = L0 + L1; L0 = L1;
	   ac[j].Ca = Ca*af[i].A*af[i].x;	/* parallel elements */
	}
	ac[ns].Rs = r1;			/* left arm of the last section	*/
	ac[ns].Ls = L1;

	if( wall == YIELDING )			/* yielding walls */
	   for(i=0, j=1; i<ns; i++, j++)
	   {  ax       = (float)(af[i].x * sqrt(af[i].A));
	      ac[j].Rw = Rw/ax;
	      ac[j].Lw = Lw/ax;
	      ac[j].Cw = Cw/ax;
	      ac[j].Gw = (float)(1.0/(ac[j].Rw + ac[j].Lw + ac[j].Cw));
	   }

/* matrix coefficients */

	eq[1].w   = ac[0].Rs + ac[0].Ls;	/* right arm */

	for(i=1, j=1; i<=ns; i++)
	{  eq[++j].w = ac[i].Ca;
	   eq[++j].w = ac[i].Rs + ac[i].Ls;
	}
	if( wall == YIELDING )
	   for(i=1; i<=ns; i++) eq[2*i].w += ac[i].Gw;
}

/*****
*	Function: clear_sources
*	Note	: Clear source terms of acoustic reactance elements.
*****/
void	clear_sources (
	int	ns,			/* # of sections */
	td_acoustic_elements  ac[] )	/* acoustic elements */
{
	int	i;

	for(i=0; i<=ns; i++)
	{  ac[i].els = 0;
	   ac[i].Ns  = 0;
	   ac[i].ica = 0;
	   ac[i].elw = 0;
	   ac[i].ecw = 0;
	   ac[i].Ud  = 0;
	}
}

/*****
*	Function: clear_pu
*	Note	: Clear pressures and volume vlocities inside VT.
*****/
void	clear_pu (
	int	ns4,			/* =2*ns+2, ou ns = # of sections */
	td_linear_equation  eq[] )		/* matrix elements */
{
	int	i;

	for(i=0; i<=ns4; i++) eq[i].x = 0;
}

/*****
*	Function: force_constants
*	Note	: Refresh current/voltage source of the reactive elements
*		  and specify force constants (eq[i].s) in the linear
*		  equation, s = wx.
*****/
void	force_constants (
	int			ns,		/* # of sections */
	td_acoustic_elements	ac[],
	td_linear_equation		eq[] )
{
	int	i, j;
	float	Uw;

/* Refresh current and voltage sources */

	ac[0].els = (float)(2.0*ac[0].Ls*eq[1].x - ac[0].els);	   /* right arm */

	for(i=1, j=1; i<=ns; i++)
	{  ac[i].ica = (float)(2.0*ac[i].Ca*eq[++j].x - ac[i].ica); /* acoustic C */
	   ac[i].els = (float)(2.0*ac[i].Ls*eq[++j].x - ac[i].els); /* acoustic L */
	}
	if( wall == YIELDING )				   /* wall imp.  */
	{  for(i=1; i<=ns; i++)
	   {  Uw = ac[i].Gw * (eq[2*i].x - ac[i].ecw + ac[i].elw);
	      ac[i].elw  = (float)(2.0*ac[i].Lw*Uw - ac[i].elw);
	      ac[i].ecw += (float)(2.0*ac[i].Cw*Uw);
	   }
	}

/* Copy force terms */

	for(i=1, j=1; i<=ns; i++)
	{  eq[++j].s = ac[i].ica + ac[i].Ud;
	   eq[++j].s = ac[i].els + ac[i].Ns;
	}
	if(wall == YIELDING)
	for(i=1; i<=ns; i++)
	   eq[2*i].s += ac[i].Gw*(ac[i].ecw - ac[i].elw);
}

/******
*	Function: elimination_t
*	Note	: a forward elimination proceduere to solve
*		  a linear algebric equation, s = wx.
******/
void	elimination_t(
	int	i0,		/* i0 = 0 for bucal and nasal tubes,	*/
				/*    = 1 for pharynx tube		*/
	int	ns3,		/* = 2*ns+1, ou ns = number of sections	*/
	td_linear_equation	eq[])

{
	int	i, i1; 

	i1 = i0 +1;
	eq[i0].W = eq[i0].w;
	eq[i0].S = eq[i0].s;

	eq[i1].W =      (float)(1.0 + eq[i0].W*eq[i1].w);
	eq[i1].S = eq[i0].S + eq[i0].W*eq[i1].s;

	for(i=i0+2; i<=ns3; i++)
	{  eq[i].W = eq[i-2].W + eq[i-1].W*eq[i].w;
	   eq[i].S = eq[i-1].S + eq[i-1].W*eq[i].s;
	}
}

/******
*	Function: substitution
*	Note	: a backward substitution proceduere to solve
*		  a linear algebric equation, s = wx.
******/
void	substitution_t(
	int	i0,		/* i0 = 0 for bucal and nasal tubes,	*/
				/*    = 1 for pharynx tube		*/
	int	ns3,		/* =2*ns+1, ou ns = number of sections	*/
	td_linear_equation	eq[] )

{
	int	i, i1;

	i1 = i0 +1;
	for(i=ns3; i>=i1; i--)
	   eq[i].x = (eq[i].S - eq[i-1].W*eq[i+1].x)/eq[i].W;

	eq[i0].x = (eq[i0].S - eq[i1].x)/eq[i0].W;
}

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
	float	temp, pi = 3.141593f;

	count_decim = 0;
	q1 = q_decim - 1;
	temp = (float)(2.0*pi/(p_decim-1));
	for( i=0; i<p_decim; i++) v_decim[i] = 0;

	cutoff = (float)(0.9*pi/deci);			/* cutoff frequency */
	for( i=0; i<q1; i++)
	{  
		hd   = (float)(sin(cutoff*(i-q1))/(pi*(i-q1)));
		h_decim[i] = (float)(hd*( 0.54 - 0.46*cos(temp*i)));
	}

	h_decim[q1] = (float)(0.5*cutoff/pi);

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

/******************( Functions called from a main )***********************/

/*****
*	Function : vtt_ini
*	Note :	Initialize the vocal-tract state. It returns the constant
*		delay due to the decimation filter.
*****/

int	vtt_ini ( void )
{
	int	i, cnst_delay;
	float	pi = 3.141593f;

	nph2 = 2*nph; nph3 = nph2+1; nph4 = nph2+2;
	nbu2 = 2*nbu; nbu3 = nbu2+1; nbu4 = nbu2+2;
	nna2 = 2*nna; nna3 = nna2+1; nna4 = nna2+2;

	deci = (int)(simfrq/smpfrq);
	dt_sim = (float)(1./simfrq);
	cnst_delay = decim_init();

/*** Coefficients for computing acoustic-aerodynamic elements ***/

/* flow registance */
	Rk = (float)(1.2*ro);		/* kinetic resistance */

	/* The Rk value depends on the cross-section shape: =1.38 for
	   the glottis (rectangular) and =1. for a supragrottal
	   constriction.  For the simplicity sake, the single value is
	   used for the two cases. */

	Rv = (float)((0.8*pi*mu)/2.0);	/* viscus resistance  */
	/* The vr value depends on the shapes.  The difference is
	   relativly small, and the single value will be used. */

/* acoustic elements */
	La = (float)((2.0/dt_sim)*(ro/2.0));	/* acoustic mass (La)		*/
	Ca = (float)((2.0/dt_sim)/(ro*c*c));	/* acoustic stiffness (1/Ca)	*/

/* walls */
	Rw = (float)(wall_resi/(2.0*sqrt(pi)));
	Lw = (float)((2.0/dt_sim)*wall_mass/(2.0*sqrt(pi)));
	Cw = (float)((dt_sim/2.0)*wall_comp/(2.0*sqrt(pi)));

/* radiation impedance; 1/G_rad and 1/S_rad in parallel */
	Grad = (float)((9.0*pi*pi)/(128.0*ro*c));	  /* conductance (G_rad) */
	Srad = (float)((dt_sim/2.0)*(3.0*pi*sqrt(pi))/(8.0*ro));/* suceptance  (S_rad) */

/* radiated sound pressure at 1 m */
	Kr = (float)(ro*simfrq/(2.0*pi*100.0));

/*** memory allocations ***/

	afph = (area_function *) calloc( nph, sizeof(area_function) );
	dph  = (area_function *) calloc( nph, sizeof(area_function) );
	acph = (td_acoustic_elements *) calloc( nph+1, sizeof(td_acoustic_elements) );
	eqph = (td_linear_equation *) calloc( 2*nph+3, sizeof(td_linear_equation) );

	afbu = (area_function *) calloc( nbu, sizeof(area_function) );
	dbu  = (area_function *) calloc( nbu, sizeof(area_function) );
	acbu = (td_acoustic_elements *) calloc( nbu+1, sizeof(td_acoustic_elements) );
	eqbu = (td_linear_equation *) calloc( 2*nbu+3, sizeof(td_linear_equation) );

	dna  = (area_function *) calloc( nna, sizeof(area_function) );
	acna = (td_acoustic_elements *) calloc( nna+1, sizeof(td_acoustic_elements) );
	eqna = (td_linear_equation *) calloc( 2*nna+3, sizeof(td_linear_equation) );

/***  Initalization of memory terms  ***/

/* current/voltage sources associated with reactances */
	clear_sources( nph, acph );
	clear_sources( nbu, acbu );
	irad_lips = 0;
	clear_sources( nna, acna );
	irad_nose = 0;

/* initial volume velocities and central pressures (the rest condition) */
	clear_pu( nph4, eqph );
	clear_pu( nbu4, eqbu );
	U1_lips = 0;
	clear_pu( nna4, eqna );
	U1_nose = 0;

/**** Acoustic and matrix elements ****/

	copy_initial_af_t();		/* copy the initial area function */
	dax();

/* pharyngeal tract */
	acou_mtrx( nph, afph, dph, acph, eqph, 0., 0.);
	Ag = nonzero_t( Ag );			/* add glottal resistance */
	eqph[1].w =  (float)(eqph[1].w +( Rv*xg/Ag + Rk*fabs(eqph[1].x) )/(Ag*Ag));

/* bucal cavity */
	acou_mtrx( nbu, afbu, dbu, acbu, eqbu, 0., 0.);

/* nasal tract */
	acou_mtrx( nna-1, afnt, dna, acna, eqna, 0., 0. );
	for(i=1; i<=nna3; i+=2) eqna[i].w += 0.1f;  /* add some extra loss */

	Rs_na = acna[nna-2].Rs;			  /* left arm of the inlet*/
	Ls_na = acna[nna-2].Ls;			  /* to the nasal tract.  */
	acou_mtrx(1, afnc, dnc, acna+nna-1, eqna+2*(nna-1), Rs_na, Ls_na);

/* Radiation loads */
	if( rad_boundary == RL_CIRCUIT )
	{  Grad_lips = Grad*afbu[0].A;		/* radiation conductance */
	   Lrad_lips = (float)(Srad*sqrt(afbu[0].A));	/* radiation suceptance  */
	   eqbu[0].w = Grad_lips + Lrad_lips;	/* rad. admitance        */

	   Grad_nose = Grad*afnt[0].A;
	   Lrad_nose = (float)(Srad*sqrt(afnt[0].A));
	   eqna[0].w = Grad_nose + Lrad_nose;
	}
	else
	{  eqbu[0].w = 5.0;			/* short circuit	 */
	   eqna[0].w = 5.0;
	}

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
	int	j;
	float	f, g, h, p, q, sound, sound_decim;

/*** compute da and dx with a new area function, and Ud=d(A*x)/dt ***/

	if( vocal_tract == TIME_VARYING)
	{  dax();
	   if( dynamic_term == ON ) Ud();
	}

/*** Simulate deci (=simfrq/smpfrq) cycles with interpolation of a and x ***/

	for(j=0; j<deci; j++)
	{

/*** solve s = Wx ***/

	   if( nasal_tract == ON )
	   {  elimination_t(1, nph3, eqph);
	      elimination_t(0, nbu3, eqbu);
	      elimination_t(0, nna3, eqna);

	      f = eqph[nph3].S/eqph[nph3].W;
	      g = eqbu[nbu3].S/eqbu[nbu3].W;
	      h = eqna[nna3].S/eqna[nna3].W;
	      p = f + g + h;
	      f = eqph[nph2].W/eqph[nph3].W;
	      g = eqbu[nbu2].W/eqbu[nbu3].W;
	      h = eqna[nna2].W/eqna[nna3].W;
	      q = f + g + h;
	      eqph[nph4].x = eqbu[nbu4].x = eqna[nna4].x = p/q;

	      substitution_t(1, nph3, eqph);
	      substitution_t(0, nbu3, eqbu);
	      substitution_t(0, nna3, eqna);
	   }
	   else
	   {  elimination_t(1, nph3, eqph);
	      elimination_t(0, nbu3, eqbu);

	      f = eqph[nph3].S/eqph[nph3].W;
	      g = eqbu[nbu3].S/eqbu[nbu3].W;
	      p = f + g;
	      f = eqph[nph2].W/eqph[nph3].W;
	      g = eqbu[nbu2].W/eqbu[nbu3].W;
	      q = f + g;
	      eqph[nph4].x = eqbu[nbu4].x = p/q;

	      substitution_t(1, nph3, eqph);
	      substitution_t(0, nbu3, eqbu);
	   }

/*** Refresh acoustic and matrix elements ***/

	   if( vocal_tract == TIME_VARYING )
	   {
/* pharyngeal tract */
	      acou_mtrx( nph, afph, dph, acph, eqph, 0., 0.);

/* bucal cavity */
	      acou_mtrx( nbu, afbu, dbu, acbu, eqbu, 0., 0.);
	      if( rad_boundary == RL_CIRCUIT )
	      { 
			  Grad_lips = Grad*afbu[0].A;
			  Lrad_lips = (float)(Srad*sqrt(afbu[0].A));
			  eqbu[0].w = Grad_lips + Lrad_lips;
	      }
	      else
		 eqbu[0].w = 10.0;		/* short circuit */
/* nasal inlet */
	      acou_mtrx(1, afnc, dnc, acna+nna-1, eqna+2*(nna-1), Rs_na, Ls_na);
	   }
/* add the glottal resistance (it alwalys time_varying) */
	   Ag = nonzero_t( Ag );
	   eqph[1].w = (float)(acph[0].Rs + acph[0].Ls
		     + (Rv*xg/Ag + Rk*fabs(eqph[1].x))/(Ag*Ag));

/*** Refresh force constants ***/

	   force_constants(nph, acph, eqph);
	   eqph[1].s = acph[0].els + H2O_bar*Psub;	/* right arm */

	   if( rad_boundary == RL_CIRCUIT )
	      irad_lips = (float)(2.0*Lrad_lips*eqbu[0].x + irad_lips);
	   force_constants(nbu, acbu, eqbu);
	   eqbu[0].s = -irad_lips;		/* rad. admitance */
	   eqbu[1].s = acbu[0].els;		/* right arm      */

	   U0_lips = U1_lips;
	   U1_lips = -eqbu[1].x;
	   sound   = U1_lips - U0_lips;

	   if( nasal_tract == ON )
	   {  
		   if( rad_boundary == RL_CIRCUIT )
			   irad_nose = (float)(2.0*Lrad_nose*eqna[0].x + irad_nose);
		   force_constants(nna, acna, eqna);
		   eqna[0].s = -irad_nose;		/* rad. admitance */
		   eqna[1].s = acna[0].els;		/* right arm      */
	
		   U0_nose = U1_nose;
		   U1_nose = -eqna[1].x;
		   sound   = sound + U1_nose - U0_nose;
	   }

/*** decimation of the radiated sound ***/

	   if( j == deci - 1 ) sound_decim = decim( 1, Kr*sound );
	   else                 	     decim( 0, Kr*sound );
	}

/*** return the radiated sound pressure ***/

	return( sound_decim );
}

/*****
*	Function : vtt_term
*	Note :	free memories
****/

void	vtt_term ( void )
{
	free( afph );
	free( dph );
	free( acph );
	free( eqph );

	free( afbu );
	free( dbu );
	free( acbu );
	free( eqbu );

	free( dna );
	free( acna );
	free( eqna );
}
