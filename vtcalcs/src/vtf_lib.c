


/***************************************************************************
*                                                                          *
*	File :	vtf_lib.c                                                  *
*	Note :	Functions for a frequency-domain simulatio of the vocal    *
*		tract with a nasal side branch. The vocal tract transfer   *
*		ratio, (U_lip + U_nostrils)/U_glottis, is calculated for   *
*		a specified frequency.					   *
*									   *
*		An extra fixed resistance is added to each of the nasal	   *
*		sections.  It affects on the appearance of the vowel	   *
*		spectral modification due to the occurence of pole-zero	   *
*		pairs.							   *
*									   *
*		The vocal tract is excited at the glottis in the current   *
*		program.  But series pressure or/and parallel current	   *
*		sources can be placed anywhere inside the vocal tract by   *
*		specifying eq[2j].s for a parallele current source and	   *
*		by eq[2j+1].s for a series pressure source.		   *
*									   *
***************************************************************************/

#include	<always.h>

//#include	<grafix.h>
#include	<complex.h>
#include	<plot_lib.h>
#include	<vtsimul.h>

/********************(stracture array definitions)***********************/

typedef struct { float	Rs,	/* series resistance			*/
			Ls,	/* series inductance (acoustic mass)	*/
			Ca,	/* parallel capacitance (compliance)	*/
			G,	/* conduntance due to heat conduction   */
			Rw,	/* wall mechanical resistance		*/
			Lw,	/* wall mass (inductance)		*/
			Cw;	/* wall compliance			*/
		}  fd_acoustic_elements;

typedef	struct { cmplx	s,	/* interlaced voltage & current sources	*/
			w,	/* interlaced admittances & imedances	*/
			x;	/* interlaced P's and U's to be solved	*/
		}  fd_linear_equation;


/*************************( global variables )*****************************/

	static	float	Rs, La, Ca, G, Rw, Lw, Cw, Rrad, Lrad;
/* pharyngeal tube */
	static	int			nph2, nph3, nph4;
	static	area_function		*afph;
	static	fd_acoustic_elements	*acph;
	static	fd_linear_equation	*eqph;
/* bucal tube	*/
	static	int			nbu2, nbu3, nbu4;
	static	area_function		*afbu;
	static	fd_acoustic_elements	*acbu;
	static	fd_linear_equation	*eqbu;
/* nasal tract */
	static	int			nna2, nna3, nna4;
	static	fd_acoustic_elements	*acna;
	static	fd_linear_equation	*eqna;
/* normalized radiation impedance */
	float	*rad_re, *rad_im;
	static	float	dka;
	static	int	n_ka;
	static	float	lip_Rrad, lip_Lrad;
	static	float	nose_Rrad, nose_Lrad;
	static	float	lip_rad_ka, lip_rad_norm;
	static	float	nose_rad_ka, nose_rad_norm;


/******************( Global variables defined in "main" )****************/

extern	int	nbu;		/* # of sections in the bucal tube	*/
extern	int	nph;		/* # of sections in the phryngeal tube	*/
extern	area_function *afvt;	/* AF from the glottis to the lips	*/
extern	float	anc;		/* nasal coupling area in cm2		*/
extern	int	nna;		/* # of sections in the nasal tract	*/
extern	area_function *afnt;	/* AF from nostrils to coupling point	*/

	/*( simulation options )*/

extern	int	nasal_tract;
extern	int	wall;
extern	int	rad_boundary;
extern	int	glt_boundary;
extern	int	source_loc;

	/*( an extra loss factor to increas heat loss in the nasal tract )*/

extern	float	extra_loss_factor;		/* in % */

	/*( physical constants )*/

extern	float	ro;		/* = 1.14e-3, air density, gm/cm**3	*/
extern	float	c;		/* = 3.5e+4, sound velocity, cm/s		*/
extern	float	eta;		/* = 1.4, adiabatic constant		*/
extern	float	cp;		/* = 0.24, specific heat, cal/gm.degree	*/
extern	float	lamda;		/* = 5.5e-5, heat conduction, cal/cm.sec.degree */
extern	float	mu;		/* = 1.86e-4, viscosity coef, dyne.sec/cm**2 */
extern	float	wall_resi;	/* = 1600., wall resistance, gm/s/cm2	*/
extern	float	wall_mass;	/* = 1.5, wall mass per area, gm/cm2	*/
extern	float	wall_comp;	/* = 3.0e+5, wall compliance		*/

extern	char	*RADIMPpath;

/***************************( Local functions )***************************/

/*****
*	Function: nonzero
*	Note	: limit x to a small nonzero positive value.
*****/
float	nonzero( float x )
{
	float	limit = 0.0001f;

	if( x >= limit ) return( x );
	else             return( limit );
}

/*****
*	Function: copy_initial_af
*	Note :	Copy an initial input vocal tract area function, af,
*		having nss ( = nph + nbu) sections onto afph and afbu
*		with with check zeros. Note that the address for the bucal
*		tube will be reversed, afbu1[0] corresponds to the lip
*		end and afbu1[nbu-1] to the nasal brantch point.
****/

void	copy_initial_af ( void )
{
	int	i, j;

	for(i=0; i<nph; i++)
	{  afph[i].A = nonzero(afvt[i].A);
	   afph[i].x = nonzero(afvt[i].x);
	}
	for(i=0, j=nph+nbu-1; i<nbu; i++, j--)
	{  afbu[i].A = nonzero(afvt[j].A);
	   afbu[i].x = nonzero(afvt[j].x);
	}
}

/*****
*	Function: acou_elements
*	Note  :	compute coefficient values so that acoustic elements
*		of the transmission line can be calculated as a function
*		of a specified frequency value.
*****/
void	acou_elements (
	int			ns,		/* # of sections */
	area_function		af[],
	fd_acoustic_elements	ac[] )
{
	float	r0, r1, L0, L1, sqa, sqad, xda;
	int	i, j;

/* acoustic elements */

	r0 = 0;
	L0 = 0;

	for(i=0, j=1; i<ns; i++, j++)
	{  sqa  = (float)sqrt(af[i].A);
	   sqad = sqa*af[i].x;
	   xda  = af[i].x / af[i].A;
	   r1   = Rs*xda/sqa;
	   L1   = La*xda;
	   ac[i].Rs = r0 + r1; r0 = r1;		/* series elements */
	   ac[i].Ls = L0 + L1; L0 = L1;
	   ac[j].Ca = Ca*af[i].A*af[i].x;	/* parallel elements */
	   ac[j].G  = G*sqad;
	}
	ac[ns].Rs = r1;			/* left arm of the last section	*/
	ac[ns].Ls = L1;

	if( wall == YIELDING )			/* yielding walls */
	   for(i=0, j=1; i<ns; i++, j++)
	   {  sqad     = (float)(af[i].x*sqrt(af[i].A));
	      ac[j].Rw = Rw/sqad;
	      ac[j].Lw = Lw/sqad;
	      ac[j].Cw = Cw*sqad;
	   }
}

/*****
*	Function: eq_elements
*	Note  :	compute the values of the force constants and coefficients
*		for the linear equation, at a given frequency.
*****/
void	eq_elements (
	float			freq,		/* frequency in Hz	*/
	int			ns,		/* number of sections	*/
	fd_acoustic_elements	ac[],
	fd_linear_equation	eq[]  )
{
	float	sq_freq;
	float	re, im, r, wl, wc, wlc, d;
	int	i, j;

	sq_freq = (float)sqrt(freq);

	re  = sq_freq*ac[0].Rs;			/* right arm */
	im  =    freq*ac[0].Ls;
	eq[1].w = complex(re, im);

	for(i=1, j=1; i<=ns; i++)
	{  if( wall == YIELDING )
	   {  r   = ac[i].Rw;			/* parallel element */
	      wl  = freq*ac[i].Lw;
	      wc  = freq*ac[i].Cw;
	      wlc = (float)(wl - 1.0/wc);
	      d   = (float)(1.0/(r*r + wlc*wlc));
	      re  = sq_freq*ac[i].G + r*d;
	      im  = freq*ac[i].Ca - wlc*d;
	   }
	   else
	   {  re  = sq_freq*ac[i].G;
	      im  = freq*ac[i].Ca;
	   }
	   eq[++j].w = complex(re, im);

	   re  = sq_freq*ac[i].Rs;		/* series element */
	   im  =    freq*ac[i].Ls;
	   eq[++j].w = complex(re, im);
	}
}

/******
*	Function: elimination
*	Note	: an in place forward elimination proceduere to solve
*		  the linear algebric equation, s = wx.
******/
void	elimination(
	int		ns3,	/* = 2*ns+1, ou ns = number of sections	*/
	fd_linear_equation  eq[] )

{
	int	i;

	eq[1].w = cmplx_add_mul( complex(1., 0.), eq[0].w, eq[1].w );
	eq[1].s = cmplx_add_mul( eq[0].s, eq[0].w, eq[1].s );

	for(i=2; i<=ns3; i++)
	{  eq[i].w = cmplx_add_mul( eq[i-2].w, eq[i-1].w, eq[i].w );
	   eq[i].s = cmplx_add_mul( eq[i-1].s, eq[i-1].w, eq[i].s );
	}
}

/******
*	Function: substitution
*	Note	: a backward substitution proceduere to solve
*		  the linear algebric equation, s = wx.
******/
void	substitution(
	int		ns3,	/* =2*ns+1, ou ns = number of sections	*/
	fd_linear_equation  eq[] )

{
	int	i;

	for(i=ns3; i>=1; i--)
	   eq[i].x = cmplx_div (
	      cmplx_sub_mul( eq[i].s, eq[i-1].w, eq[i+1].x ), eq[i].w
				);

	eq[0].x = cmplx_div( cmplx_sub( eq[0].s, eq[1].x ), eq[0].w );
}

/*****
*	Function : read_rad
*	Note :	read normalized radiation impedance from the current
*		directory.
*****/
void	read_rad (void)
{
	FILE	*in;

	in = fopen(RADIMPpath, "rb");
	if (in == NULL)
	{
		printf("Could not read file [%s]",RADIMPpath);
		return;
	}

	fread(&dka, sizeof(float), 1, in);
	fread(&n_ka, sizeof(int), 1, in);
	rad_re = (float *) calloc(n_ka, sizeof(float));
	rad_im = (float *) calloc(n_ka, sizeof(float));
	fread(rad_re, sizeof(float), n_ka, in);
	fread(rad_im, sizeof(float), n_ka, in);
	fclose(in);
}

/*****
*	Function : radimp
*	Note :	returns complex radiation impedance at a given freq..
*****/

cmplx	radimp (float ka)
{
	float	re, im, x1, x2, y1, y2;
	int	k;

	if( ka >= dka*n_ka )
	{  re = rad_re[n_ka - 1];
	   im = rad_im[n_ka - 1];
	}
	else
	{  k = (int)(ka/dka);
	   if( k == 0 ) { x1 = 0; y1 = 0;}
	   if( k >= 1 ) { x1 = rad_re[k - 1]; y1 = rad_im[k - 1];}
	   x2 = rad_re[k];
	   y2 = rad_im[k];
	   re = x1 + (ka - k*dka)*(x2 - x1)/dka;
	   im = y1 + (ka - k*dka)*(y2 - y1)/dka;
	}
	return ( complex(re, im) );
}


/******************( Functions called from a main )***********************/

/*****
*	Function : vtf_ini
*	Note :	Compute acoustic element coefficients from which their
*		values are determined as a function of the frequency (Hz).
*****/

void	vtf_ini ( void )
{
	float	omega;
	int	i;
	float	pi = 3.141593f;

	nph2 = 2*nph; nph3 = nph2+1; nph4 = nph2+2;
	nbu2 = 2*nbu; nbu3 = nbu2+1; nbu4 = nbu2+2;
	nna2 = 2*nna; nna3 = nna2+1; nna4 = nna2+2;

/*** Coefficients for computing acoustic elements ***/

	omega = (float)(2.0*pi);

/* series registance due to visscus loss */
	Rs  = (float)(sqrt(0.5*pi*omega*ro*mu));

/* acoustic elements */
	La = (float)(omega*ro/2.0);	/* acoustic mass (La)		*/
	Ca = (float)(omega/(ro*c*c));	/* acoustic stiffness (1/Ca)	*/

/* parallel conductance due to heat conduction loss */
	G =  (float)((eta-1.)*sqrt(pi*omega*lamda/(ro*cp))/(ro*c*c));

/* walls */
	Rw = (float)(wall_resi/(2.0*sqrt(pi)));
	Lw = (float)(wall_mass*sqrt(pi));
	Cw = (float)(wall_comp*sqrt(pi));

/* radiation impedance */
	Rrad = (float)(2.0*pi*ro/c);
	Lrad = (float)(16.0*ro/(3.0*sqrt(pi)));

/* read normalized radiation impedance from the file */
	if( rad_boundary == BESSEL_FUNCTION ) read_rad();

/*** memory allocations ***/

	afph = (area_function *) calloc( nph, sizeof(area_function) );
	acph = (fd_acoustic_elements *) calloc( nph+1, sizeof(fd_acoustic_elements));
	eqph = (fd_linear_equation *) calloc( 2*nph+3, sizeof(fd_linear_equation) );

	afbu = (area_function *) calloc( nbu, sizeof(area_function) );
	acbu = (fd_acoustic_elements *) calloc( nbu+1, sizeof(fd_acoustic_elements) );
	eqbu = (fd_linear_equation *) calloc( 2*nbu+3, sizeof(fd_linear_equation) );

	acna = (fd_acoustic_elements *) calloc( nna+1, sizeof(fd_acoustic_elements) );
	eqna = (fd_linear_equation *) calloc( 2*nna+3, sizeof(fd_linear_equation) );

/**** Acoustic and matrix elements ****/

	copy_initial_af();		/* copy the VT area function */

/* pharyngeal tract */
	acou_elements( nph, afph, acph );

/* bucal cavity */
	acou_elements( nbu, afbu, acbu );

	lip_Rrad = Rrad;
	lip_Lrad = (float)(Lrad/sqrt(afbu[0].A));

	lip_rad_ka = (float)(2.0*sqrt(pi*afbu[0].A)/c);
	lip_rad_norm = ro*c/afbu[0].A;

/* nasal tract */
	if( nasal_tract == ON )
	{  afnt[nna-1].A = nonzero( anc );	/* nasal coupling */
	   acou_elements( nna, afnt, acna );
	   for(i=0; i<=nna; i++)
	      acna[i].G *= extra_loss_factor;	/* ad.hoc.extra loss */

	   nose_Rrad = Rrad;
	   nose_Lrad = (float)(Lrad/sqrt(afnt[0].A));

	   nose_rad_ka = (float)(2.0*sqrt(pi*afnt[0].A)/c);
	   nose_rad_norm = ro*c/afnt[0].A;
	}
}

/*****
*	Function: vtf_sim
*	Note :	compute vocal tract transfer, U/Ug, at a specified
*		frequency (in Hz) in the frequency domain simulation.
*		Returns the tranfer ratio.
*****/

float	vtf_sim( float freq )
{
	int	i, j;
	cmplx	f, g, h, p, q;
	cmplx	Uout;


/*** pharyngeal cavity ***/

	eq_elements( freq, nph, acph, eqph );

/* glottal boundary condition, open or close */
	if(source_loc == 0) glt_boundary = CLOSE;	/* by default */
	if(glt_boundary == CLOSE) eqph[0].w = complex( 0., 0.);
	else			  eqph[0].w = complex( 1.0e+8, -1.0e+8 );

/*** bucal cavity ***/

	eq_elements( freq, nbu, acbu, eqbu );

/* lip end boundary condition */
	switch( rad_boundary )
	{
	case BESSEL_FUNCTION :			/* circular piston set on */
	   eqbu[0].w = cmplx_inv (		/* an infinite buffle.	  */
	      cmplx_scl( lip_rad_norm, radimp(freq*lip_rad_ka ) ) );
	   break;

	case RL_CIRCUIT :			/* R-L circuit approximation */
	   eqbu[0].w = cmplx_inv (
	      complex( freq*freq*lip_Rrad, freq*lip_Lrad ) );
	   break;

	case SHORT_CIRCUIT :			/* nearly short circuit */
	   eqbu[0].w = complex(1.0e+8, -1.0e+8);
	   break;
	}

/*** Unit source inside the pharyngeal and bucal cavities ***/

	for(i=0; i<=nph3; i++) eqph[i].s = complex( 0., 0.);
	for(i=0; i<=nbu3; i++) eqbu[i].s = complex( 0., 0.);

	if(source_loc <= nph)
	{  if(source_loc == 0)  eqph[0].s = complex( 1., 0. );	/* flow */
	   else  eqph[2*source_loc+1].s = complex( 10., 0. );	/* pressure */
	}
	else
	{  j = (nph+nbu) - source_loc;
	   eqbu[2*j+1].s = complex(100., 0. );
	}

/*** nasal tract ***/

	if( nasal_tract == ON )
	{  eq_elements( freq, nna, acna, eqna );
	   for(i=0; i<=nna3; i++) eqna[i].s = complex( 0., 0.);

	/* nostrils boundary condition */

	   switch( rad_boundary )
	   {
	   case BESSEL_FUNCTION :
	      eqna[0].w = cmplx_inv (
		 cmplx_scl( nose_rad_norm, radimp(freq*nose_rad_ka) ) );
	      break;

	   case RL_CIRCUIT :
	      eqna[0].w = cmplx_inv (
		 complex( freq*freq*nose_Rrad, freq*nose_Lrad ) );
	      break;

	   case SHORT_CIRCUIT :
	      eqna[0].w = complex(1.0e+8, -1.0e+8);
	      break;
	   }
	}

/* solve s = wx */
	if( nasal_tract	== ON )
	{  elimination(nph3, eqph);
	   elimination(nbu3, eqbu);
	   elimination(nna3, eqna);

	   f = cmplx_div( eqph[nph3].s, eqph[nph3].w );
	   g = cmplx_div( eqbu[nbu3].s, eqbu[nbu3].w );
	   h = cmplx_div( eqna[nna3].s, eqna[nna3].w );
	   p = cmplx_add( f, cmplx_add( g, h ) );

	   f = cmplx_div( eqph[nph2].w, eqph[nph3].w );
	   g = cmplx_div( eqbu[nbu2].w, eqbu[nbu3].w );
	   h = cmplx_div( eqna[nna2].w, eqna[nna3].w );
	   q = cmplx_add( f, cmplx_add( g, h ) );

	   eqph[nph4].x = eqbu[nbu4].x = eqna[nna4].x = cmplx_div( p, q );

	 /*substitution(nph3, eqph);*/
	   substitution(nbu3, eqbu);
	   substitution(nna3, eqna);
	}
	else
	{  elimination(nph3, eqph);
	   elimination(nbu3, eqbu);

	   f = cmplx_div( eqph[nph3].s, eqph[nph3].w );
	   g = cmplx_div( eqbu[nbu3].s, eqbu[nbu3].w );
	   p = cmplx_add( f, g );

	   f = cmplx_div( eqph[nph2].w, eqph[nph3].w );
	   g = cmplx_div( eqbu[nbu2].w, eqbu[nbu3].w );
	   q = cmplx_add( f, g );

	   eqph[nph4].x = eqbu[nbu4].x = cmplx_div( p, q );

	 /*substitution(nph3, eqph);*/
	   substitution(nbu3, eqbu);
	}

/*** return the transfer ratio ***/

	Uout = eqbu[1].x;
	if( nasal_tract == ON ) Uout = cmplx_add( Uout, eqna[1].x );
	return( abs_cmplx( Uout) );
}

/*****
*	Function: vtf_UPdistribution
*	Note :	compute volume velocity and pressure distribution inside
*		the vocal tract for a given frequency (in Hz), in the
*		frequency domain simulation.
*****/

void	vtf_UPdistribution (
	float	freq,		/* frequency in Hz */
	float	*U,		/* volume velocity distribution */
	float	*P )		/* presure distribution */
{
	int	i, j;
	cmplx	f, g, h, p, q;



/*** pharyngeal cavity ***/

	eq_elements( freq, nph, acph, eqph );

/* glottal boundary condition, open or close */
	if(source_loc == 0) glt_boundary = CLOSE;	/* by default */
	if(glt_boundary == CLOSE) eqph[0].w = complex( 0., 0.);
	else			  eqph[0].w = complex( 1.0e+8, -1.0e+8 );

/*** bucal cavity ***/

	eq_elements( freq, nbu, acbu, eqbu );

/* lip end boundary condition */
	switch( rad_boundary )
	{
	case BESSEL_FUNCTION :			/* circular piston set on */
	   eqbu[0].w = cmplx_inv (		/* an infinite buffle.	  */
	      cmplx_scl( lip_rad_norm, radimp(freq*lip_rad_ka ) ) );
	   break;

	case RL_CIRCUIT :			/* R-L circuit approximation */
	   eqbu[0].w = cmplx_inv (
	      complex( freq*freq*lip_Rrad, freq*lip_Lrad ) );
	   break;

	case SHORT_CIRCUIT :			/* nearly short circuit */
	   eqbu[0].w = complex(1.0e+8, -1.0e+8);
	   break;
	}

/*** Unit source inside the pharyngeal and bucal cavities ***/

	for(i=0; i<=nph3; i++) eqph[i].s = complex( 0., 0.);
	for(i=0; i<=nbu3; i++) eqbu[i].s = complex( 0., 0.);

	if(source_loc <= nph)
	{  if(source_loc == 0)  eqph[0].s = complex( 1., 0. );	/* flow */
	   else  eqph[2*source_loc+1].s = complex( 10., 0. );	/* pressure */
	}
	else
	{  j = (nph+nbu) - source_loc;
	   eqbu[2*j+1].s = complex(10., 0. );
	}

/*** nasal tract ***/

	if( nasal_tract == ON )
	{  
		eq_elements( freq, nna, acna, eqna );
	   for(i=0; i<=nna3; i++) eqna[i].s = complex( 0., 0.);

	/* nostrils boundary condition */

	   switch( rad_boundary )
	   {
	   case BESSEL_FUNCTION :
	      eqna[0].w = cmplx_inv (
		 cmplx_scl( nose_rad_norm, radimp(freq*nose_rad_ka) ) );
	      break;

	   case RL_CIRCUIT :
	      eqna[0].w = cmplx_inv (
		 complex( freq*freq*nose_Rrad, freq*nose_Lrad ) );
	      break;

	   case SHORT_CIRCUIT :
	      eqna[0].w = complex(1.0e+8, -1.0e+8);
	      break;
	   }
	}

/* solve s = wx */
	if( nasal_tract	== ON )
	{  elimination(nph3, eqph);
	   elimination(nbu3, eqbu);
	   elimination(nna3, eqna);

	   f = cmplx_div( eqph[nph3].s, eqph[nph3].w );
	   g = cmplx_div( eqbu[nbu3].s, eqbu[nbu3].w );
	   h = cmplx_div( eqna[nna3].s, eqna[nna3].w );
	   p = cmplx_add( f, cmplx_add( g, h ) );

	   f = cmplx_div( eqph[nph2].w, eqph[nph3].w );
	   g = cmplx_div( eqbu[nbu2].w, eqbu[nbu3].w );
	   h = cmplx_div( eqna[nna2].w, eqna[nna3].w );
	   q = cmplx_add( f, cmplx_add( g, h ) );

	   eqph[nph4].x = eqbu[nbu4].x = eqna[nna4].x = cmplx_div( p, q );

	   substitution(nph3, eqph);
	   substitution(nbu3, eqbu);
	   substitution(nna3, eqna);
	}
	else
	{  elimination(nph3, eqph);
	   elimination(nbu3, eqbu);

	   f = cmplx_div( eqph[nph3].s, eqph[nph3].w );
	   g = cmplx_div( eqbu[nbu3].s, eqbu[nbu3].w );
	   p = cmplx_add( f, g );

	   f = cmplx_div( eqph[nph2].w, eqph[nph3].w );
	   g = cmplx_div( eqbu[nbu2].w, eqbu[nbu3].w );
	   q = cmplx_add( f, g );

	   eqph[nph4].x = eqbu[nbu4].x = cmplx_div( p, q );

	   substitution(nph3, eqph);
	   substitution(nbu3, eqbu);
	}

/*** distibutions ***/

	for(i=0; i<=nph; i++)
	{  U[i-1] = abs_cmplx( eqph[2*i+1].x );
	   P[i-1] = abs_cmplx( eqph[2*i].x );
	}
	for(i=0; i<=nbu; i++)
	{  U[nph+nbu-i] = abs_cmplx( eqbu[2*i-1].x );
	   P[nph+nbu-i] = abs_cmplx( eqbu[2*i].x );
	}
}

/*****
*	Function : vtf_term
*	Note :	free memories
****/

void	vtf_term ( void )
{
	free( afph );
	free( acph );
	free( eqph );

	free( afbu );
	free( acbu );
	free( eqbu );

	free( acna );
	free( eqna );

	if( rad_boundary == BESSEL_FUNCTION )
	{ 
	//	free( rad_re );
	//	free( rad_im );
	}
}



/******
*	Function : calplot_tf_FA
*	Note :	Compute Ul/Ug transfer ratio of a given VT area function
*		and plot.  Formant frrequencies and amplitudes are
*		estimated by a peak picking algorithm with a variable
*		frequency step size scheme.  The transfer is plotted
*		but not returned to calling program.
*****/
void	calplot_tf_FA (
	int	nfrmmax,	/* maximum number of F-peaks to be detected */
	float	*frm,		/* F1, F2, .... frequency in Hz */
	float	*amp,		/* A1, A2, .... amplitude in dB */
	int	*nfrms,		/* number of F-peaks detected   */
	float *tf,   /* transfer function output */
	int *ncount)	/* no of points */
{
	float	df, dfmin = .5f, dfmax = 40.f;			/* in Hz */
	float	freq, freq2, mag, mag1, mag2, mag3;
	float	fmin = 0.f, ftic = 0.2f, fmax = 5.f;		/* in kHz */
	int	i, j, count = 0;

#ifdef PLOTME
	def_vp (vp.n, vp.x, vp.y, vp.x+vp.w, vp.y+vp.h, 1);
	open_vp(vp.n);
	clear_vp(vp.n);
	g_frame ( vp.n,
		 1, title, "Frequency (kHz)","Relative mag. (dB)",
		 fmin,  40.,  fmax, -40.0,
		 0, 5, fmin,  ftic, fmax,
		 0, 1, -40.0,  10.0, 40.0 );
#endif

	j  = 0;
	df = dfmax;

	vtf_ini();
	for( freq=df, i=0; freq <= 1000.*fmax; freq=freq+df, i++ )
	{ 
		mag = (float)(20.*log10( max(vtf_sim(freq), 0.01) ));
		if( i == 0 )
		{ 
			tf[count] = mag;
			count++;
#ifdef PLOT_ME 
			mov( datolo_x(vp.n, freq/1000.), datolo_y(vp.n, mag) );
#endif
			mag1 = mag;
		}
		if( i == 1 )
		{ 
			tf[count] = mag;
			count++;
#ifdef PLOT_ME
			drw( datolo_x(vp.n, freq/1000.), datolo_y(vp.n, mag) );
#endif
			mag2 = mag;
			freq2 = freq;
		}
		if( i > 1 )
		{ 
			tf[count] = mag;
			count++;
#ifdef PLOT_ME
			drw( datolo_x(vp.n, freq/1000.), datolo_y(vp.n, mag) );
#endif
			mag3 = mag;
			/* detect spectral peak */
			if( mag1<mag2 && mag2>=mag3 && j<nfrmmax )
			{
				frm[j]   = freq2;		/* Fn in Hz    */
				amp[j++] = mag2;		/* An in dB    */
			}
			/* adapt frequency step size */
			if ( fabs(mag3 - mag2) > 1.f )
			{
				df = df - 8;
				if( df <= dfmin )
					df = dfmin;
			}
			else
			{
				df = df + 8;
				if( df >= dfmax ) 
					df = dfmax;
			}
			mag1 = mag2;
			mag2 = mag3;
			freq2 = freq;
		}
	}
	*nfrms = j;
	*ncount = count;
	vtf_term();
}


/******
*	Function : calplot_tf_FBA
*	Note :	Compute Ul/Ug transfer ratio of a given VT area function
*		and plot.  Formant frrequencies and amplitudes are
*		estimated by a peak picking algorithm with a variable
*		frequency step size scheme.  The transfer is plotted
*		but not returned to calling program.
*	Modif.: David Gesbert added bandwidth detection routine (16/09/92).
*	Modif.: Maeda modified it with a binary search argorithm (23/09/92).
*****/
void	calplot_tf_FBA (
	int	nfrmmax,	/* maximum number of F-peaks to be detected */
	float	*frm,		/* F1, F2, .... frequency in Hz */
	float	*bw,		/* B1, B2, .... frequency in Hz */
	float	*amp,		/* A1, A2, .... amplitude in dB */
	int	*nfrms,		/* number of F-peaks detected   */
	float *tf,   /* transfer function output */
	int *ncount)	/* no of points */
{
	float	df, dfmin = .5f, dfmax = 40.f;			/* in Hz */
	float	freq, freq2, mag, mag1, mag2, mag3;
	float	fmin = 0.f, ftic = 0.2f, fmax = 5.f;		/* in kHz */
	int	i, j,count = 0;

	/* variables used in David' algorithm */
/****
	int	change_direction, bw_not_found;
	float	pas, mag0;
****/
	/* variables used in binary search algorithm */

	float	dA = 0.5f;					/* in dB */
	float	A3dB, left_half_bw, right_half_bw;

/* Graphics frame intialization */

#ifdef PLOTME
	def_vp (vp.n, vp.x, vp.y, vp.x+vp.w, vp.y+vp.h, 1);
	open_vp(vp.n);
	clear_vp(vp.n);
	g_frame ( vp.n,
		 1, title, "Frequency (kHz)","Relative mag. (dB)",
		 fmin,  40.,  fmax, -40.0,
		 0, 5, fmin,  ftic, fmax,
		 0, 1, -40.0,  10.0, 40.0 );
#endif

/* Transfer ratio calculation */

	j  = 0;
	df = dfmax;

	vtf_ini();
	for( freq=df, i=0; freq <= 1000.*fmax; freq=freq+df, i++ )
	{
		mag = (float)(20.*log10( max(vtf_sim(freq), 0.01) ));
		if( i == 0 )
		{
			tf[count] = mag;
			count++;
#ifdef PLOTME
			mov( datolo_x(vp.n, freq/1000.), datolo_y(vp.n, mag) );
#endif
			mag1 = mag;
		}
		if( i == 1 )
		{
			tf[count] = mag;
			count++;
#ifdef PLOTME
			drw( datolo_x(vp.n, freq/1000.), datolo_y(vp.n, mag) );
#endif
			mag2 = mag;
			freq2 = freq;
		}
		if( i > 1 )
		{ 
			tf[count] = mag;
			count++;
#ifdef PLOTME
			drw( datolo_x(vp.n, freq/1000.), datolo_y(vp.n, mag) );
#endif
			mag3 = mag;
			/* detect spectral peak */
			if( mag1<mag2 && mag2>=mag3 && j<nfrmmax )
			{
				frm[j]   = freq2;		/* Fn in Hz    */
				amp[j++] = mag2;		/* An in dB    */
			}
			/* adapt frequency step size */	     
			if ( fabs(mag3 - mag2) > 1.f )
			{
				df = df - 8;
				if( df <= dfmin ) df = dfmin;
			}
			else
			{ 
				df = df + 8;
				if( df >= dfmax ) df = dfmax;
			}
			mag1 = mag2;
			mag2 = mag3;
			freq2 = freq;
		}
	}
	*nfrms = j;
	*ncount = count;

/****** -3db formant bandwidth detection (David Gesbert) ******/
/****
	for(j=0;j< *nfrms;j++)
	{
	   change_direction=0;
	   pas=-4*dfmin;
	   bw_not_found=1;

	   do
	   {
	      if(change_direction) pas=-pas;
	      freq=frm[j];
	      mag0=20*log10(max(vtf_sim(freq),0.01));
	      mag=mag0;

	      do
	      {
		 freq+=pas;
		 mag1=mag;
		 mag=20*log10(max(vtf_sim(freq),0.01));
		 if(mag<mag0-3) {bw_not_found=0; break;}
	      }
	      while(((mag1>=mag)&&(fabs(freq-frm[j])>20)) ||
		     ((mag<mag0+0.1) && (fabs(freq-frm[j])<=20)));

	      if (change_direction && bw_not_found) {freq=frm[j]; break;}
	      change_direction=1;
	   }
	   while(bw_not_found);

	   bw[j]=2*fabs(freq-frm[j]);
	}
****/

/*** -3dB bandwidth detection with binary search algorithm ***/

	for(j=0; j<*nfrms; j++)
	{
	   A3dB = (float)(amp[j] - 3.0f);


	/* -3dB point on the left skirt of the (j+1)-th formant */

	   if( j==0 )
		   df = (float)(frm[j]/2.);
	   else	
		   df = (float)((frm[j] - frm[j-1])/2.);

	   freq = frm[j] - df;
	   mag = (float)(20*log10(max(vtf_sim(freq),0.01)));
	   if( mag > A3dB ) 
		   left_half_bw = 0;
	   else
	   {  
		   for(i=0; i<20; i++)
		   {
			   df /= 2;
			   if( df<dfmin )
				   break;
			   mag = (float)(20*log10(max(vtf_sim(freq),0.01)));
			   
			   if( mag > A3dB+dA )	  freq -= df;
			   else if( mag < A3dB-dA ) freq += df;
			   else			  break;
		   }
		   left_half_bw = frm[j] - freq;
	   }

	/* -3dB point on the right skirt of the (j+1)-th formant */

	   if( j==*nfrms-1 ) 
		   df = (float)((1000.*fmax - frm[j])/2.);
	   else		   
		   df = (float)((frm[j+1] - frm[j])/2.);

	   freq = frm[j] + df;
	   mag = (float)(20*log10(max(vtf_sim(freq),0.01)));
	   if( mag > A3dB )
		   right_half_bw = 0;
	   else
	   {  
		   for(i=0; i<20; i++)
		   {
			   df /= 2;
			   if( df<dfmin ) break;

			   mag = (float)(20*log10(max(vtf_sim(freq),0.01)));
			   
			   if( mag > A3dB+dA )	  freq += df;
			   else if( mag < A3dB-dA ) freq -= df;
			   else			  break;
		   }
		   right_half_bw = freq - frm[j];
	   }
	   if( left_half_bw == 0. )  right_half_bw *= 2.;
	   if( right_half_bw == 0. ) left_half_bw  *= 2.;
	   bw[j] = left_half_bw + right_half_bw;
	}
	vtf_term();
}


/******
*	Function : calmultiplot_tf_FA
*	Note :	Compute Ul/Ug transfer ratio of a given VT area function
*		and plot.  Formant frequencies and amplitudes are
*		estimated by a peak picking algorithm with a variable
*		frequency step size scheme.  The transfer is plotted
*		but not returned to calling program.
*****/
void	calmultiplot_tf_FA (
	int	entry,		/* 0, 1, 2, ...	*/
	char	*main_title,
	char	*vowel_code,
	int	nfrmmax,	/* maximum number of F-peaks to be detected */
	float	*frm,		/* F1, F2, .... frequency in Hz */
	float	*amp,		/* A1, A2, .... amplitude in dB */
	int	*nfrms,		/* number of F-peaks detected */
	vp_lo_frame  vp )
{
	float	dfmin = .5f, dfmax = 40.f, df;		/* in Hz */
	float	frq2, freq, mag, mag1, mag2, mag3;
	float	fmin = 0.f, ftic = 0.2f, fmax = 5.f;	/* in kHz */
	int	i, j;

#ifdef PLOTME
	char	txt[80];

/* Initialize graphics frame */
	if( entry == 0 )
	{  def_vp (vp.n, vp.x, vp.y, vp.x+vp.w, vp.y+vp.h, 1);
	   open_vp(vp.n);
	   clear_vp(vp.n);
	   g_frame ( vp.n,
		 12, main_title, "Frequency (kHz)","Magnitude (dB)",
		 fmin,  40.,  fmax, -40.0,
		 0, 5, fmin,  ftic, fmax,
		 0, 1, -40.0,  10.0, 40.0 );
	}

	textstyle( 12, HORIZ_DIR);
	sprintf( txt, "Vowel=[%s]", vowel_code);
	plottextxy( vp.w-22., vp.h-18., txt );
	sprintf( txt, "NC   =%3.1f", anc );
	plottextxy( vp.w-22., vp.h-14., txt );

/* Set line style */
	if( entry == 0 ) setlinestyle( SOLID_LINE, 0, NORM_WIDTH );
	if( entry >= 1 ) setlinestyle( DOTTED_LINE, 0, NORM_WIDTH );
#endif

/* Compute and plot vocal tract transfer ratio */
	j  = 0;
	df = dfmax;

	vtf_ini();
	for( freq=df, i=0; freq <= 1000.*fmax; freq=freq+df, i++ )
	{ 
		mag = (float)(20.*log10( max(vtf_sim(freq), 0.01) ));
		if( i == 0 )
		{
#ifdef PLOTME
			mov( datolo_x(vp.n, freq/1000.), datolo_y(vp.n, mag) );
#endif
			mag1 = mag;
		}
		if( i == 1 )
		{
#ifdef PLOTME
			drw( datolo_x(vp.n, freq/1000.), datolo_y(vp.n, mag) );
#endif
			mag2 = mag;
			frq2 = freq;
		}
		if( i > 1 )
		{
#ifdef PLOTME
			drw( datolo_x(vp.n, freq/1000.), datolo_y(vp.n, mag) );
#endif
			mag3 = mag;
			/* detect spectral peak */
			if( mag1<mag2 && mag2>=mag3 && j<nfrmmax )
			{
				frm[j]   = frq2;		/* Fn in Hz    */
				amp[j++] = mag2;		/* An in dB    */
			}
			/* adapt frequency step size */
			if ( fabs(mag3 - mag2) > 1.f )
			{
				df = df - 8;
				if( df <= dfmin ) df = dfmin;
			}
			else
			{
				df = df + 8;
				if( df >= dfmax ) df = dfmax;
			}
			mag1 = mag2;
			mag2 = mag3;
			frq2 = freq;
		}
	}
	*nfrms = j;
	vtf_term();
#ifdef PLOTME
	setlinestyle( SOLID_LINE, 0, NORM_WIDTH );
#endif
}
