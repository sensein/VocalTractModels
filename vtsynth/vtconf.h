/******
*	File :	vtconf.h
*	Note :	The second version of "vtconfig.h".
*		Acoustic and geometrical configurations for the vocal
*		tract calculations.  It should be declared in a "main".
*
*		To be used with new programs. (18/3/4)
******/

#pragma warning (disable : 4244 4101 4305)
#include	"vtsimul.h"

/******************( Global variables defined in "main" )****************/


/* Specific to time domain calculations */

float	simfrq = 30000.0;	/* simulation frequency in Hz	*/
float	smpfrq = 10000.0;	/* sampling freq. in Hz		*/

/*( Used in the time and frequency domain calculations )*/

float	Psub = 8.;		/* subglottal air pressure in cmH2O	*/
float	lung_vol = 4000;	/* lung volumes in cm3			*/
int	ntr = 20;		/* # of tracheal sections		*/
area_function	aftr[20] = {
	{3.,1.},{3.,1.},{3.,1.},{3.,1.},{3.,1.},
	{3.,1.},{3.,1.},{3.,1.},{3.,1.},{3.,1.},
	{3.,1.},{3.,1.},{3.,1.},{3.,1.},{3.,1.},
	{3.,1.},{3.,1.},{3.,1.},{3.,1.},{3.,1.} };

	float	Ag = 0.0;		/* grottis area in cm**2		*/
	float	xg = 0.3;		/* glottis thickness in cm (=0.3 cm)	*/
	float	lg = 1.2;		/* fold length in cm			*/

	int	nvt = 17;		/* # of vocal tract sections		*/
	area_function	*afvt;		/* AF from the glottis to the lips	*/

	int	nbp = 9;	/* nasal branch point (= # of pharynx sections) */

	float	Anc = 0.;		/* area of nasal coupling section in cm2 */
	int	nnc = 3;		/* # of coupling sections		*/
	int	nnt;			/* # of sections in entire nasal tract	*/
	area_function	*afnt;		/* AF from coupling point to nostrils	*/
	/* Note : areas the first ncc sections	*/
	/*	vary as a function of Anc.	*/
	/* side cavity fh = about 500 Hz (Helmholtz reso*/
	float	sinus_vol = 25.;	/* volume of the nasal side cavity (cm3)*/
	float	sinus_apt = 0.1;	/* aperture of the sinus canal (cm2)	*/
	float	sinus_len = 0.5;	/* length of the sinus canal (cm)	*/
	int	sinus_pos = 7;		/* sinus position in sections from nbp	*/

	/*********************( source related variables )***********************õ

	/** specific to the frequecy-domain calculation **/

	int	source_loc = 0;		/* source location in VT section number	*/
	/* Note: source_loc = 0 corresponds to the glottis inlet.	*/
	int	source_typ = FLOW;		/* or PRESSURE			*/

	/** specific to time_domain calculations **/

	int	vocal_tract  = STATIONARY;	/* or TIME_VARYING		*/
	int	dynamic_term = OFF;		/* or ON			*/

	int	noise_source = OFF;		/* or ON			*/
	float	noiseAmp = 100.;		/* amp. bound for noise source	*/
	int	noiseSourceLoc = 1;		/* 1 section downstream from Ac */

	int	nAcc;
	/* note: consonantal constrition section address to be used
	in the function below for Ac plot. The addess is determined
	as the minimum area by checking the area function of the
	consonant.  The simulation program (vtt_sim) also detects the
	minimum area to determine the source location, but this is done
	for evry computatioan cycles.  The detected position thus can
	vary along the tract during a VCV. nAcc is used to return Uac
	at the fixed section.
	*/


	/************************( simulation options )**************************/

	int	nasal_tract  = OFF;		/* or ON			*/
	int	sinus_cavity = OFF;		/* or ON			*/
	int	subGLTsystem = OFF;		/* or ON			*/
	int	wall         = RIGID;		/* or RIGID (YIELDING)		*/
	int	rad_boundary = RL_CIRCUIT;	/* SHORT_CIRCUIT, or BESSEL_FUN	*/
	int	wall_radiation= OFF;		/* radiation from VT walls	*/
	int	glt_boundary = CLOSE;		/* or OPEN			*/

	/************( an extra heat loss factor for the nasal tract )***********/

	float	extra_loss_factor = 50.;

	/************************( physical constants )**************************/

	float	ro    = 1.14e-3;	/* air density, gm/cm**3		*/
	float	c     = 3.5e+4;		/* sound velocity, cm/s			*/
	float	eta   = 1.4;		/* adiabatic constant			*/
	float	cp    = 0.24;		/* specific heat, cal/gm.degree		*/
	float	lamda = 5.5e-5;		/* heat conduction, cal/cm.sec.degree	*/
	float	mu    = 1.86e-4;	/* viscosity coef, dyne.sec/cm**2	*/
	float	wall_resi = 1600.;	/* wall resistance, gm/s/cm2		*/
	float	wall_mass = 1.5;	/* wall mass per area, gm/cm2		*/
	float	wall_comp = 3.0e+5;	/* wall compliance			*/
	float	H2O_bar= 980.39;
	float	Kc = 1.42;		/* Ishizaka's for turbulent flow	*/
	/* =.875, Van der Berg's constant	*/

	float	trachea_resi = 500.;	/* trachea wall resistance per area	*/
	float	trachea_mass = 0.3;	/* trachea wall mass per area, gm/cm2	*/
	float	trachea_comp = 1.0e+6;	/* trachea wall compliance		*/
