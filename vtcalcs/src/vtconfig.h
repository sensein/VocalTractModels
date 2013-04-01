

/******
*	File :	vtconfig.h
*	Note :	Acoustic and geometrical configurations for the vocal
*		tract calculations.  It should be declared in a "main".
******/

#include	<vtsimul.h>

/******************( Global variables defined in "main" )****************/


	/* Specific to time domain calculations */

float	simfrq = 30000.0f;	/* simulation frequency in Hz	*/
float	smpfrq = 10000.0f;	/* sampling freq. in Hz		*/

	/*( Used in the time and frequency domain calculations )*/

float	Psub = 8.f;		/* subglottal air pressure in cmH2O	*/
float	Ag = 0.0f;		/* grottis area in cm**2		*/
float	xg = 0.3f;		/* glottis thickness in cm (=0.3 cm)	*/
float	lg = 1.2f;		/* fold length in cm			*/
float	Kc = 1.42f;		/* Ishizaka's for turbulent flow	*/
				/* =.875, Van der Berg's constant	*/

int	nbu = 8;		/* # of sections in the bucal tube	*/
int	nph = 9;		/* # of sections in the phryngeal tube	*/
int	nvt;			/* = nbu + nph, # of VT sections	*/
area_function	*afvt;		/* AF from the glottis to the lips	*/

float	anc = 0.f;		/* nasal coupling area in cm2		*/
int	nna;			/* # of sections in the nasal tract	*/
area_function	*afnt;		/* AF from nostrils to coupling point	*/

int	nss;			/* total # of sections in VT system	*/

/************************( simulation options )**************************/

int	nasal_tract  = OFF;		/* or ON			*/
int	wall         = YIELDING;	/* or RIGID			*/
int	rad_boundary = RL_CIRCUIT;	/* SHORT_CIRCUIT, or BESSEL_FUN	*/
int	glt_boundary = CLOSE;		/* or OPEN			*/

int	source_loc = 0;		/* source location in VT section number	*/
int	source_typ = FLOW;		/* or PRESSURE			*/

	/* specific to time_domain calculations */

int	vocal_tract  = STATIONARY;	/* or TIME_VARYING		*/
int	dynamic_term = OFF;		/* or OFF			*/

/************( an extra heat loss factor for the nasal tract )***********/

float	extra_loss_factor = 50.f;

/************************( physical constants )**************************/

float	ro    = 1.14e-3f;	/* air density, gm/cm**3		*/
float	c     = 3.5e+4f;		/* sound velocity, cm/s			*/
float	eta   = 1.4f;		/* adiabatic constant			*/
float	cp    = 0.24f;		/* specific heat, cal/gm.degree		*/
float	lamda = 5.5e-5f;		/* heat conduction, cal/cm.sec.degree	*/
float	mu    = 1.86e-4f;	/* viscosity coef, dyne.sec/cm**2	*/
float	wall_resi = 1600.f;	/* wall resistance, gm/s/cm2		*/
float	wall_mass = 1.5f;	/* wall mass per area, gm/cm2		*/
float	wall_comp = 3.0e+5f;	/* wall compliance			*/
float	H2O_bar= 980.39f;
