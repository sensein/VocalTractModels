

/******
*	File :	vtsimul.h
*	Note :	An include file for acoustic simulation of the vocal
*		tract and articulatory model, to be included in a simulation
*		library.
******/

/*******************( Macros : option parameter values )******************/

#define	OFF		0x00
#define	ON		0x01
#define	RIGID		0
#define	YIELDING	1
#define	FLOW		0
#define	PRESSURE	1
#define CLOSE		0
#define OPEN		1
#define	SHORT_CIRCUIT	0
#define RL_CIRCUIT	1
#define	BESSEL_FUNCTION	2

#define	TIME_VARYING	1
#define STATIONARY	0

#define	UNIFORM_TUBE	0	/* VT area function specification type */
#define	FROM_FILE	1
#define T2MODEL		2
#define	P3MODEL		3
#define	ART_MODEL	4

/***********************( structure definitions )************************/

typedef	struct{ 			/* Uniform acoustic tube	*/
		float	area;		/* cross-sectional area, cm2	*/
		float	length;		/* length, cm			*/
	      }	Utube_par;

typedef	struct{				/* Two tube model		*/
		float	A1;		/* A1: 1-st section area, cm2	*/
		float	x1;		/* x1: 1-st section length, cm	*/
		float	A2;		/* A2: 2-nd section area, cm2	*/
		float	x2;		/* x2: 2-nd section length, cm	*/
	      } T2model_par;
typedef	struct{				/* Three parameter model	*/
		float	At;		/* Constriction area, cm2	*/
		int	Xt;		/* Constriction location, cm	*/
		float	Al;		/* Lip opening apperture, cm2	*/
	      } P3model_par;

/*****************************( prototypes )*****************************/

/* Frequency domain calculation */

void	vtf_ini ( void );
float	vtf_sim ( float freq );
void	vtf_term( void );

#ifdef PLOTME
void	calplot_tf_FA
	   (char *title, int frmmax, float *frm, float *amp,
	    int *nfrms, vp_lo_frame vp);
void	calplot_tf_FBA
	   (int nfrmmax, float *frm, float *bw, float *amp,int *nfrms,float *tfmag,float *tffreq,int *ncount);


void	calmultiplot_tf_FA
	   (int entry, char *title, char *vowel_code,
	    int frmmax, float *frm, float *amp, int *nfrms, vp_lo_frame vp);
#endif

void	calplot_tf_FA
	   (int frmmax, float *frm, float *amp,int *nfrms,float *tf,int *ncount);
void	calplot_tf_FBA
	   (int nfrmmax, float *frm, float *bw, float *amp,int *nfrms,float *tfmag,float *tffreq,int *ncount);

void	calmultiplot_tf_FA
	   (int entry, char *title, char *vowel_code,
	    int frmmax, float *frm, float *amp, int *nfrms, vp_lo_frame vp);

/* Time domain calculation */

int	vtt_ini ( void );
float	vtt_sim ( void );
void	vtt_term( void );
