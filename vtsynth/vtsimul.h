

/******
*	File :	vtsimul.h
*	Note :	An include file for acoustic simulation of the vocal
*		tract and articulatory model, to be included in a simulation
*		library.
******/

/*******************( Macros : option parameter values )******************/
#pragma warning (disable : 4244 4101 4305)

#define	OFF		0
#define	ON		1

#define AREA_MIN	0.001	/* minimum cross_sectional area cm**2 */

#define	TRACHEA		0	/* tract (VT section) correpondence */
#define	GLOTTIS		1
#define	PHARYNX		2
#define	MOUTH		3
#define	NOSE		4	/* fixed nasal tract	*/
#define	VOCAL_TRACT	5	/* PHARYNX + MOUTH = 5	*/
#define	N_PORT		6	/* nasal coupling section */

#define	RIGID		0	/* wall type		*/
#define	YIELDING	1
#define	FLOW		0	/* source type		*/
#define	PRESSURE	1
#define CLOSE		0	/* glottal condition	*/
#define OPEN		1
#define	SHORT_CIRCUIT	0	/* radiation load type	*/
#define RL_CIRCUIT	1
#define	BESSEL_FUNCTION	2

#define	TIME_VARYING	1	/* VT time property	*/
#define STATIONARY	0

#define	UNIFORM_TUBE	0	/* VT area function specification type */
#define	FROM_FILE	1
#define T2MODEL		2
#define	P3MODEL		3
#define	ART_MODEL	4

/***********************( structure definitions )************************/

#if !defined(AREA_FUNCTION)
#define AREA_FUNCTION area_function
typedef struct{float A, x;}AREA_FUNCTION;
#endif

typedef	struct{	
	int		Loc,	/* =0, in phrx+mouth, =1 mouth, =2 nose */
			N;	/* section number, noise source at the exit */
	float	A,	/* cross_sectional area (cm2)		*/
			x;	/* length (cm)				*/
}	constriction;

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

/*
void	calplot_tf_FA
(char *title, int frmmax, float *frm, float *amp,
int *nfrms, vp_lo_frame vp);
void	calplot_tf_FBA
(char *title, int frmmax, float *frm, float *bw, float *amp,
int *nfrms, vp_lo_frame vp );
void	calmultiplot_tf_FA
(int entry, char *title, char *vowel_code,
int frmmax, float *frm, float *amp, int *nfrms, vp_lo_frame vp);
*/

/* Time domain calculation */

int	vtt_ini ( void );
float	vtt_sim ( void );
void	vtt_term( void );
