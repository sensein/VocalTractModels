

/*****
*	File :	vsyn_lib.h
*	Note :	Prototypes etc..
*****/

#define	DACscale 4096.0		/* scale-up signal to fit 16 bits integer */
#define gltDACscale 100.*DACscale	/* sacle up for glottal source */

float	glottal_area( char model, char mode, float Ap, int *t0 );
void	vowel_synthesis( FILE *sig_file, FILE *glt_file, int buflen );

int	vtt_ini( void );
float	vtt_sim( void );
void	vtt_term( void );

