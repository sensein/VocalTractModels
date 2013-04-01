/*******
*	File :	plot_lib.h
*	Note :	Definition of structures and prototypes.
********/

typedef struct{ float A, x; } area_function;

typedef struct{ int	n;			/* viewport number	*/
		float	x,			/* left-top position in	*/
			y,			/* logical units	*/
			w,			/* width and height	*/
			h;
		} vp_lo_frame;

void	read_af(char *af_path, int *nss, area_function **af);

/*
// All this gets converted to Matlab
void	plot_af	(int ns, area_function *af, vp_lo_frame vp);
void	plot_af_source_loc
		(char *title, int src_loc, int ns, area_function *af,
		 vp_lo_frame vp );
void	plot_signal
		(int n, float *t, float *s, vp_lo_frame vp);
void	plot_tf	(int np, float *f, float *tf, vp_lo_frame vp);
void	plot_FA	(int nfrms, float *frm, float *amp,
		 vp_lo_frame vp);
void	plot_FBA(int nfrms, float *frm, float *bw, float *amp,
		 vp_lo_frame vp);
void	plot_PUR (float Psub, float Udc, float Rg_v, float Rg_k,
		  float Rc_k, vp_lo_frame vp);
void	plot_glottal_area
		(int n, float tmin, float tmax,
		 float *t, float *s, vp_lo_frame vp );
void	plot_speech_signal
		(int n, float tmin, float tmax,
		 float *t, float *s, vp_lo_frame vp );
*/
