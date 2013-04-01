/***************************************************************************
*									   *
*	File :	plot_lib.c						   *
*	Note :	plotting functions for vocal tract calculation programs	   *
*									   *
***************************************************************************/

#include	<always.h>

#include	<plot_lib.h>
//#include	<grafix.h>


/*****
*	Function: read_af
*	Note :	Read area function data (for vocal trac and nasal tract )
*		from a text file.
*****/
void	read_af (
	char	*af_path,		/* complete file specification	*/
	int	*inss,			/* # of sections		*/
	area_function	**af )		/* area function		*/
{
	FILE	*in;
	area_function	*p;
	char	af_title[6];
	float	a, x;

	in = fopen(af_path, "rt");
	if (NULL == in)
	{
		printf("Could not read file [%s]",af_path);
		return;
	}

	fscanf(in, "%s", af_title );
	fscanf(in, "%d", inss);
	*af = (area_function *) calloc( *inss, sizeof(area_function) );
	fscanf(in, "%f", &x);
	for(p = *af; p<*af+*inss; p++)
	{  fscanf(in, "%f", &a);
	   p->A = a;
	   p->x = x;
	}
	fclose ( in );
}


#ifdef PLOTME
/*****
*       Function : plot_af
*	Note	 : Plot the area function of the vocal tract on a predefined
*		   frame and a viewport ( vp.n ).
*****/

void	plot_af ( int ns, area_function *af,  vp_lo_frame vp )
{
	int	i;
	float	*A_plt, *x_plt, x;

	def_vp(vp.n, vp.x, vp.y, vp.x + vp.w, vp.y + vp.h, 1);
	open_vp(vp.n);
	clear_vp(vp.n);
	g_frame ( vp.n,
		  1, "AREA FUNCTION", "distance from glottis (cm)",
		  "Area (cm~2)",
		  0.,  15.,  19., 0.0,
		  0, 5, 0.0,  1.0, 19.,
		  0, 2, 0.0,  1.0, 15.0 );

	A_plt = (float *) calloc(2*ns, sizeof(float));
	x_plt = (float *) calloc(2*ns, sizeof(float));
	x = 0;
	for( i=0; i<ns; i++)
	{  A_plt[2*i] = (float) max( 0., af[i].A);
	   x_plt[2*i] = x;
	   x = x + af[i].x;
	   A_plt[2*i + 1] = (float) max(0., af[i].A);
	   x_plt[2*i + 1] = x;
	}
	plotdata_xy(vp.n, 2*ns, x_plt, A_plt);
	free ( A_plt );
	free ( x_plt );
}

/*****
*       Function : plot_af_source_loc
*	Note	 : Plot the area function of the vocal tract and a source
*		   location on a predefined frame and a viewport ( vp.n ).
*****/

void	plot_af_source_loc (
	char		*title,		/* main title			*/
	int		srcloc,		/* source location in section #	*/
	int		ns,		/* # of sections		*/
	area_function	*af,		/* area function		*/
	vp_lo_frame	vp )		/* viewport specification	*/
{
	int	i;
	float	*A_plt, *x_plt, x;

	def_vp(vp.n, vp.x, vp.y, vp.x + vp.w, vp.y + vp.h, 1);
	open_vp(vp.n);
	clear_vp(vp.n);
	g_frame ( vp.n,
		  1, title, "Distance from glottis (cm)",
		  "Area (cm~2)",
		  0.,  15.,  19., 0.0,
		  0, 5, 0.0,  1.0, 19.,
		  0, 2, 0.0,  1.0, 15.0 );

	A_plt = (float *) calloc(2*ns, sizeof(float));
	x_plt = (float *) calloc(2*ns, sizeof(float));
	x = 0;
	for( i=0; i<ns; i++)
	{  A_plt[2*i] = (float) max( 0., af[i].A);
	   x_plt[2*i] = x;
	   x = x + af[i].x;
	   A_plt[2*i + 1] = (float) max(0., af[i].A);
	   x_plt[2*i + 1] = x;
	}
	plotdata_xy(vp.n, 2*ns, x_plt, A_plt);
/* put a mark at the source location */
	i = max( 0, 2*srcloc - 1);
	diamond_lo( datolo_x( vp.n, x_plt[i] ),
		    datolo_y( vp.n, A_plt[i]/2.), 1.0 );
	free ( A_plt );
	free ( x_plt );
}

/*****
*	Function : plot_signal
*	Note :	Plot signal on a fixed format graphics frame, i.e.,
*		time (ms) vs. relative pressure.
*****/
void	plot_signal (
	int	n,		/* number of samples */
	float	*t,		/* time coordinate value in ms */
	float	*s,		/* signal */
	vp_lo_frame vp )	/* frame in logical units */
{
	float	tmax;

	tmax = t[n-1];
	def_vp(vp.n, vp.x, vp.y, vp.x+vp.w, vp.y+vp.h, 1);
	clear_vp( vp.n );
	open_vp(vp.n);
	g_frame (vp.n,
	   1, "Time response", "time (ms)", "pressure (relative)",
	   0., 20., tmax, -20.,
	   0, 1, 0., 10., tmax,
	   0, 1, -20., 5., 20.);
	plotdata_xy(vp.n, n, t, s);
}

/*****
*	function : plot_tf
*	Note :	Plot VT transfer function
****/
void	plot_tf (
	int	n,		/* number of spectral samples */
	float	*f,		/* frequencies in Hz */
	float	*tf,		/* transfer magnitude in dB */
	vp_lo_frame vp )	/* viewport spec */
{
	float	fmax;

	fmax = f[n-1];
	def_vp(vp.n, vp.x, vp.y, vp.x+vp.w, vp.y+vp.h, 1);
	open_vp(vp.n);
	clear_vp(vp.n);
	g_frame (vp.n,
	   1, "Vocal tract transfer (U/Ug)", "frequency (Hz)",
	      "magnitude (dB)",
	   0., 40., fmax, -40.0,
	   0, 1, 0., 1000., fmax,
	   0, 1, -40., 10., 40.);
	plotdata_xy(vp.n, n, f, tf);
}

/******
*	Function : plot_FA
*	Note :	Plot formant Frequencie and Amplitude values
*		in graphics mode.
******/
void	plot_FA (
	int	nfrms,		/* number of formants */
	float	*frm,		/* formant frequencies in Hz */
	float	*amp,		/* formant amplitudes in dB */
	vp_lo_frame vp )	/* view port spec. */
{
	int	i;
	char	txt[40];

	def_vp(vp.n, vp.x, vp.y, vp.x + vp.w, vp.y + vp.h, 1);
	open_vp(vp.n);
	clear_vp(vp.n);
	plot_vp(vp.n);
	sprintf(txt, "  Frq(Hz) A(dB)");
	plottextxy( 1., 1., txt );
	for(i=0; i<nfrms; i++)
	{  sprintf(txt, "F%1d %5d %5d",
			i+1, (int)(frm[i]+.5), (int)(amp[i]+.5) );
	   plottextxy( 1., 2.5*i + 4., txt);
	}
}


/******
*	Function : plot_FBA
*	Note :	Plot formant Frequencie, Bandwidth and Amplitude values
*		in graphics mode.
******/
void	plot_FBA (
	int	nfrms,		/* number of formants */
	float	*frm,		/* formant frequencies in Hz */
	float	*bw,		/* formant bandwidth in Hz */
	float	*amp,		/* formant amplitudes in dB */
	vp_lo_frame vp )	/* view port spec. */
{
	int	i;
	char	txt[40];

	def_vp(vp.n, vp.x, vp.y, vp.x + vp.w, vp.y + vp.h, 1);
	open_vp(vp.n);
	clear_vp(vp.n);
	plot_vp(vp.n);
	sprintf(txt, "  Frq(Hz) Bw(Hz) A(dB)");
	plottextxy( 1., 1., txt );
	for(i=0; i<nfrms; i++)
	{  if( bw[i] > 0. )
	      sprintf( txt, "F%1d %5d %5d %5d", i+1,
		      (int)(frm[i]+.5), (int)(bw[i]+.5), (int)(amp[i]+.5) );
	   else
	      sprintf( txt, "F%1d %5d     * %5d", i+1,
		      (int)(frm[i]+.5), (int)(amp[i]+.5) );
	   plottextxy( 1., 2.5*i + 4., txt);
	}
}


/******
*	Function : plot_PUR
*	Note :	Plot DC air flow, and resistances in graphics mode.
******/
void	plot_PUR (
	float	Psub,		/* Subglottal air pressure in cmH2O */
	float	Udc,		/* air flow in cm3/s */
	float	Rg_v,		/* in acoustic ohm */
	float	Rg_k,		/* in acoustic ohm */
	float	Rc_k,		/* in acoustic ohm */
	vp_lo_frame vp )	/* view port spec. */
{
	char	txt[40];

	def_vp(vp.n, vp.x, vp.y, vp.x + vp.w, vp.y + vp.h, 1);
	open_vp(vp.n);
	clear_vp(vp.n);
	plot_vp(vp.n);

	sprintf( txt, "Psub = %8.1f cmH2O", Psub );
	plottextxy ( 1., 1.0, txt);
	sprintf( txt, "Udc  = %8.1f cm3/s", Udc );
	plottextRxy( 0., 2.5, txt);
	sprintf( txt, "Rg_v = %8.1f Ohms", Rg_v );
	plottextRxy( 0., 2.5, txt);
	sprintf( txt, "Rg_k = %8.1f Ohms", 2.*Rg_k );
	plottextRxy( 0., 2.5, txt);
	sprintf( txt, "Rc_k = %8.1f Ohms", 2.*Rc_k );
	plottextRxy( 0., 2.5, txt);

}

/*****
*	Function : plot_glottal_area
*	Note :	Plot glottal area variation on a fixed format graphics
*		frame, i.e., time (ms) vs. area (cm2).
*****/
void	plot_glottal_area (
	int	n,		/* number of samples */
	float	tmin,		/* left time coordinate in ms */
	float	tmax,		/* right time coordinate in ms */
	float	*t,		/* time coordinate value in ms */
	float	*s,		/* glottal area (cm2) */
	vp_lo_frame vp )	/* frame in logical units */
{
	def_vp(vp.n, vp.x, vp.y, vp.x + vp.w, vp.y + vp.h, 1);
	clear_vp( vp.n );
	open_vp( vp.n );
	g_frame (vp.n,
	   1, "Glottal area", "time (ms)", "area (cm2)",
	   tmin, 0.3, tmax, 0.,
	   0, 1, tmin, 10., tmax,
	   0, 1, 0.0, 0.1, 0.3);
	plotdata_xy(vp.n, n, t, s);
}


/*****
*	Function : plot_speech_signal
*	Note :	Plot speech signal on a fixed format graphics
*		frame, i.e., time (ms) vs. relative pressure.
*****/
void	plot_speech_signal (
	int	n,		/* number of samples */
	float	tmin,		/* left time coordinate in ms */
	float	tmax,		/* right time coordinate in ms */
	float	*t,		/* time coordinate value in ms */
	float	*s,		/* glottal area (cm2) */
	vp_lo_frame vp )	/* frame in logical units */
{
	def_vp(vp.n, vp.x, vp.y, vp.x + vp.w, vp.y + vp.h, 1);
	clear_vp( vp.n );
	open_vp( vp.n );
	g_frame (vp.n,
	   1, "Speech signal", "time (ms)", "pressure (relative)",
	   tmin, 3., tmax, -3.,
	   0, 1, tmin, 10., tmax,
	   0, 1, -3., 1., 3.);
	plotdata_xy(vp.n, n, t, s);
}

#endif
