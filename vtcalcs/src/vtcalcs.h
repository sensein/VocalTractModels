#ifndef _VTCALCS_H
#define _VTCALCS_H

/**********( global variables seen by functions in "vtcalcs.c" )*********/

/* default values for the models */
	Utube_par	UTpar = {6.0, 17.5};	/* uniforme tube	*/
	int		nss_uni = 17;		/* # of sections	*/
	char		AFfilename[9] = "iy";	/* area function file	*/
	T2model_par	T2par = {1.5, 8.5, 8.f, 8.5};	/* 2 tube model	*/
	P3model_par	P3par = {2.6f, 8, 4.f};	/* 3 paramter model	*/
	float		AMpar[7]		/* articulatory model	*/
			   ={0., 0., 0., 0., 0., 0., 0.};

	char	save_file_code[8] = "VOWEL";

	int	save_file_count   = 0;
	int	sig_save_flag     = OFF;
	FILE	*sig_flux, *glt_flux;

/* Footnotes */
	char	*note_calc =
	   "Calculate the vocal tract transfer and formant frequencies";
	char	*note_synt =
	   "Synthesize a stationary vowel using the specified area function";
	char	*note_save =
	   "Save model parameters, area function, etc. on .DOC, .SIG files";
	char	*note_vowel =
	   "Enter name; iy, ey, eh, ah, aa, ao, oh, uw, iw, ew, oe, ...etc.";

/*********************( grobal plotting variables )************************/

	int	nfrms, nfrmmax = 8;	/* number of formant frequencies  */
	float	frm[8], bw[8], amp[8];	/* formant freqs. and amplitude   */
	int ntf;
	float tfmag[8192];
	float tffreq[8192];

	vp_lo_frame vp1 = { 1, 0.5, 5.0, 49., 25. }; 	/* area_function  */
	vp_lo_frame vp2 = { 2, 0.5, 32., 49., 37. };    /* VT transfer	  */
	vp_lo_frame vp3 = { 3, 50., 42., 29., 23. };	/* formants	  */
	vp_lo_frame vp4 = { 4, 55., 5.,  35., 35. };	/* VT profile	  */

	char	*UmUg_title = "(Um + Un)/Ug TRANSFER RATIO";
	char	*UmPt_title = "(Um + Un)/Pt TRANSFER RATIO";

/****************************( global flags )******************************/

/**	static	int	playback_flag = ON;  **/

extern float	*rad_re, *rad_im;

#define max(a,b) ((a>b)?a:b)

#define min(a,b) ((a<b)?a:b)

#endif
