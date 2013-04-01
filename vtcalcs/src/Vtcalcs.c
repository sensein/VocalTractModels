


/***************************************************************************
*                                                                          *
*	File :	vtcalcs.c                                                  *
*	Author : S. Maeda						   *
*	Note :	Various calculations of the vocal tract acoustics.	   *
*		Tract configurations, physical variables, and VT area	   *
*		function are modified by menu operations.  The modifiable  *
*		items are :						   *
*		1) Physical properties of air; density, sound velocity.	   *
*		2) Mechanical properties of tract walls ; mass etc..	   *
*		3) Simulation configurations ; nasal tract on/off,	   *
*		   radiation impedance representation, current/pressure	   *
*		   source.						   *
*		4) Area function specification ;			   *
*		   a) Uniform tube.					   *
*		   b) Arbitrary shape specified in a .ARE file.		   *
*		   c) 2 tube model					   *
*		   d) Fant's 3 parameter model.				   *
*		   e) Articulatory model with 7 parameters.		   *
*									   *
*		By menu commands, the program computes the VT transfer	   *
*		ratio and formant frequencies (Hz), bandwidth (Hz), and	   *
*		amplitudes (dB), or synthesized the corresponding	   *
*		stationary vowel (or nasalized vowel).  A frequency	   *
*		domain and a time domain simulation method are used	   *
*		respectively in the transfer caculations and in the	   *
*		synthesis.						   *
*									   *
*	NOTE:	This program uses a TI TMS320C30 EVM card to output the	   *
*		synthesized vowel signal.  Before starting this program	   *
*		an EVM routine, "io.out" must be loaded by executing	   *
*		"evmload" as						   *
*			>evmload -p 320 io.out				   *
*		The port address, 0x0320 can be specified by DOS setup	   *
*		as "SET D_OPTIONS= -p 320.  In this case, it is loaded by  *
*			>evmload io.out					   *
*               The following DOS commands should be included in	   *
*		"autoexec.bat" in order to set EVM environment;		   *
*									   *
*			SET D_DIR=C:\C30\SHELL				   *
*			SET D_SRC=d:\io\evm;d:\io\pc			   *
*			SET D_OPTIONS= -P 320				   *
*			SET C_DIR=C:\C30\TOOLS				   *
*			d:\io\evm\evmreset				   *
*									   *
*		In order to change the address, the following 2 things	   *
*		have to be done:					   *
*			1) Set appropriately the dip switch on the card.   *
*			   (See Page 3 of Technical Reference)		   *
*			2) "SET D_OPTIONS= -p xxx" in autoexec.bat	   *
*									   *
*									   *
*	NOTE :	In this version, the signal playback is suppressed.	   *
*		To reactivate playback, recover instruction lines around   *
*		variable "playback_flag" by the deletion of the comment	   *
*		markers.						   *
***************************************************************************/

#include	<always.h>

//#include	<grafix.h>
#include	<complex.h>
#include	<plot_lib.h>
#include	<menu.h>
#include	<lam_lib.h>
#include	<vtconfig.h>
#include	<vsyn_lib.h>
#include	<mydir.h>
#include	<vtcalcs.h>

/***********************( prototype definitions )**************************/

void	change_tract_config ( int left, int top );
void	change_physi_consts ( int left, int top );
void	change_save_options ( int left, int top );
void	uniform_tube ( int left, int top );
void	af_from_file( int left, int top, char *af_device );
void	T2model( int left, int top );
void	P3model( int left, int top );
void	ARTmodel ( int left, int top );
void	vowel_synt ( void );
void	save_par_sig( int model_type  );


/********************************( main )**********************************/

void main(void)
{
	static	char	*nt_filename = "nt.are",
			af_path[40];
	static	menu_state	ms;		/* menu state */
	static	menu_items	mm[4];		/* main menue */
	static	menu_state	ss;		/* sub_menu state */
	static	menu_items	sb[5];		/* sub_menu items */
	static	void	*imgbuf;
	int	i;

/* menu items */
	char *mainmenu[4] = {	"VT calculation",
				"Configuration",
				"Physical consts.",
				"Save options "		};
	char *modelmenu[5]= {	"uniform tube",
				"from file",
				"2 tube model",
				"3 parameter model",
				"articu. model"		};

/* read nasal tract area function from the file */

	strcpy(af_path, VTAFdir);
	strcat(af_path, nt_filename);
	read_af( af_path, &nna, &afnt );

/* Initialization */
	enter_graphics_mode( BGIdir );
	read_rad();

/* Initialize EVM */
/****
	if( init_evm() ) playback_flag = OFF;
	if( playback_flag )
	{  init_AIC( (int)smpfrq/1000, OFF, ON, LINE_V );
	   send_command_and_value(S, ON);
	}
****/
	cleardevice();

/* Menu specification */
	ms.n = 4;
	for(i=0; i<ms.n; i++) strcpy( mm[i].mssg, mainmenu[i]);

	ss.n = 5;
	for(i=0; i<ss.n; i++) strcpy( sb[i].mssg, modelmenu[i]);

/*** Main loop ***/

	for(;;)
	{  menu_main( &ms, mm );

	   if( ms.p == -1 ) break;		/* Esc hit, quit program */

	   switch( ms.p )
	   {
	   case 1 :				/* change loss conditions */
	      change_tract_config( mm[ms.p].left, mm[ms.p].bottom );
	      break;

	   case 2 :		/* change values of physical constants */
	      change_physi_consts( mm[ms.p].left, mm[ms.p].bottom );
	      break;

	   case 3 :				/* change save options */
	      change_save_options( mm[ms.p].left, mm[ms.p].bottom );
	      break;

	   case 0 :				/* vocal tract calculation */
	   {
	      ms.p = ms.p_old;			/* save old main menu pointer */
	      for(;;)
	      {  menu_horiz( 0, 0, &ss, sb, &imgbuf );

		 if( ss.p == -1 )	/* Esc hit */
		 {  ss.p = ss.p_old;	/* save old submenu pointer */
		    break;
		 }

		 switch( ss.p )
		 {
		 case 0	:			/* uniform tube */
		    uniform_tube(lotopix_x(55.), lotopix_y(10.));
		    break;
		 case 1 :			/* file specified area */
		    af_from_file(lotopix_x(55.), lotopix_y(10.), VTAFdir);
		    break;
		 case 2 :			/* 2-tube model */
		    T2model(lotopix_x(55.), lotopix_y(10.));
		    break;
		 case 3 :			/* 3-parameter model */
		    P3model(lotopix_x(55.), lotopix_y(10.));
		    break;
		 case 4	:			/* a-model specified area */
		    ARTmodel(lotopix_x(80.), lotopix_y(42.));
		    break;
		 default : break;
		 }
	      }					/* end of submenus */
	      break;
	   }					/* end of main case : 0 */
	   default : break;
	   }					/* end of main switch */
	}					/* end of main menus */

	quit();
}


/***************************************************************************
*                                                                          *
*                           local functions                                *
*                                                                          *
***************************************************************************/

/*****
*	Function : change_tract_config
*	Note :	modifies wall characteristics, yielding or rigid;
*		radiation boundary, Bessel function, RL-circuit, or
*		short circuit, at the lip or nostril opening; and
*		nasal tract, ON or OFF.
*****/

void	change_tract_config (
	int	left,		/* left coordinate value of menu box */
	int	top )		/* top                               */
{
	static	menu_state	s;
	static	menu_items	m[4];
	static	void		*imgbuf;
	int	i;

/* menu message specifications */
	char	*conf_menu[4] =	{
			"radiation load = %s",
			"wall           = %s",
			"nasal tract    = %s",
			"glottis        = %s"  };

	char	*conf_note[4] =	{
		"Enter 'r' or 'CR' to toggle between ON and OFF",
		"Enter 'w' or 'CR' to toggle between YIELDING and RIGID",
		"Enter 'n' or 'CR' to toggle between OFF and ON",
		"Don't touch !!"  };

/* initialize menu state and compose menu item massage */

	i = 0;

	if(rad_boundary==RL_CIRCUIT) sprintf(m[i].mssg, conf_menu[i],"ON");
	if(rad_boundary==SHORT_CIRCUIT)sprintf(m[i].mssg, conf_menu[i],"OFF");
	strcpy (m[i++].footnote, conf_note[i]);

	if( wall == YIELDING )sprintf(m[i].mssg, conf_menu[i], "YIELDING");
	else		      sprintf(m[i].mssg, conf_menu[i], "RIGID   ");
	strcpy(m[i++].footnote, conf_note[i]);

	if( nasal_tract == ON)sprintf(m[i].mssg, conf_menu[i], "ON ");
	else		      sprintf(m[i].mssg, conf_menu[i], "OFF");
	strcpy(m[i++].footnote, conf_note[i]);

	if(glt_boundary==CLOSE)sprintf(m[i].mssg, conf_menu[i], "CLOSE");
	else		    sprintf(m[i].mssg, conf_menu[i],    "OPEN ");
	strcpy(m[i++].footnote, conf_note[i]);

	s.n     = i;

/* prompt modifications by vertical menue */

	for(;;)
	{  menu_vert( left, top, &s, m, &imgbuf );
	   if( s.p == -1 )			/* Esc hit */
	   {  s.p = s.p_old;			/* save old pointer */
	      return;
	   }

	   switch( s.p )
	   {
	   case 0 :			/* toggle ON/OFF */
	      if(rad_boundary == RL_CIRCUIT)
	      {  rad_boundary = SHORT_CIRCUIT;
		 sprintf(m[s.p].mssg, conf_menu[s.p], "OFF");
	      }
	      else
	      {  rad_boundary = RL_CIRCUIT;
		 sprintf(m[s.p].mssg, conf_menu[s.p], "ON ");
	      }
	      break;

	   case 1 :			/* toggle wall = YIELDING/RIDGID */
	      if( wall == YIELDING )
	      {  wall = RIGID;
		 sprintf(m[s.p].mssg, conf_menu[s.p], "RIGID   ");
	      }
	      else
	      {  wall = YIELDING;
		 sprintf(m[s.p].mssg, conf_menu[s.p], "YIELDING");
	      }
	      break;

	   case 2 :			/* toggle nasal_tract = ON/OFF */
	      if( nasal_tract == ON )
	      {  nasal_tract = OFF;
		 sprintf(m[s.p].mssg, conf_menu[s.p], "OFF");
	      }
	      else
	      {  nasal_tract = ON;
		 sprintf(m[s.p].mssg, conf_menu[s.p], "ON ");
	      }
	      break;

	   case 3 :			/* toggle glt_boundary: CLOSE/OPEN */
	      if( glt_boundary == CLOSE )
	      {  glt_boundary = OPEN;
		 sprintf(m[s.p].mssg, conf_menu[s.p], "OPEN ");
	      }
	      else
	      {  glt_boundary = CLOSE;
		 sprintf(m[s.p].mssg, conf_menu[s.p], "CLOSE ");
	      }
	      break;

	   default : break;
	   }
	}
}

/*****
*	Function : change_physi_consts
*	Note :	modifies values of physical constants in simulation.
*		For example, "helium speech" can be simulated with
*		c = 84600 cm/s and rho = 0.00432 g/cm3, at t = 37 °C
*		and pressure 570 mBar.
*****/

void	change_physi_consts (
	int	left,		/* left coordinate value of menu box */
	int	top )		/* top                               */
{
	static	menu_state	s;
	static	menu_items	m[5];
	static	void		*imgbuf;
	int	i, max_len = 7;
	char	res [8];

	char	*physimenu[5] = {
			"air density (gm/cm3)  = %7.5f",
			"sound velocity (cm/s) = %7.1f",
			"wall resistance (gm/s/cm2) = %7.1f",
			"wall mass (gm/cm2)    = %7.1f",
			"wall compliance       = %7.1f"		};
	char	*physinote[5] = {
			"e.g., density of Helium = 0.00432 gm/cm3",
			"e.g., sound velocity in He = 84600 cm/s",
			"Mechanical resistance of vocal tract walls",
			"Wall mass per unit area",
			"Compliance of vocal tract walls"};

/* initialize menu state and compose menu item massage */

	i = 0;
	sprintf(m[i].mssg, physimenu[i], ro );
	strcpy (m[i++].footnote, physinote[i]);
	sprintf(m[i].mssg, physimenu[i], c );
	strcpy (m[i++].footnote, physinote[i]);
	sprintf(m[i].mssg, physimenu[i], wall_resi );
	strcpy (m[i++].footnote, physinote[i]);
	sprintf(m[i].mssg, physimenu[i], wall_mass );
	strcpy (m[i++].footnote, physinote[i]);
	sprintf(m[i].mssg, physimenu[i], wall_comp );
	strcpy (m[i++].footnote, physinote[i]);
	s.n = i;

/* prompt modifications by vertical menue */

	for(;;)
	{  menu_vert( left, top, &s, m, &imgbuf );

	   if( s.p == -1 )			/* Esc hit */
	   {  s.p = s.p_old;			/* save old pointer */
	      return;
	   }

	   get_str(m[s.p].pos, m[s.p].bottom, "new value = ", max_len, res);
	   switch( s.p )
	   {
	   case 0 :
	      sscanf (res, "%f", &ro);
	      sprintf(m[s.p].mssg, physimenu[s.p], ro );
	      break;
	   case 1 :
	      sscanf (res, "%f", &c);
	      sprintf(m[s.p].mssg, physimenu[s.p], c );
	      break;
	   case 2 :
	      sscanf (res, "%f", &wall_resi);
	      sprintf(m[s.p].mssg, physimenu[s.p], wall_resi );
	      break;
	   case 3 :
	      sscanf (res, "%f", &wall_mass);
	      sprintf(m[s.p].mssg, physimenu[s.p], wall_mass );
	      break;
	   case 4 :
	      sscanf (res, "%f", &wall_comp);
	      sprintf(m[s.p].mssg, physimenu[s.p], wall_comp );
	      break;
	   default : break;
	   }
	}
}


/******
*	Function : change_save_options
*	Note :	Filename consists of a fix code followed by an integer,
*		which is incremented every time save file(s) is (are)
*		created.  The code can be change in this menu.  With the
*		new code, the counter is zeroed.  The code length is less
*		or equal to 7.  If the code length is 5, then the maximum
*		1000 files can be created, i.e., from xxxxx0.DOC (SIG)
*		to xxxxx999.DOC, under the same code, "xxxxx".
*		In addition to the code, directories for save files and
*		for area_function file can be changed.
******/
void	change_save_options (
	int	left,		/* see af_uniform_tube */
	int	top )
{
	static	menu_state	s;
	static	menu_items	m[5];
	static	void		*imgbuf;
	int	i, max_len = 29;

	char	*savemenu[5] = {
			"Save file code = %s",
			"Save directory = %s",
			"temp_sig_path  = %s",
			"temp_glt_path  = %s",
			"Area directory = %s"	};

/* initialize menu state and compose menu item massage */

	i = 0;
	sprintf(m[i++].mssg, savemenu[i], save_file_code );
	sprintf(m[i++].mssg, savemenu[i], SIGdir );
	sprintf(m[i++].mssg, savemenu[i], tempSIGpath );
	sprintf(m[i++].mssg, savemenu[i], tempGLTpath );
	sprintf(m[i++].mssg, savemenu[i], VTAFdir );
	s.n = i;

/* prompt modifications by vertical menue */

	for(;;)
	{  menu_vert( left, top, &s, m, &imgbuf );

	   if( s.p == -1 )			/* Esc hit */
	   {  s.p = s.p_old;			/* save old pointer */
	      return;
	   }
	   switch( s.p )
	   {
	   case 0 :
	      get_str( m[s.p].pos, m[s.p].bottom, "new code = ",
		       max_len, save_file_code );
	      sprintf(m[s.p].mssg, savemenu[s.p], save_file_code );
	      save_file_count = 0;
	      break;

	   case 1 :
	      get_str( m[s.p].pos, m[s.p].bottom, "New directory = ",
		       max_len, SIGdir );
	      sprintf(m[s.p].mssg, savemenu[s.p], SIGdir );
	      break;

	   case 2 :
	      get_str( m[s.p].pos, m[s.p].bottom, "New path = ",
		       max_len, tempSIGpath );
	      sprintf(m[s.p].mssg, savemenu[s.p], tempSIGpath );
	      break;
	   case 3 :
	      get_str( m[s.p].pos, m[s.p].bottom, "New path = ",
		       max_len, tempGLTpath );
	      sprintf(m[s.p].mssg, savemenu[s.p], tempGLTpath );
	      break;
	   case 4 :
	      get_str( m[s.p].pos, m[s.p].bottom, "New directory = ",
		       max_len, VTAFdir );
	      sprintf(m[s.p].mssg, savemenu[s.p], VTAFdir );
	      break;

	   default : break;
	   }
	}
}

/******
*	Function : uniform_tube
*	Note :	A uniform tube is specified by two variables, A (cm**2)
*		and the length (cm) and the number of sections.
*****/
void	uniform_tube (
	int	left,		/* left-top coordinate of mssg. box */
	int	top  )          /* in pixels                        */
{
	static	menu_state	s;
	static	menu_items	m[7];
	static	void		*imgbuf;
	int	max_len = 5;
	char	res [6];
	int	nss_max = 100;		/* maximum number of sections */
	float	x;
	int	i;

	char	*UTmenu[7] = {
			"Calculate",
			"Synthesize",
			"Keep (save)",
			"Cross area (cm2) = %6.2f",
			"Tube length (cm) = %6.1f",
			"Number of sections = %4d",
			"Nasal couping (cm2)=%5.1f"	};


/* initialize menu state and compose menu item massage */

	nss = nss_uni;			/* initial default value */

	for(i=0; i<3; i++) strcpy ( m[i].mssg, UTmenu[i] );
	sprintf( m[i++].mssg, UTmenu[i], UTpar.area );
	sprintf( m[i++].mssg, UTmenu[i], UTpar.length );
	sprintf( m[i++].mssg, UTmenu[i], nss );
	if( nasal_tract == ON ) sprintf( m[i++].mssg, UTmenu[i], anc );
	s.n = i;

	i = 0;
	strcpy ( m[i++].footnote, note_calc );
	strcpy ( m[i++].footnote, note_synt );
	strcpy ( m[i++].footnote, note_save );

/* acclocate memory for area function */

	afvt = (area_function *) calloc( nss_max, sizeof(area_function) );

/* prompt modifications by vertical menue */
	for(;;)
	{
/* specify uniform area function */

	   x = UTpar.length/(float) nss;
	   nph = 9.0/x;		/* nasal brantch point is 9 cm above glottis */
	   nbu = nss - nph;
	   for(i=0; i<nss; i++)
	   {  afvt[i].A = UTpar.area;
	      afvt[i].x = x;
	   }
	   if( nasal_tract == ON )
	   {  anc = (float) min( anc, afvt[nph].A );
	      afvt[nph].A -= anc;
	      sprintf( m[6].mssg, UTmenu[6], anc );
	   }

	   plot_af ( nss, afvt, vp1 );

/* get menu command and process */

	   menu_vert( left, top, &s, m, &imgbuf );

	   if( s.p == -1 )			/* Esc hit */
	   {  s.p = s.p_old;			/* save old pointer */
	      sig_save_flag = OFF;
	      clear_vp( vp1.n );
	      clear_vp( vp2.n );
	      clear_vp( vp3.n );
	      free( afvt );
	      return;
	   }

/* modify uniform tube parameters */

	   if( s.p >= 3 )
	   {  sig_save_flag = OFF;
	      get_str(m[s.p].pos, m[s.p].bottom,"new value = ", max_len, res);
	   }

	   switch( s.p )
	   {
	   case 0 :
	      s.p = s.p_old;			/* save old pointer */

	      if(source_loc == 0)
		calplot_tf_FBA(UmUg_title, nfrmmax, frm, bw, amp, &nfrms, vp2);
	      else
		calplot_tf_FBA(UmPt_title, nfrmmax, frm, bw, amp, &nfrms, vp2);

	      plot_FBA( nfrms, frm, bw, amp, vp3 );
	      break;

	   case 1 :
	      s.p = s.p_old;
	      vowel_synt();
	      break;

	   case 2 :				/* save parameters */
	      s.p = s.p_old;
	      save_par_sig( UNIFORM_TUBE );
	      break;

	   case 3 :
	      sscanf (res, "%f", &UTpar.area);
	      sprintf( m[s.p].mssg, UTmenu[s.p], UTpar.area );
	      break;

	   case 4 :
	      sscanf (res, "%f", &UTpar.length);
	      sprintf(m[s.p].mssg, UTmenu[s.p], UTpar.length);
	      break;

	   case 5 :
	      sscanf (res, "%d", &nss_uni);
	      if( nss_uni > nss_max ) nss_uni = nss_max;
	      nss = nss_uni;
	      sprintf( m[s.p].mssg, UTmenu[s.p], nss );
	      break;

	   case 6 :
	      sscanf (res, "%f", &anc);
	      sprintf( m[s.p].mssg, UTmenu[s.p], anc );
	      break;

	   default : break;
	   }
	}
}


/******
*	Function : af_from_file
*	Note :	An area function is read from the file.
******/
void	af_from_file (
	int	left,		/* see af_uniform_tube */
	int	top,
	char	*af_dir  )	/* device and directories to .ARE files */
{
	char	af_path[40];
	FILE	*in;

	static	menu_state	s;
	static	menu_items	m[5];
	static	void		*imgbuf;
	int	max_len = 8;
	char	res [9];
	char	af_title[20];
	float	x;
	int	i, nss_max = 60;

	char	*FFmenu[5] = {
			"Calculate",
			"Synthesize",
			"Keep (save)",
			"File name = %s",
			"Nasal coupling(cm2) = %4.1f" };


/* initialize menu state and compose menu item massage */

	i = 0;
	strcpy (m[i++].mssg, FFmenu[i] );
	strcpy (m[i++].mssg, FFmenu[i] );
	strcpy (m[i++].mssg, FFmenu[i] );
	sprintf(m[i++].mssg, FFmenu[i], AFfilename );
	if( nasal_tract == ON ) sprintf(m[i++].mssg, FFmenu[i], anc);
	s.n = i;

	i = 0;
	strcpy (m[i++].footnote, note_calc );
	strcpy (m[i++].footnote, note_synt );
	strcpy (m[i++].footnote, note_save );
	strcpy (m[i++].footnote, note_vowel );

/* allocate memories for afvt */

	afvt = (area_function *) calloc(nss_max, sizeof(area_function) );

/* prompt modifications by vertical menue */

	for(;;)
	{
/* read area function file and plot it */

	   strcpy( af_path, af_dir );
	   strcat( af_path, AFfilename );
	   strcat( af_path, ".ARE" );
	   in = fopen( af_path, "rt" );
	   fscanf(in, "%s", af_title );
	   fscanf(in, "%d", &nss);
	   if( nss > nss_max ) nss = nss_max;
	   nph = 9;
	   if( nss < 11 ) nph = 5;
	   nbu = nss - nph;
	   fscanf(in, "%f", &x);
	   for(i=0; i<nss; i++)
	      { fscanf(in, "%f", &afvt[i].A); afvt[i].x = x; }
	   fclose ( in );

	   if( nasal_tract == ON )
	   {  anc = (float) min( anc, afvt[nph].A );
	      afvt[nph].A -= anc;
	      sprintf( m[4].mssg, FFmenu[4], anc );
	   }

	   plot_af ( nss, afvt, vp1 );

/* Get menu command and process */

	   menu_vert( left, top, &s, m, &imgbuf );

	   if( s.p == -1 )			/* Esc hit */
	   {  s.p = s.p_old;			/* save old pointer */
	      sig_save_flag = OFF;
	      clear_vp( vp1.n );
	      clear_vp( vp2.n );
	      clear_vp( vp3.n );
	      free( afvt );
	      return;
	   }

	   switch( s.p )
	   {
	   case 0 :				/* transfer calculation */
	      s.p = s.p_old;
	      if(source_loc == 0)
		calplot_tf_FBA(UmUg_title, nfrmmax, frm, bw, amp, &nfrms, vp2);
	      else
		calplot_tf_FBA(UmPt_title, nfrmmax, frm, bw, amp, &nfrms, vp2);
	      plot_FBA( nfrms, frm, bw, amp, vp3 );
	      break;

	   case 1 :				/* vowel synthesis */
	      s.p = s.p_old;
	      vowel_synt();
	      break;

	   case 2 :				/* save area_function etc */
	      s.p = s.p_old;
	      save_par_sig( FROM_FILE );
	      break;

	   case 3 :				/* get new area function */
	      sig_save_flag = OFF;
	      do
	      {  get_str( m[s.p].pos, m[s.p].bottom, "new name = ",
			  max_len, AFfilename );
		 strcpy( af_path,af_dir );
		 strcat( af_path, AFfilename );
		 strcat( af_path, ".ARE" );
	      }while( access( af_path, 0 ) == -1 );
	      sprintf(m[s.p].mssg, FFmenu[s.p], AFfilename );
	      break;

	   case 4 :				/* nasal coupling */
	      sig_save_flag = OFF;
	      get_str(m[s.p].pos, m[s.p].bottom,"new value = ", max_len, res);
	      sscanf (res, "%f", &anc);
	      sprintf( m[s.p].mssg, FFmenu[s.p], anc );
	      break;

	   default : break;
	   }
	}
}


/******
*	Function : T2model
*	Note :	An area function is specified by 2 tubes,
*
*		1) cross-sectional area of the 1-st section (A1) in cm2,
*		2) length of the 1-st section in cm,
*		3) area of the 2-nd section (A2) in cm2,
*		4) length of the 2-nd section in cm.
*
*****/
void	T2model (
	int	left,		/* left-top coordinate of mssg. box */
	int	top  )          /* in pixels                        */
{
	static	menu_state	s;
	static	menu_items	m[8];
	static	void		*imgbuf;
	int	max_len = 5;
	static	char	res [6];
	int	max_vtlen = 40;
	int	n1, n2;
	float	dx1, dx2, ratio;
	int	i;

	char	*T2menu[8] = {
				"Calculate",
				"Synthesize",
				"Keep (save)",
				"A1 area (cm2)  = %4.1f",
				"x1 length (cm) = %4.1f",
				"A2 area (cm2)  = %4.1f",
				"x2 length (cm) = %4.1f",
				"Nasal couping (cm2)= %5.2f"	};


/* initialize menu state and compose menu item massage */

	for(i=0; i<3; i++) strcpy ( m[i].mssg, T2menu[i] );
	sprintf( m[i++].mssg, T2menu[i], T2par.A1 );
	sprintf( m[i++].mssg, T2menu[i], T2par.x1 );
	sprintf( m[i++].mssg, T2menu[i], T2par.A2 );
	sprintf( m[i++].mssg, T2menu[i], T2par.x2 );
	if( nasal_tract == ON ) sprintf( m[i++].mssg, T2menu[i], anc );
	s.n = i;

	i = 0;
	strcpy ( m[i++].footnote, note_calc );
	strcpy ( m[i++].footnote, note_synt );
	strcpy ( m[i++].footnote, note_save );

/* allocate memory for area function */

	afvt = (area_function *) calloc( max_vtlen, sizeof(area_function) );

	for(;;)
	{
/* Compute an area function and plot it */

	   n1 = ceil( T2par.x1 );
	   n2 = ceil( T2par.x2 );
	   nss = n1 + n2;
	   if( nss > max_vtlen )
	   {  ratio = (float) nss / (float) max_vtlen;
	      n1 *= ratio;
	      n2 *= ratio;
	      nss = n1 + n2;
	   }
	   dx1 = T2par.x1/n1;
	   dx2 = T2par.x2/n2;
	   nph = 9;
	   nbu = nss - nph;

	   for(i=0; i<n1; i++)
	   {  afvt[i].x = dx1;
	      afvt[i].A = T2par.A1;
	   }
	   for(i=n1; i<nss; i++)
	   {  afvt[i].x = dx2;
	      afvt[i].A = T2par.A2;
	   }

	   if( nasal_tract == ON )
	   {  anc = (float) min( anc, afvt[nph].A );
	      afvt[nph].A -= anc;
	      sprintf( m[7].mssg, T2menu[7], anc );
	   }

	   plot_af ( nss, afvt, vp1 );

/* Get menu command and process */

	   menu_vert( left, top, &s, m, &imgbuf );

	   if( s.p == -1 )			/* Esc hit */
	   {  s.p = s.p_old;			/* save old pointer */
	      sig_save_flag = 0;
	      clear_vp( vp1.n );
	      clear_vp( vp2.n );
	      clear_vp( vp3.n );
	      free( afvt );
	      return;
	   }

/* modify uniform tube parameters */

	   if( s.p >= 3 )
	   {  sig_save_flag = OFF;
	      get_str(m[s.p].pos,m[s.p].bottom,"new value = ",max_len,res);
	   }

	   switch( s.p )
	   {
	   case 0 :
	      s.p = s.p_old;			/* save old pointer */
	      if(source_loc == 0)
		calplot_tf_FBA(UmUg_title, nfrmmax, frm, bw, amp, &nfrms, vp2);
	      else
		calplot_tf_FBA(UmPt_title, nfrmmax, frm, bw, amp, &nfrms, vp2);
	      plot_FBA( nfrms, frm, bw, amp, vp3 );
	      break;

	   case 1 :
	      s.p = s.p_old;
	      vowel_synt();
	      break;

	   case 2 :
	      s.p = s.p_old;
	      save_par_sig( T2MODEL );
	      break;

	   case 3 :
	      sscanf (res, "%f", &T2par.A1);
	      sprintf(m[s.p].mssg, T2menu[s.p], T2par.A1);
	      break;

	   case 4 :
	      sscanf (res, "%f", &T2par.x1);
	      sprintf(m[s.p].mssg, T2menu[s.p], T2par.x1);
	      break;

	   case 5 :
	      sscanf (res, "%f", &T2par.A2);
	      sprintf(m[s.p].mssg, T2menu[s.p], T2par.A2);
	      break;

	   case 6 :
	      sscanf (res, "%f", &T2par.x2);
	      sprintf(m[s.p].mssg, T2menu[s.p], T2par.x2);
	      break;

	   case 7 :
	      sscanf (res, "%f", &anc);
	      sprintf( m[s.p].mssg, T2menu[s.p], anc );
	      break;

	   default : break;
	   }
	}
}


/******
*	Function : P3model
*	Note :	An area function is specified by Fant's 4 tube model with
*		3 parameters,
*
*		1) location of the tongue-body center (Xt) in cm,
*		2) cross-sectional area (At) of the tounge-body section,
*		3) the mouth-opening area (Al) in cm2.
*
*		Note in the original model that the 3rd parameter is
*		the length over area ratio. But we assume, here, the lip
*		tube length is fixed to 1 cm.  Ref : Gunnar Fant (1960).
*		Acoustic Theory of Speech Production, Mouton, The Hague,
*		page 74.
*****/
void	P3model (
	int	left,		/* left-top coordinate of mssg. box */
	int	top  )          /* in pixels                        */
{
	static	menu_state	s;
	static	menu_items	m[7];
	static	void		*imgbuf;
	int	max_len = 6;
	static	char	res [7];
	float	area = 8.0;
	int	half_tng_len = 3, lip_len = 1, vt_len = 16;
	float	x;
	int	i;

	char	*P3menu[7] = {
			"Calculate",
			"Synthesize",
			"Keep (save)",
			"Cross area, At (cm2)   = %6.2f",
			"Position, Xt (cm)      = %6d",
			"Lip aperture, Al (cm2) = %6.2f",
			"Nasal couping (cm2)    = %6.2f"  };

	char	*P3note[3] = {
	   "Cross-sectional area of middle constriction (length = 6 cm)",
	   "Center position of middle constriction, measured from glottis",
	   "Lip aperture (length = 1 cm)"	};

/* initialize menu state and compose menu item massage */

	for(i=0; i<3; i++) strcpy ( m[i].mssg, P3menu[i] );
	sprintf( m[i++].mssg, P3menu[i], P3par.At);
	sprintf( m[i++].mssg, P3menu[i], P3par.Xt );
	sprintf( m[i++].mssg, P3menu[i], P3par.Al);
	if( nasal_tract == ON ) sprintf( m[i++].mssg, P3menu[i], anc );
	s.n = i;

	i = 0;
	strcpy ( m[i++].footnote, note_calc );
	strcpy ( m[i++].footnote, note_synt );
	strcpy ( m[i++].footnote, note_save );
	for( ; i<6; i++) strcpy ( m[i].footnote, P3note[i-3] );

/* prompt modifications by vertical menue */

	nss = vt_len + lip_len;
	nph = 9;
	nbu = nss - nph;
	afvt = (area_function *) calloc( nss, sizeof(area_function) );

	for(;;)
	{
/* Compute an area function and plot it */

	   for(i=0; i<vt_len+1; i++) afvt[i].x = 1.0;
	   for(i=0; i<P3par.Xt-half_tng_len; i++) afvt[i].A = area;
	   for(i=P3par.Xt-half_tng_len; i<P3par.Xt+half_tng_len; i++)
		  afvt[i].A = P3par.At;
	   for(i=P3par.Xt+half_tng_len; i<vt_len; i++) afvt[i].A = area;
	   afvt[vt_len].A = P3par.Al;
	   if( nasal_tract == ON )
	   {  anc = (float) min( anc, afvt[nph].A );
	      afvt[nph].A -= anc;
	      sprintf( m[6].mssg, P3menu[6], anc );
	   }

	   plot_af ( nss, afvt, vp1 );

/* Get menu command and process */

	   menu_vert( left, top, &s, m, &imgbuf );

	   if( s.p == -1 )			/* Esc hit */
	   {  s.p = s.p_old;			/* save old pointer */
	      sig_save_flag = 0;
	      clear_vp( vp1.n );
	      clear_vp( vp2.n );
	      clear_vp( vp3.n );
	      free( afvt );
	      return;
	   }

/* modify uniform tube parameters */

	   if( s.p >= 3 )
	   {  sig_save_flag = OFF;
	      get_str(m[s.p].pos,m[s.p].bottom,"new value = ",max_len,res);
	   }

	   switch( s.p )
	   {
	   case 0 :
	      s.p = s.p_old;			/* save old pointer */
	      if(source_loc == 0)
		calplot_tf_FBA(UmUg_title, nfrmmax, frm, bw, amp, &nfrms, vp2);
	      else
		calplot_tf_FBA(UmPt_title, nfrmmax, frm, bw, amp, &nfrms, vp2);
	      plot_FBA( nfrms, frm, bw, amp, vp3 );
	      break;

	   case 1 :
	      s.p = s.p_old;
	      vowel_synt();
	      break;

	   case 2 :
	      s.p = s.p_old;
	      save_par_sig( P3MODEL );
	      break;

	   case 3 :
	      sscanf (res, "%f", &P3par.At);
	      sprintf(m[s.p].mssg, P3menu[s.p], P3par.At);
	      break;

	   case 4 :
	      sscanf (res, "%d", &P3par.Xt);
	      P3par.Xt=min(vt_len - half_tng_len,max(P3par.Xt,half_tng_len));
	      sprintf(m[s.p].mssg, P3menu[s.p], P3par.Xt);
	      break;

	   case 5 :
	      sscanf (res, "%f", &P3par.Al);
	      sprintf(m[s.p].mssg, P3menu[s.p], P3par.Al);
	      break;

	   case 6 :
	      sscanf (res, "%f", &anc);
	      sprintf(m[s.p].mssg, P3menu[s.p], anc );
	      break;

	   default : break;
	   }
	}
}

/******
*	Function : ARTmodel
*	Note :	An area function is specified from the linear articulatory
*		model by specifing parameters
******/
void    ARTmodel (
	int	left,		/* menu position in pixels */
	int	top )
{
	static	menu_state	s;
	static	menu_items	m[12];
	static	void	*imgbuf;
		int	max_len = 6;
		char	res [7];

	static	char	*ARTmenu[12] = {
		"Calculate",
		"Syntesize",
		"Keep (save)",
		"Vowel = %s",
		"jaw   = %5.2f",
		"tongue= %5.2f",
		"shape = %5.2f",
		"apex  = %5.2f",
		"lip_ht= %5.2f",
		"lip_pr= %5.2f",
		"larynx= %5.2f",
		"nasal = %5.2f"	};

	static	char	*ARTnote[9] = {
		"Vowel = iy, ey, eh, ah, aa, ao, oh, uw, iw, ew, and oe",
		"Jaw position",
		"Tongue dorsum position",
		"Tongue dorsum shape",
		"Tongue apex position",
		"Lip height (aperture)",
		"Lip protrusion",
		"Larynx height",
		"Nasal coupling (cm2)"  };

	static	char	vowelcode[11][3]
		= { "iy", "ey", "eh", "ah", "aa", "ao", "oh", "uw",
		    "iw", "ew", "oe" };
	static	float	vowelpar[11][7]
		= {
		   {  0.5, -2.0, 1.0, -2.0,  1.0, -1.0, 0.0 },	/* iy */
		   {  0.0, -1.0, 1.0, -2.0,  1.0, -1.0, 0.0 },	/* ey */
		   { -1.0,  0.0, 1.0, -2.0,  1.0, -0.5, 0.0 },	/* eh */
		   { -1.5,  0.5, 0.0, -0.5,  0.5, -0.5, 0.0 },	/* ah */
		   { -1.5,  2.0, 0.0, -0.5,  0.5, -0.5, 0.0 },	/* aa */
		   { -0.4,  3.0, 1.5,  0.0, -0.3,  0.0, 0.0 },	/* ao */
		   {  -.7,  3.0, 1.5,  0.0, -0.6,  0.0, 0.0 },	/* oh */
		   {  0.5,  2.0, 1.5, -2.0, -1.0,  1.5, 0.0 },	/* uw */
		   {  0.5, -1.0, 1.0, -2.0, -0.5,  1.0, 0.0 },	/* iw */
		   {  0.0, -0.2, 1.0, -1.5, -0.25, 0.5, 0.0 },	/* ew */
		   { -1.0, -0.5, 0.5, -2.0,  0.2, -0.5, 0.0 }	/* oe */
		  };
		float	x;
		int	i, p;
		int	ns0 = 29;
	static	area_function	*af0;

/* Initialization */

	read_model_spec();

	af0 = (area_function *) calloc( ns0, sizeof(area_function) );
	nph = 9;
	nbu = 8;
	nss = nbu + nph;
	afvt  = (area_function *) calloc( nss, sizeof(area_function) );

	def_vp (vp4.n, vp4.x, vp4.y, vp4.x+vp4.w, vp4.y+vp4.h, 1);
	convert_scale( vp4.n );
	semi_polar();
/***	plot_semi_polar( vp4.n );  ***/

/* compose menu item massage */
	for(i=0; i<3; i++) 	strcpy (m[i].mssg, ARTmenu[i] );
				sprintf(m[i++].mssg, ARTmenu[i], "  " );
	for( ; i<11; i++)	sprintf(m[i].mssg, ARTmenu[i], AMpar[i-4] );
	if( nasal_tract == ON ) sprintf(m[i++].mssg, ARTmenu[i], anc );
	s.n = i;

	i = 0;
	strcpy(m[i++].footnote, note_calc );
	strcpy(m[i++].footnote, note_synt );
	strcpy(m[i++].footnote, note_save );
	for( ; i<12; i++) strcpy(m[i].footnote, ARTnote[i-3] );

/*** Loop ***/

	for(;;)
	{
/* Compute VT profile and area function, and plot them */

	   lam( AMpar );				/* profile	*/
	   clear_vp( vp4.n );
	   open_vp( vp4.n );
	   plot_vp(vp4.n);
	   plot_lam( vp4.n );

	   sagittal_to_area( &ns0, af0 );		/* area function */
	   appro_area_function( ns0, af0, nss, afvt);
	   if( nasal_tract == ON )
	   {  anc = (float) min( anc, afvt[nph].A );
	      afvt[nph].A -= anc;
	      sprintf(m[11].mssg, ARTmenu[11], anc );
	   }
	   plot_af ( nss, afvt, vp1 );

/* Get menu command and process */

	   menu_vert( left, top, &s, m, &imgbuf );

	   if( s.p == -1 )
	   {  s.p = s.p_old;		/* save old pointer and quit	*/
	      sig_save_flag = OFF;
	      clear_vp( vp1.n );
	      clear_vp( vp2.n );
	      clear_vp( vp3.n );
	      clear_vp( vp4.n );
	      free( af0 );
	      free( afvt );
	      return;
	   }

	   if( s.p >= 3 ) sig_save_flag = OFF;

	   switch( s.p )
	   {
	   case 0:			/* compute transfer function	*/
	      s.p = s.p_old;
	      if(source_loc == 0)
		calplot_tf_FBA(UmUg_title, nfrmmax, frm, bw, amp, &nfrms, vp2);
	      else
		calplot_tf_FBA(UmPt_title, nfrmmax, frm, bw, amp, &nfrms, vp2);
	      plot_FBA( nfrms, frm, bw, amp, vp3 );
	      break;

	   case 1 :
	      s.p = s.p_old;
	      vowel_synt();
	      break;

	   case 2:			/* save parameters etc.		*/
	      s.p = s.p_old;
	      save_par_sig( ART_MODEL );
	      break;

	   case 3:			/* default para. for selected vowel */
	      get_str(m[s.p].pos,m[s.p].bottom,"vowel code = ",max_len,res);
	      p = 0;
	      while( strcmp(res, vowelcode[p++]) && p<=11 );
	      if( p-- > 11 ) break;
	      sprintf(m[s.p].mssg, ARTmenu[s.p], res);
	      for(i=0; i<7; i++) AMpar[i] = vowelpar[p][i];
	      for(i=4; i<11; i++)
		 sprintf(m[i].mssg, ARTmenu[i], AMpar[i-4] );
	      break;

	   case 4: case 5: case 6: case 7: case 8: case 9: case 10:
	      sprintf(m[3].mssg, ARTmenu[3], "  ");
	      get_str(m[s.p].pos,m[s.p].bottom,"new value = ",max_len,res);
	      sscanf (res, "%f", &AMpar[s.p-4] );
	      sprintf(m[s.p].mssg, ARTmenu[s.p], AMpar[s.p-4]);
	      break;

	   case 11 :
	      get_str(m[s.p].pos, m[s.p].bottom,"new value = ",max_len, res);
	      sscanf( res, "%f", &anc );
	      sprintf(m[s.p].mssg, ARTmenu[s.p], anc );
	      break;

	   default : break;
	   }
	}
}

/*****
*	Function : vowel_synt
*	Note :	Stnthesize a stationary vowel with the current area
*		function.  The duration is about 300 ms, and a predefined
*		rise-fall F0 pattern is used.  The vowel signal and glottal
*		signal is stored in temporary files.  To transfer onto a
*		permanent file, select menu "save".
*****/
void	vowel_synt( void )
{
	int	buflen = 256;
	char	key;

	if( (sig_flux = fopen(tempSIGpath, "w+b")) == NULL )
	{  prompt(50.,30.,"Can't open temporary SIG file. Change path",ON);
	   return;
	}
	if( (glt_flux = fopen(tempGLTpath, "w+b")) == NULL )
	{  prompt(50.,30.,"Can't open temporary GLT file. Change path",ON);
	   return;
	}

	vowel_synthesis(sig_flux, glt_flux, buflen);

	fclose( sig_flux );
	fclose( glt_flux );
	sig_save_flag = ON;

/* Playback */
/****
	if( playback_flag )
	{
	   for(;;)
	   {  key = prompt(50.,30.,
	      "Hit 's' to play signal, 'g' glottal source, 'Esc' to return",
	      OFF);
	      if( key == ESC ) break;
	      if( key == 'g')  play_file( tempGLTpath );
	      else	       play_file( tempSIGpath );
	   }
	}
****/
}

/******
*	Function : save_par_sig
*	Note :	Save the current tract calculation conditions in a text
*		file.  The file name consists of the save_file_code
*		followed by the current value of save_file_count.
*		After the saving, the file counter is incremented by 1.
*		The catalogue code is .DOC for area function etc., and
*		.SIG for synthesized signal
******/
void	save_par_sig (
	int	model_type  )	/* the origin of area function :	*/
				/* UNIFORM_TUBE (=0), FROM_FILE (=1),	*/
				/* T2MODEL (=2), P3MODEL (=3),		*/
				/* and ART_MODEL (=4).			*/
{
	#define BLOCKSIZE 256
	char	save_path[80], number[4];
	FILE	*out, *in;
	int	*buf;
	int	i, np;

	strcpy(save_path, SIGdir);
	strcat(save_path, "\\");
	strcat(save_path, save_file_code);
	itoa(save_file_count, number, 10);
	strcat(save_path, number);
	strcat(save_path, ".DOC");

	if(access(save_path, 0) == 0)
	{ prompt( 50., 30., "File already exists ! Change code.", ON);
	  return;
	}

/* open file and save */
	out = fopen(save_path, "wt");

	fprintf(out, "%s\n", save_path);

/* Configuration options */
	fprintf(out, "\nVocal tract configuration\n");

	if(rad_boundary == SHORT_CIRCUIT)
	   fprintf(out, "\tRadiation load = OFF\n");
	if(rad_boundary == RL_CIRCUIT)
	   fprintf(out, "\tRadiation load = ON\n");

	if(wall == RIGID)
	   fprintf(out, "\tWall           = RIGID\n");
	if(wall == YIELDING)
	   fprintf(out, "\tWalls          = YIELDING\n");

	if(nasal_tract == OFF)
	   fprintf(out, "\tNasale tract   = OFF\n");
	if(nasal_tract == ON)
	   fprintf(out, "\tNasale tract   = ON\n");

	if(glt_boundary == CLOSE)
	   fprintf(out, "\tGlottis        = CLOSE\n");
	if(glt_boundary == OPEN)
	   fprintf(out, "\tGlottis        = OPEN\n");

/* Physical constants */
	fprintf(out, "\nPhysical constants\n");
	fprintf(out, "\tAir density     = %f\n", ro);
	fprintf(out, "\tSound velocity  = %f\n", c);
	fprintf(out, "\tWall mass       = %f\n", wall_mass);
	fprintf(out, "\tWall resistance = %f\n", wall_resi);
	fprintf(out, "\tWall compliance = %f\n", wall_comp);

/* Model type: the source of the area function */
	fprintf(out, "\nModel and parameters\n");
	switch(model_type)
	{
	case UNIFORM_TUBE :
	   fprintf(out, "\tmodel_type = UNIFORM_TUBE\n");
	   fprintf(out, "\tArea    = %f\n", UTpar.area);
	   fprintf(out, "\tLength  = %f\n", UTpar.length);
	   break;
	case FROM_FILE :
	   fprintf(out, "\tmodel_type    = FROM_FILE\n");
	   fprintf(out, "\tAFdevice   = %s\n", VTAFdir);
	   fprintf(out, "\tAFfilename = %s\n", AFfilename);
	   break;
	case T2MODEL:
	   fprintf(out, "\tmodel_type = T2MODEL\n");
	   fprintf(out, "\tA1 area (cm2)  = %f\n", T2par.A1);
	   fprintf(out, "\tx1 length (cm) = %d\n", T2par.x1);
	   fprintf(out, "\tA2 area (cm2)  = %f\n", T2par.A2);
	   fprintf(out, "\tx2 length (cm) = %d\n", T2par.x2);
	   break;
	case P3MODEL :
	   fprintf(out, "\tmodel_type = P3MODEL\n");
	   fprintf(out, "\tConstriction area     = %f\n", P3par.At);
	   fprintf(out, "\tConstriction Position = %d\n", P3par.Xt);
	   fprintf(out, "\tLip openning area     = %f\n", P3par.Al);
	   break;
	case ART_MODEL :
	   fprintf(out, "\tmodel_type = ART_MODEL\n");
	   fprintf(out, "\tArticulatory parameters\n");
	   for(i=0; i<7; i++) fprintf(out, "\t%f", AMpar[i]);
	   fprintf(out, "\n");
	   break;
	}

/* Area function */
	fprintf(out, "\nArea function: section area (cm2) and length (cm)\n");
	fprintf(out, "\tNumber of sections = %d\n", nss);
	for(i=0; i<nss; i++)
	   fprintf(out, "\t%3d\t%f\t%f\n", i, afvt[i].A, afvt[i].x );

/* Nasal coupling magnitude */
	if(nasal_tract == ON)
	   fprintf(out, "\nNasal coupling (cm2) = %f\n", anc);

/* Formant frequencies, bandwidth and amplitudes */
	fprintf(out, "\nFormant frequencies (Hz), bandwidth (Hz) and amplitudes (dB)\n");
	fprintf(out, "\tNumber of detected formants = %d\n", nfrms);
	for(i=0; i<nfrms; i++)
   fprintf(out, "\t%3d\t%7.1f\t%7.1f\t%7.1f\n", i+1, frm[i], bw[i], amp[i]);

/* Close .DOC file */
	fclose (out);

/* If exists, save the corresponding signal */

	if( sig_save_flag == ON )
/* open files */
	{  strcpy(save_path, SIGdir);
	   strcat(save_path, "\\");
	   strcat(save_path, save_file_code);
	   itoa(save_file_count, number, 10);
	   strcat(save_path, number);
	   strcat(save_path, ".SIG");
	   if( (out = fopen(save_path, "wb")) == NULL)
	   {  prompt(50.,30.,"Can't open SIG file. Do something", ON);
	      return;
	   }
	   if( (in  = fopen(tempSIGpath, "rb")) == NULL )
	   {  prompt(50.,30.,"Can't find temporary SIG file. Change path",ON);
	      return;
	   }
/* copy */
	   buf = (int *) calloc(BLOCKSIZE, sizeof(int));
	   do
	   {  np = fread (buf, sizeof(int), BLOCKSIZE, in);
		   fwrite(buf, sizeof(int), BLOCKSIZE, out);
	   }  while(!feof(in) && np == BLOCKSIZE);

/* clean up */
	   free(buf);
	   fclose(out);
	   fclose(in);
	   sig_save_flag = 0;
	}

/* incremant file counter and quit */
	save_file_count++;

}
