/************************************************************************
*	File : vCv.c														*
*	Note : Synthesis with VT area-functions as the input. The temporal	*
*	patterns of vocal tract variables, area function and glottal area,	*
*	are specified by the following four input data files.               *
*																		*
*	Input files:														*
*		f_name1.AF	VT-AF specification (cm2, cm)						*
*		f_name2.AG0	glottal area baseline (cm2)							*
*		f_name2.AGP	glottal oscilation peak (cm2)						*
*		f_name2.F0	fundamental frequency (Hz)							*
*																		*
*	Patterns are specified by a series of triplets, {time (ms),         *
*   target_value (cm2 for area and Hz for F0), transition_type}. The    *
*   transtion type is either LIN or COS.  By default, the initial       *
*   triplet has time = 0 (ms) and type = ***. The last 3 files should   *
*   have the same name, since they describe together glottal source.    *
*   Moreover, the glottal source can be common to a different VCV       *
*   sequences. In this case, the file name could be like 'vsv'.         *
*                                                                       *
*	The tract shape at each sample (with the rate of 10 kHz) is			*
*	calculated by interpolations, and signal is synthesized using a		*
*	time-doamin acoustic-tube simulation. The signal will be written	*
*	with MicroSoft WAVE format.											*
*																		*
*	Output file:														*
*		f_name1.WAV, and f_name.SIG in option                           *
*                                                                       *
*   The file name is borrowed from the TP area file name.               *
*                                                                       *
************************************************************************/

//#include "stdafx.h"

#include	<stdio.h>
#include	<string.h>
#include <stdlib.h>
#include <conio.h>

#include	"vtconf.h"
#include	"iowav.h"
#include	"vcv_lib.h"

#define	NOISE_SOURCE_LEVEL 2.0e-8;
extern	char	*stopConsonants; //="dfghjkpqstvxz";

	int		sigFlag=1;			/* if = 1, *.WAV + *.SIG file is created */
	char	TPfileDirectory[_MAX_PATH];
	char	*rootDir = ""; //"D:\\Projects\\vt2k2\\";
	char	*WavDirectory = ""; //"D:\\Projects\\vt2k2\\vkv\\";
	char	*BGIdir = "C:\\TC\\";

int	AFSynthesis(char* GLTPpath, char* AFTPpath, char* wavPath, int sigFlag);

extern	double	Udc,		/* dc flow cm3/s */
		Ug,					/* glottal flow */
		Uac;				/* flow at supra-glottal constrction */

/******************************( main )**********************************/

static  AgMem = 0.;			/* memory for -6dB/Oct smoothing filter */

void main(void)
{
	char	AFTPfilename[10], AFTPpath[_MAX_PATH];
	char	GLTPfilename[10], GLTPpath[_MAX_PATH];
	char	CMDfilename[10], CMDpath[_MAX_PATH];
	char	subDirectory[10];
	char	wavPath[_MAX_PATH];
	FILE	*in;
	int		CMDflag, flag;
	int		i;

	float	noiseSourceLevel = 0.5; /* 1.5; */

/* Initialize graphics */
//	enter_graphics_mode(BGIdir);

/* Get file names */
	printf("Enter the subdirectory name; <ex., vkv> = ");
	scanf("%s", subDirectory);
	strcpy(TPfileDirectory, rootDir);
	strcat(TPfileDirectory, subDirectory);
	strcat(TPfileDirectory, "\\");

	printf("Enter glottal temporal-patern filename; <ex., vcv> = ");
	scanf("%s", GLTPfilename);

	nbp = 9;

	printf("Command file <y/n>?");
	puts("");
	if(_getch() == 'y')
	{ CMDflag = ON;
	  printf("Enter the filename <asa> = ");
	  scanf("%s", CMDfilename);
	  strcpy(CMDpath, TPfileDirectory);
	  strcat(CMDpath, CMDfilename);
	  strcat(CMDpath, ".CMD");
	  if((in = fopen(CMDpath, "rt")) == NULL)
	  { puts("No command file found");
	    _getch();
	    exit(-1);
	  }
	}
	else
	  CMDflag = OFF;

	for(;;)
	{ if(CMDflag == OFF)
	  { printf("Enter AF temporal-pattern fname <aakk05uw>, q to quit = ");
	    scanf("%s", AFTPfilename);
	    if(strlen(AFTPfilename) == 1)break;
	  }
	  else
	  { do
		  if((flag=fscanf(in, "%s", AFTPfilename)) == EOF)break;
		while(!strncmp(AFTPfilename, "//", 2));
		if(flag == EOF)break;
	  }

	/* filepath for *.WAV and input data files */

	  strcpy(wavPath, WavDirectory);
	  /**strcat(wavPath, subDirectory);
	  strcat(wavPath, "\\");**/
	  strcat(wavPath, AFTPfilename);	/* add .WAV, and .SIG later */

	  strcpy(GLTPpath, TPfileDirectory);
	  strcat(GLTPpath, GLTPfilename);

	  strcpy(AFTPpath, TPfileDirectory);
	  strcat(AFTPpath, AFTPfilename);

	/* simulation options */

	  vocal_tract = TIME_VARYING;
	  nasal_tract = OFF;
	  wall = YIELDING;
	  rad_boundary = RL_CIRCUIT;
	  wall_radiation = OFF;
	  noise_source = ON;		/* or OFF */
	  noiseAmp = noiseSourceLevel*NOISE_SOURCE_LEVEL;

	/* synthesis */
//	  cleardevice();

	  if(AFSynthesis(GLTPpath, AFTPpath, wavPath, sigFlag) == 0)
	  { puts("Can't find the path");
	    _getch();
	    exit(-1);
	  }
	  if(CMDflag == OFF)_getch();
//	  cleardevice();
	}								/** the end of for loop */
	if(CMDflag == ON)free(in);
	exit(0);
}


/****
*	function : AFSynthesis
*	note : VCV systhesis using the interpolation for the generation of
*		smoothly time-varying VT area function and glottal area wave.
*
*		The VT area function is specified by area-functions
*		and the sections, {A, x}. It can vary asynchronously
*		by increasing control points.
*
*		Hit any key to get out of this function after examining
*		the display.
****/
int	AFSynthesis(
	char*	GLTPpath,	/* glottal parameters	*/
	char*	AFTPpath,	/* file path for AF temp. pattern spec. */
	char*	wavPath,	/* path for signal file, .WAV will be added */
	int		sigFlag)	/* if = 1, *.SIG file is also created */
{
	DWORD	sampCount, xPoint;
	float	F0;
	int		T0;
	interpoMem	TAg0, TAgP, TF0;
	interpoMemVTAF	TAF;

	int		i, plotlen = 4500;	/* or 500 */
	float	Ag0, AgP;
	float	dt, tmin, tmax;

/*	vp_lo_frame	vp1 = {1, 5.0,  0., 90., 14.},
				vp2 = {2, 5.0, 15., 90., 14.},
				vp3 = {3, 5.0, 30., 90., 14.},
				vp4 = {4, 5.0, 45., 90., 14.},
				vp5 = {5, 5.0, 60., 90., 14.};
*/

	FILE	*wavOut, *sigOut;
	int		sig;
	float	fsig;
	/*int	nAcc;	(defined in vtconf.h) */
	/* note: consonantal constrition section address to be used
	   in the plotting function below for Ac. The addess is determined
	   as the minimum area by checking the area function of the
	   consonant.  The simulation program (vtt_sim) also detects the
	   minimum area to determine the source location, but this is done
	   for evry computatioan cycles.  The detected position thus can
	   vary along the tract during a VCV. nAcc is used to return Uac
	   at the fixed section.
	*/
	int		dummy;
	float	fDummy;
	char	AFfilename[15], sDummy[40], cCode[2];
	char    AFpath[_MAX_PATH], sigPath[_MAX_PATH];
	area_function	*af;
	int		nss;
	float	amin;

	strcpy(sigPath, wavPath);
	
	/* create .WAV file */
	//strcat(wavPath, ".WAV");
	//if((wavOut = newWavFile(wavPath)) == NULL) return(0);

	if(sigFlag==1)
	{
		strcat(sigPath, ".SIG");
		sigOut=fopen(sigPath, "wb");
	}

/* compose file names for temporal patterns */
	strcpy(TAg0.path, GLTPpath);
	strcat(TAg0.path, ".AG0");

	strcpy(TAgP.path, GLTPpath);
	strcat(TAgP.path, ".AGP");

	strcpy(TF0.path, GLTPpath);
	strcat(TF0.path, ".F0");

	strcpy(TAF.path, AFTPpath);
	strcat(TAF.path, ".AF");

	TAg0.endFlag = OFF;
	TAg0.cpCount = 0;
	TAgP.endFlag = OFF;
	TAgP.cpCount = 0;
	TF0.endFlag = OFF;
	TF0.cpCount = 0;
	TAF.endFlag = OFF;
	TAF.cpCount = 0;
	sampCount = 0;
	xPoint = 0;

/* detect constriction section address, nAcc, used for display */

	TAF.in = fopen(TAF.path, "rt");
	fscanf(TAF.in, "%d", &TAF.NControlPoints);
	fscanf(TAF.in, "%d", &dummy);

	for(nAcc=0, i=0; i<TAF.NControlPoints; i++)
	{ fscanf(TAF.in, "%f", &fDummy);		/* dummy read */
	  fscanf(TAF.in, "%s", AFfilename);
	  fscanf(TAF.in, "%s", sDummy);			/* dummy read */
	  cCode[0] = AFfilename[0];				/* copy the 1st letter */
	  cCode[1] = '\0';
	  strlwr(cCode);
	  if(strpbrk(cCode, stopConsonants))
	  { strcpy(AFpath, TPfileDirectory);
	    strcat(AFpath, AFfilename);
	    if(!readVTAFalloc(AFpath, sDummy, &nss, &af))
	    { printf("%s not found\n", AFpath);
	      _getch();
	      exit(-1);
	    }
	    for(amin=999., i=0; i<nss; i++)
	    if(amin >= af[i].A) {amin=af[i].A; nAcc=i;}
	    free(af);
	    break;
	  }
	}
	fclose(TAF.in);

/* synthesis and plotting loop */

	dt = 1000.0/smpfrq;			/* ms */
	do
	{
	/* time axis and graphics frame */
	  tmin = dt*sampCount;
	  tmax = tmin +  dt*plotlen;
//	  setcolor(WHITE);

	  //plotAreaFunctionFrame("AF(cm2)", 0., 20., vp1);

	  //plotAcAgFrame(" ", "Ag & Ac", tmin, tmax, vp2);
	  /*settextjustify(LEFT_TEXT, BOTTOM_TEXT);*/
	  /*plottextxy(datolo_x(vp2.n, 420.), datolo_y(vp2.n, 1.), "(a)");*/

	  /*plotF0Frame(" ", "F0(Hz)", tmin, tmax, vp1);*/
	  /*plotGlottalAreaFrame(" ", "AgP(cm2)", tmin, tmax, vp2);*/
	  /*plotGlottalAreaFrame(" ", "Ag0(cm2)", tmin, tmax, vp3);*/
	  //plotGlottalAreaFrame(" ", "Ag(cm2)", tmin, tmax, vp4);

	  //plotUdcFrame(" ", "Uc (cm3/s)", tmin, tmax, vp3);
	  /*settextjustify(LEFT_TEXT, BOTTOM_TEXT);*/
	  /*plottextxy(datolo_x(vp3.n, 420.), datolo_y(vp3.n, 500.), "(b)");*/

	  /*plotUdcFrame(" ", "Ug (cm3/s)", tmin, tmax, vp4);*/

	  //plotSpeechSignalFrame("time (ms)", "Signal", tmin, tmax, vp5);/* vp5 */
	  /*settextjustify(LEFT_TEXT, BOTTOM_TEXT);*/
	  /*plottextxy(datolo_x(vp4.n, 420.), datolo_y(vp4.n, 1.3), "(c)");*/

	  for(i=0; i<plotlen; i++)
	  {
	  /* compute glottal area by interpolations */
	    Ag0 = interpolX(sampCount, &TAg0);
	    AgP = interpolX(sampCount, &TAgP);
	    F0  = interpolX(sampCount, &TF0);
	    if(TF0.Ttype == RESET) xPoint = TF0.nextInterPoint;
	    if(xPoint == sampCount)
	    { T0 = 0;
	      if(F0 > 0.)T0 = (int)(smpfrq/F0 + .5);	/* refresh pitch */
	      xPoint += T0;
	    }
		Ag = FantGlottalArea(AgP, &T0) + Ag0;
		Ag = (Ag + 0.8*AgMem)/1.8; AgMem = Ag;		/* smoothing */

	  /* compute VT area-function by interpolation */
	    interpolVTAF(sampCount, &TAF);
	    if(sampCount == 0)vtt_ini();	/* initialize simulater */
	    fsig = vtt_sim();
		sig = (int)(DACscale*fsig+0.5);
		//fwrite(&sig, sizeof(int), 1, wavOut);

		if(sigFlag==1)fwrite(&fsig, sizeof(float), 1, sigOut);

	  /* plot data */
		//plotAF(vp1.n, nvt, afvt);
		//plotAcAgData(vp2.n, dt*sampCount, afvt[nAcc].A, Ag);

		if(TF0.Ttype != RESET)
		{ /*plotData(vp1.n, dt*sampCount, F0);*/
		  /*plotData(vp2.n, dt*sampCount, AgP);*/
		}
		/*plotData(vp3.n, dt*sampCount, Ag0);*/
		//plotData(vp4.n, dt*sampCount, Ag);

		/*plotData(vp3.n, dt*sampCount, (float)Uac);*/
		//plotData(vp3.n, dt*sampCount, (float)Udc);
		/*plotData(vp4.n, dt*sampCount, (float)Ug);*/
		//plotData(vp5.n, dt*sampCount, fsig);     /* vp4.n */

	    sampCount++;
	    if(TAg0.endFlag == ON && TAgP.endFlag == ON)break;
	   }
	}while(TAg0.endFlag == OFF || TAgP.endFlag == OFF);

/* termination procedurs */
	vtt_term();
	free(afvt);
	//closeWavFile(wavOut, (long)smpfrq, 16, 1, sampCount-1);
	if(sigFlag==1)fclose(sigOut);
	return(1);
}
