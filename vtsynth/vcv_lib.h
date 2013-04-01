/*****
*	File :	vcv_lib.h
*	Note :	Prototypes etc..
*****/
//#include <dir.h>
#pragma warning (disable : 4244 4101 4305 4018)

#define NSSMAX	40
#define	DACscale 28000.0		/* scale-up signal to fit 16 bits integer */
#define gltDACscale 100.*DACscale	   /* sacle up for glottal source */
#define	DWORD unsigned int

#ifndef _MAX_PATH
#define _MAX_PATH 256
#endif

#if !defined(AREA_FUNCTION)
  #define AREA_FUNCTION area_function;
  typedef struct{float A, x;}AREA_FUNCTION;
#endif

typedef	enum{RESET, LINEAR, COSINE, EXPO}transiType;
	/**
	RESET  : reset interpolation (pad zeros the segment);	"SET"
	LINEAR : linear transitions;							"LIN"
	COSINE : cosine transitions;							"COS"
	EXPO   : exponential (1st order) transition;			"EXP"
	**/

typedef	struct{
	char	path[_MAX_PATH];
	FILE*	in;
	int	endFlag;
	int	cpCount, NControlPoints;
	transiType Ttype;
	DWORD	interPoint, nextInterPoint;
	float *time;
	float *value;
	int   *type;
	float	a, b, w;
	}interpoMem;

typedef	struct{
	char	path[_MAX_PATH];
	FILE*	in;
	int	endFlag;
	int	cpCount, NControlPoints;
	transiType Ttype;
	DWORD	interPoint, nextInterPoint;
	float	aA, bA, w, ax, bx;
	}interpoMemTube;

typedef	struct{
	char	path[_MAX_PATH];
	FILE*	in;
	int	endFlag;
	int	cpCount, NControlPoints;
	transiType Ttype;
	DWORD	interPoint, nextInterPoint;
	float	aA[NSSMAX], bA[NSSMAX], w, ax[NSSMAX], bx[NSSMAX];
	int N;
	float dx;
	float **af;
	float *time;
	int *nloc;
	int *type;
	}interpoMemVTAF;

/* Plototypes */
int	readVTAFalloc(char* AFpath, char* AFtitle, int* nss,
	area_function **af);
void	readVTAF(char* AFpath, char* AFtitle, int nvt, area_function *af);
float	FantGlottalArea(float Ap, int* T0);
float	MaedaGlottalArea(float Ap, int* T0);

float	interpolX(DWORD sampCount, interpoMem* m);
area_function	interpolAx(DWORD sampCount, interpoMemTube* m);
void	interpolVTAF(DWORD sampCount, interpoMemVTAF* m);
void	FourTubeAF(area_function T1,area_function T2,
		   area_function T3, area_function T4);

/*
void	plotData(int vpn, float t, float y);
void	plotAcAgData(int vpn, float t, float Ac, float Ag);
void	plotAF(int vpn, int ns, area_function *af);

void	plotGlottalAreaFrame(char* xtitle, char* ytitle, float, float,
		vp_lo_frame);
void	plotF0Frame(char* xtitle, char* ytitle, float, float, vp_lo_frame);
void	plotUdcFrame(char* xtitle, char* ytitle, float, float, vp_lo_frame);
void	plotAcFrame(char* xtitle, char* ytitle, float tmin, float tmax,
		vp_lo_frame vp);	// frame in logical units 
void	plotAcAgFrame(char* xtitle, char* ytitle, float tmin, float tmax,
		vp_lo_frame vp);	// frame in logical units 
void	plotSpeechSignalFrame(char* xtitle, char* ytitle, float, float,
		vp_lo_frame);
void	plotAreaFunctionFrame(char* ytitle, float, float, vp_lo_frame);
*/

int	vtt_ini(void);
float	vtt_sim(void);
void	vtt_term(void);
