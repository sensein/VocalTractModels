/******
*	iowav.h
*	note:	definitions and type declarations for input/output
*		of .WAV files with MicroSoft multimedia format.
*****/

#pragma warning (disable : 4244 4101 4305)

typedef unsigned char       BYTE;
typedef unsigned short      WORD;
typedef unsigned long       DWORD;

#define HDReminderLEN_IN_BYTES  36		/* so-called 'chunk' */
#define WAVE_FORMAT_PCM     1			/* linear PCM codeing */
#define PACK4CHARS(ch0,ch1,ch2,ch3)    \
	((DWORD)(BYTE)(ch0)|((DWORD)(BYTE)(ch1)<<8) |    \
	((DWORD)(BYTE)(ch2)<<16) | ((DWORD)(BYTE)(ch3)<<24))

	/* waveform format structure (information common to all formats) */
typedef struct waveformat_tag {
    WORD    wFormatTag;        /* coding type */
    WORD    nChannels;         /* number of channels (i.e. mono, stereo, etc.) */
    DWORD   nSamplesPerSec;    /* sample rate */
    DWORD   nAvgBytesPerSec;   /* for buffer estimation */
    WORD    nBlockAlign;       /* block size of data */
} WAVEFORMAT;

	/* waveform format structure for PCM data */
typedef struct pcmwaveformat_tag {
    WAVEFORMAT  wf;
    WORD        wBitsPerSample;		/* 16 bits DA-AD */
} PCMWAVEFORMAT;

	/* file head structure */
typedef struct {
		DWORD Riff;		/* .WAV file code = RIFF */
		DWORD Ident1;		/* len (in bytes) of the reminder */
		DWORD Wave;		/* debut du chunk propre a WAVE */
		DWORD fmt;
		DWORD Ident2;		/* len (in bytes) of format */
		PCMWAVEFORMAT format;
		DWORD Data;
		DWORD Ident3;		/* len (in bytes) of signal */
		} HEADER;

void	sig_to_wav(char *sigpath, long smp_freq, int smp_bits, int nchannels,
	char *wavpath);
void	wav_to_sig(char *wavpath, char *sigpath);
FILE	*newWavFile(char *wavpath);
void	closeWavFile(FILE *wav, long smp_freq, int smp_bits, int nchannels,
	DWORD sampCount);
void	readWavFileHead(FILE *wav, long *smp_freq, int *smp_bits,
	int *nchannels, DWORD *smp_count);
