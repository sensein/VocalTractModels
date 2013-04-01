/*****
*	file: iowav_lb.c
*	note :	functions related sigal input/output using MicoSoft
*		multimedia tools
******/
//#include "stdafx.h"

#include <stdlib.h>
#include <io.h>
#include <conio.h>
#include <stdio.h>
#include <string.h>
#include "iowav.h"

typedef union{DWORD val; char text[5];} ID;

void sig_to_wav(
	char	*sigpath,		/* path to .SIG file */
	long	samp_freq,		/* in Hz */
	int	samp_bits,		/* 16 bits */
	int	nchannels,		/* =1 for mono, =2 for stereo */
	char	*wavpath)		/* path to .WAV file */
{

	FILE	*fsig, *fwav;
	BYTE	buf[1];
	int	block_align;
	long	siglen_in_bytes;
	HEADER	myhead[1];

	if ((fsig=fopen(sigpath,"rb"))==NULL)
	{ puts("can't find .SIG file");
	  exit(-1);
	}
	if ((fwav=fopen(wavpath,"wb"))==NULL)
	{ puts("can't create .WAV file");
	  exit(-1);
	}

	siglen_in_bytes=_filelength(_fileno(fsig));

	block_align = (WORD)samp_bits*nchannels/8;
	myhead->Riff   = PACK4CHARS('R','I','F','F');
	myhead->Ident1 = siglen_in_bytes + HDReminderLEN_IN_BYTES;
	myhead->Wave   = PACK4CHARS('W','A','V','E');
	myhead->fmt    = PACK4CHARS('f','m','t',' ');
	myhead->Ident2 = sizeof(PCMWAVEFORMAT);

	myhead->format.wf.wFormatTag      = WAVE_FORMAT_PCM;
	myhead->format.wf.nChannels       = nchannels;
	myhead->format.wf.nSamplesPerSec  = samp_freq;
	myhead->format.wf.nAvgBytesPerSec = samp_freq*block_align;
	myhead->format.wf.nBlockAlign     = block_align;
	myhead->format.wBitsPerSample     = samp_bits;

	myhead->Data   = PACK4CHARS('d','a','t','a');
	myhead->Ident3 = siglen_in_bytes;

	fwrite(myhead,sizeof(HEADER),1,fwav);

	while (!feof(fsig))
	{ fread(buf,sizeof(BYTE),1,fsig);
	  fwrite(buf,sizeof(BYTE),1,fwav);
	}
	fclose(fsig);
	fclose(fwav);
}

void wav_to_sig(
	char	*wavpath,		/* path to .SIG file */
	char	*sigpath)		/* path to .WAV file */
{

	FILE	*fsig, *fwav;
	long	smp_freq;		/* in Hz */
	int	smp_bits;		/* 16 bits */
	int	nchannels;		/* =1 for mono, =2 for stereo */
	DWORD	i, smp_count;		/* number of samples */

	BYTE	buf[1];
	int	block_align;
	long	siglen_in_bytes;
	HEADER	myhead[1];

	if ((fwav=fopen(wavpath,"rb"))==NULL)
	{ puts("can't find .WAV file");
	  exit(-1);
	}
	if ((fsig=fopen(sigpath,"wb"))==NULL)
	{ puts("can't create .SIG file");
	  exit(-1);
	}

/* read HEADER */

	readWavFileHead(fwav, &smp_freq, &smp_bits, &nchannels, &smp_count);		/* number of samples */

/* copy signal on .WAV to .SIG without head */
	/*while (!feof(fwav))*/
	for(i=0; i<2*smp_count; i++)
	{ fread(buf,sizeof(BYTE),1,fwav);
	  fwrite(buf,sizeof(BYTE),1,fsig);
	}
	fclose(fwav);
	fclose(fsig);

}

static	HEADER	wavHead[1];

/*****
*	function: newWavFile
*	note:	Create .WAV file and reserve space for the head. If fail,
*		the argument "wav" points NULL.
*
*		In the calling program, signal samples can be written by :
*
*		FILE	*wav;
*		int	buf;
*		DWORD	smp_count=0;
*		   .......
*		wav = newWavFile(wavPath);
*		   .......
*		buf[0] = sig;
*		fwrite(&buf,sizeof(int),1,wav);
*		smp_cout++;
*		   .......
*		closeWavFile(wav, smp_freq, smp_bits, nchannels, smp_count);
*		   .......
*
*****/
FILE *newWavFile(			/* return file flux */
	char	*wavPath)		/* path to .WAV file */
{
	FILE	*wav;

/* create a .WAV file */
	if ((wav=fopen(wavPath,"wb"))==NULL)return(NULL);

/* dummy write of the head */
	fwrite(wavHead,sizeof(HEADER),1,wav);
	return(wav);
}

/*****
*	function: closeWavFile
*	note:	write head and close .WAV file.
*****/

void closeWavFile(
	FILE	*wav,
	long	samp_freq,		/* in Hz */
	int	samp_bits,		/* 16 bits */
	int	nchannels,		/* =1 for mono, =2 for stereo */
	DWORD	smp_count)		/* number of samples in WORDs */
{

	WORD	block_align;

	wavHead->Riff   = PACK4CHARS('R','I','F','F');
	wavHead->Ident1 = 2*smp_count + HDReminderLEN_IN_BYTES;
	wavHead->Wave   = PACK4CHARS('W','A','V','E');
	wavHead->fmt    = PACK4CHARS('f','m','t',' ');
	wavHead->Ident2 = sizeof(PCMWAVEFORMAT);

	block_align = samp_bits*nchannels/8;
	wavHead->format.wf.wFormatTag      = WAVE_FORMAT_PCM;
	wavHead->format.wf.nChannels       = nchannels;
	wavHead->format.wf.nSamplesPerSec  = samp_freq;
	wavHead->format.wf.nAvgBytesPerSec = samp_freq*block_align;
	wavHead->format.wf.nBlockAlign     = block_align;
	wavHead->format.wBitsPerSample     = samp_bits;

	wavHead->Data   = PACK4CHARS('d','a','t','a');
	wavHead->Ident3 = 2*smp_count;			/* in bytes */

	rewind(wav);
	fwrite(wavHead,sizeof(HEADER),1,wav);
	fclose(wav);
}


/*****
*	function: readWavFileHead
*	note:	read the header of .WAV file, and return samping frequency,
*		etc. to the calling program.  Also checks some of the
*		characteristics. After return to the calling program,
*		the signal can be read by "fread".
*****/

void readWavFileHead(
	FILE	*wav,
	long	*smp_freq,		/* in Hz */
	int	*smp_bits,		/* 16 bits */
	int	*nchannels,		/* =1 for mono, =2 for stereo */
	DWORD	*smp_count)		/* number of samples */
{
	ID	fileType, wavFormat;

	fread(wavHead,sizeof(HEADER),1,wav);

/* check if RIFF */
	fileType.val = wavHead->Riff;
	strcpy(fileType.text+4, "\0");
	if(strcmp(fileType.text, "RIFF") != 0)
	{ printf("FileType = %s, (!=RIFF)\n", fileType.text);
	  _getch();
	  exit(-1);
	}
/* check if WAVE format */
	wavFormat.val = wavHead->Wave;
	strcpy(wavFormat.text+4, "\0");
	if(strcmp(wavFormat.text, "WAVE") != 0)
	{ printf("Format = %s, (!= WAVE)\n", wavFormat.text);
	  _getch();
	  exit(-1);
	}

/* check if linear PCM */
	if(wavHead->format.wf.wFormatTag != WAVE_FORMAT_PCM)
	{ puts("Not a linear PCM");
	  _getch();
	  exit(-1);
	}

	*nchannels = wavHead->format.wf.nChannels;
	*smp_freq  = wavHead->format.wf.nSamplesPerSec;
	*smp_bits  = wavHead->format.wBitsPerSample;
	*smp_count = wavHead->Ident3/2;			/* # of samples */
}