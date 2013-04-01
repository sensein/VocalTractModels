% track_synth: synthesizes speech from a set of 
% formant tracks - frequency-domain version of
% Klatt cascade synthesizer
% Peter Assmann - Mar 24 1993
% Usage: y=track_synth(dur,nmspf,rate,F1,F2,F3,F4,F5,AV,F0,sc,B1,B2,B3,B4,B5);
%  dur: duration of signal in ms
%  nmspf: number of ms per frame
%  rate: sample rate in Hz
% [parameter vectors, each of dimension 1 x nframes]
%  nframes=round(dur/nmspf)]:
%  F1: First formant frequency in Hz
%  F2: Second formant frequency in Hz
%  F3: Third formant frequency in Hz
%  F4: Fourth formant frequency in Hz
%  F5: Fifth formant frequency in Hz
%  B1-B5: formant bandwidths in Hz
%  AV: amplitude of voicing in dB
%  F0: Fundamental frequency in Hz
%  sc:  scale factor (default 120) 
% Peter Assmann - Sept 26, 1994
function y=track_synth(dur,nmspf,rate,F1,F2,F3,F4,F5,AV,F0,sc,B1,B2,B3,B4,B5);
if ~exist('dur'), dur=200; end;
if ~exist('nmspf'), nmspf=5; end;
if ~exist('rate'), rate=8000; end;
if ~exist('sc'), sc=120; end;
npts=dur*(rate/1000);
nsampf=nmspf*(rate/1000);
nframes=round(dur/nmspf);
% synthesize waveform
T=1/rate;
o=ones(nsampf,1);
if ~exist('F1') F1=[o*500]; end;
if ~exist('F2') F2=[o*1500]; end;
if ~exist('F3') F3=[o*2500]; end;
if ~exist('F4') F4=[o*3500]; end;
if ~exist('F5') F5=[o*4500]; end;
if ~exist('AV') AV=[o*60]; end;
if ~exist('F0') F0=[o*100]; end;
if ~exist('B1'), B1=[o*90]; end;
if ~exist('B2'), B2=[o*110]; end;
if ~exist('B3'), B3=[o*170]; end;
if ~exist('B4'), B4=[o*250]; end;
if ~exist('B5'), B5=[o*300]; end;
%
if length(F1)~=nframes, F1=splint(F1,nframes);end;
if length(F2)~=nframes, F2=splint(F2,nframes);end;
if length(F3)~=nframes, F3=splint(F3,nframes);end;
if length(F4)~=nframes, F4=splint(F4,nframes);end;
if length(F5)~=nframes, F5=splint(F5,nframes);end;
if length(AV)~=nframes, AV=splint(AV,nframes);end;
if length(F0)~=nframes, F0=splint(F0,nframes);end;
if length(B1)~=nframes, B1=splint(B1,nframes);end;
if length(B2)~=nframes, B2=splint(B2,nframes);end;
if length(B3)~=nframes, B3=splint(B3,nframes);end;
if length(B4)~=nframes, B4=splint(B4,nframes);end;
if length(B5)~=nframes, B5=splint(B5,nframes);end;
y=zeros(npts,1);
ps=zeros(nsampf,(rate/2)/min(F0)+1);
psi=zeros(1,(rate/2)/min(F0)+1);
for iframe=1:nframes,
 ft=[F1(iframe) F2(iframe) F3(iframe) F4(iframe) F5(iframe)];
 bt=[B1(iframe) B2(iframe) B3(iframe) B4(iframe) B5(iframe)];
 f0=F0(iframe);
% set number of formants based on sample rate
% exclude formants with frequencies <30 Hz
% (Note: nft must be set within frame loop!)
 if(rate>=8000) nft=4; end;
 if(rate>=10000) nft=5; end;
 i=find(ft(1:nft)<30); 
 n1=(iframe-1)*nsampf; n2=n1+nsampf-1; % frame sample indices
 nft=nft-length(i); bt(i)=[]; ft(i)=[];
 if nft>0, 
  [amp,ph]=get_vspe(f0,ft,bt,nft,rate,0); 
  nh=length(amp);  % number of harmonics to be synthesized
  k=(o*(0:nh-1))*2*pi*f0*T;
  ps=cumsum(k);
  y(n1+1:n2+1)=sum((o*amp'.*cos(ps-k+o*ph'+o*psi(1:nh)))')'...
  *10.^(AV(iframe)/20)*sc;
 else,  % no formants > 30 Hz
  y(n1+1:n2+1)=zeros(size(n1+1:n2+1)); % generate silence
 end;
 psi(1,1:nh)=psi(1,1:nh)+(ps(nsampf,1:nh)); 
end;
% if last frame exceeds npts, trim off any extra samples:
if n2+1>npts, y(npts+1:n2+1)=[]; end;
