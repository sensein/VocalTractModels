% GET_VSPE - function to generate vowel spectrum
% Usage: [amp,phase,freq]=get_vspe(f0,ft,bt,nft,rate,plotflag);
% Peter Assmann - Nov 21 1991
function [amp,phase,f]=get_vspe(f0,ft,bt,nft,rate,plotflag);
ft_neutral_vowel=[500 1500 2500 3500 4500];
bt_neutral_vowel=[100 100 100 100 100];
if(~exist('nft')) nft=5; end;
if(~exist('rate')) rate=10000; end;
if(~exist('f0')) f0=100; end;
if(~exist('ft')) ft=ft_neutral_vowel; end;
if(~exist('bt')) bt=bt_neutral_vowel; end;
if(~exist('plotflag')) plotflag=0; end;
nf=length(ft); nb=length(bt);
if(nf<nft) ft(nf+1:nft)=ft_neutral_vowel(nf+1:nft); end;
if(nb<nft) bt(nb+1:nft)=bt_neutral_vowel(nb+1:nft); end;
if length(f0)==1,  % f0: fundamental frequency
 f=0:f0:rate/2; 
else,  % f0: vector of frequencies
 f=f0; 
end;
tf=ones(size(f));
T=1/rate;
z=exp(i*2*pi*f*T);
% Combined glottal source/radiation function
% Glottal source spectrum has slope of -12 dB/oct
% Radiation function imparts a +6 dB/oct boost: 1-z.^(-1)
fgp=0; bgp=100;
cp=-exp(-2*pi*bgp*T);
bp=2*exp(-pi*bgp*T).*cos(2*pi*fgp*T);
ap=1-bp-cp;
tf=tf.*(ap./(1-bp.*z.^(-1)-cp.*z.^(-2))).*(1-z.^(-1));
% Glottal zero to shape the spectrum (Klatt, 1980: p.976)
fgz=1500; bgz=6000;
cz=-exp(-2*pi*bgz*T);
bz=2*exp(-pi*bgz*T).*cos(2*pi*fgz*T);
az=1-bz-cz;
cz=-cz/az;
bz=-bz/az;
az=1/az;
tf=tf.*(az./(1-bz.*z.^(-1)-cz.*z.^(-2)));
% vocal tract transfer function
% add contributions of individual formants
c=-exp(-2*pi*bt*T);
b=2*exp(-pi*bt*T).*cos(2*pi*ft*T);
a=1-b-c;
for ift=1:nft
 tf=tf.*(a(ift)./(1-b(ift).*z.^(-1)-c(ift).*z.^(-2)));
end
amp=abs(tf)';
phase=angle(tf)';
if plotflag==0 return; end;
clg,hold off;
%to avoid compilation errors
subplot(211),comb(f,20*log10(amp./max(amp)+.001)+60),
axis([0 rate/2 0 70]);
xlabel('Frequency (Hz)'); ylabel('Amplitude (dB)'); 
title('Amplitude spectrum'); grid;
subplot(212),comb(f,180/pi.*phase);
axis([0 rate/2 -180 180]);
xlabel('Frequency (Hz)'); ylabel('Phase (deg)'); 
title('Phase spectrum'); grid;
subplot(111),hold off;
