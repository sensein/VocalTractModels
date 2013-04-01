function [F]=frmnts1(a,fs)
% FRMNTS1   Determine the formants given the LPC coefficients
%   The function determines the formans by taking the roots of the
%   LPC polynomial. For further details, see digital speech processing
%   by Rabiner and Shafer
%   Input: a    : LPC coefficients
%          fs   : Sampling frequency

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/private/frmnts1.m 5     10/19/01 11:19a Satra $

% $NoKeywords: $

% Setup globals
global RELEASE

const=fs/(2*pi);
rts=roots(a);
k=1;
save = [];

for i=1:length(a)-1
    re=real(rts(i)); 
    im=imag(rts(i));
    
    formn=const*atan2(im,re);        %formant frequencies
    bw=-2*const*log(abs(rts(i)));    %formant bandwidth
    
    % If the bandwidths and formants are reasonable, save them   
    if formn>90 & bw <700 & formn<4000
        save(k)=formn;
        bandw(k)=bw;
        k=k+1;
    end
    
end

[y, ind]=sort(save);
F = [NaN NaN NaN];
F(1:min(3,length(y))) = y(1:min(3,length(y)));