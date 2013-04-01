function [F,B,A,Af,P1,P2] = convert1(vtm,X)
% CONVERT1  Converts articulatory parameters to acoustic measures
%   This function takes as input:
%       vtm - vocal tract object
%       X   - Matrix of articulatory sets, one per column
%   And provides as output:
%       F   - Formants
%       B   - Bandwidths
%       A   - Amplitude
%       Af  - Area functions
%       P1,P2 - Plotting parameters

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/convert1.m 4     10/18/01 2:45p Satra $

% $NoKeywords: $

% Setup globals
global RELEASE

% Initialize the outputs data to NaNs
F = NaN*zeros(5,1);
B = NaN*zeros(5,1);
A = NaN*zeros(5,1);

% Cycle through each input set
for i=1:size(X,2),
    [Ft,Bt,At,Afunc,Tfunc,P1{i},P2t] = doAM(X(:,i));
    Af{i} = Afunc;
    F(1:length(Ft),i)=Ft; 
    B(1:length(Bt),i)=Bt; 
    A(1:length(At),i)=At;
    P2(:,i) = P2t';
end
