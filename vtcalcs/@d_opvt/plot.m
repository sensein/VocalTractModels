function plot(vtm,X,AV)
% PLOT  Plots the shape of the vocal tract
%   This function plots the shape of the vocal tract corresponding
%   to the provided articulatory positions.
%   X is the matrix of articulatory inputs
%   AV is a vector signifying voicing. Whenever AV is non-zero, a
%   signal is drawn on the body to signify that the current utterance 
%   is being voiced.

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/plot.m 5     10/24/01 11:40a Satra $

% $NoKeywords: $

% Setup globals
global RELEASE

% In order to speed up plotting this sin function has to be taken 
% out of the loop
fs = 10000;
t = [0:1/fs:0.01]';
y = sin(2*pi*t*500);

% scale and shift the plot to the appropriate location and size
AVval(:,1) = 275+50/0.01*t; 
AVval(:,2) = 525+25*y;

% Cycle through each input
for i=1:size(X,2),
    [Ft,Bt,At,Af{i},Tfm{i},Tff{i},P1{i},P2t] = doAM(X(:,i));

    % Draw the figure
    diva_drawman(P1{i},P2t,X(1,i),AV(i),AVval);
    drawnow;
end