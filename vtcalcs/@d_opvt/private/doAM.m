function [F,B,A,Af,Tf,P1,P2] = doAM(th)
% DOAM  A wrapper around AMgetdata to provide external parameters
%   This function serves to add parameters from the file vtsimul
%   and also converts the data.TC (tract constants) and data.PC 
%   (physical constants) values to an array form using gettc and 
%   getpc.

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/private/doam.m 3     10/18/01 2:45p Satra $

% $NoKeywords: $

% Setup globals
global RELEASE

%load miscellaneous data structures
vtsimul;

% reshape to ensure correct input form
th = th(:)';
data.AMpar.ampar = th;
[Af,Tf,F,B,A,P1,P2] = AMgetdata(...
    gettc(data.TC),...
    getpc(data.PC),...
    [data.AMpar.ampar,data.AMpar.anc]);