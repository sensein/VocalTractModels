function [sig] = vtsynth(nss,Ag0,AgP,F0,TAF1,TAF2,dx);
% VTSYNTH   Mex function which synthesizes the sounds
%   The parameters for this function are as follows:
%       - Number of cross sectional areas 
%       - Ag0 and AgP are vocal tract controlparameters.
%       - F0 is the formant contour
%       - TAF1 is the time points for the data
%       - TAF2 is the area function
%       - dx is the thickness of each area function
%   
%   The typical form of the parameters Ag0, AgP and F0 are of the form:
%   [timestamp1 value1 interpolation_scheme1;timestamp2 value2 interpolation_scheme2; ...]'
%   NOTE: Please note the transpose at the end of the above line -------------------------^
%       as opposed to diva_synth2 where the structure is without the transpose
%   where interpolation_scheme takes values: 
%           0 - set
%           1 - linear interpolation
%           2 - exponential interpolation
%   
%   See synth1 in the @d_opvt directory for an example of how to
%   use this function. Otherwise see the vt2k folder for the source
%   of the DLL and a demo routine

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/private/vtsynth.m 3     10/19/01 11:19a Satra $

% $NoKeywords: $

% Setup globals
global RELEASE
