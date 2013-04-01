function sig = diva_synth2(data)
% DIVA_SYNTH2   Synthesizes the sound
%   This core function relies on a mex file to provide the output signal
%   Input:
%   data: The structure data must have the following fields
%       Af,Ag0,AgP,F0, duration, dx

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/private/diva_synth2.m 7     11/02/01 10:33a Mshiffer $

% $NoKeywords: $

% Setup globals
global RELEASE

% Determine the dimensionality of the area functions
n = size(data.Af,2);

% Initialize all parameters for vtsynth - [note the transpose]
Ag0 = data.Ag0';
AgP = data.AgP';
F0  = data.F0';

TAF1 = [linspace(0,data.duration,n)' [0;2 * ones(n-1,1)] -1*ones(n,1)]';
TAF2 = data.Af;

% call vtsynth to generate the signal
if isfield(data,'fs'),
    [sig] = vtsynth(size(data.Af,1),Ag0,AgP,F0,TAF1,TAF2,data.dx,data.fs);
else,
    [sig] = vtsynth(size(data.Af,1),Ag0,AgP,F0,TAF1,TAF2,data.dx);
end;

% In order to reduce clicks we are multiplying the signal by a voicing
% component defined by the minimum of the area function at each point
%AV=interp1(linspace(1,length(sig),n),min(data.Af)>0.1,1:length(sig))';
%sig = sig.*AV;
