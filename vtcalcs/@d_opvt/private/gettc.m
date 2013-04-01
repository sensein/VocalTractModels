function tccfg = gettc(TCcfg)
%	GETPC Converts the data structure for tract configuration
%		into a row vector

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/private/gettc.m 2     10/18/01 2:45p Satra $

% $NoKeywords: $

% Setup globals
global RELEASE

tccfg = [...
      TCcfg.rad_boundary TCcfg.wall...
      TCcfg.nasal_tract TCcfg.glt_boundary];
