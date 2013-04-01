function pccfg = getpc(PCcfg)
%	GETPC Converts the data structure for physical constants
%		into a row vector

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/private/getpc.m 2     10/18/01 2:45p Satra $

% $NoKeywords: $

% Setup globals
global RELEASE

pccfg = [...
      PCcfg.ro PCcfg.c PCcfg.wall_resi ...
      PCcfg.wall_mass PCcfg.wall_comp];
