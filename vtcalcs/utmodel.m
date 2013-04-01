function utmodel
%	UTMODEL Starts the Unifrom tube model
%		This function essenitally defines the menu items
%		and calls up the menu program.

% Copyright (c) 1999 Satrajit S. Ghosh (satra@bu.edu)
% Department of Cognitive and Neural Systems, Boston University

% $Revision: 1.00$ $Date:Fri Oct  1 17:02:05 EDT 1999$

% Bug fixes

% Modifications

obj1 = findobj('Tag','Vtcalcs');
data = get(obj1,'Userdata');
obj = findobj('Tag','AM');
axes(obj);
cla
set(obj,'visible', 'off');

menustr = {...
      'Calculate' 'Synthesize' 'Keep (save)' ...
      'Cross area (cm2)' ...
      'Tube length (cm)' ...
      'Number of sections' ...
      'Nasal coupling (cm2)' ...
      'Close'};

data.UTmpar.nss_max = 100;
set(obj1,'Userdata',data);
h = menuvt(menustr,'utcb');
utcb(0);
