function changepc
%	CHANGEPC Changes the physical constant definitions
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
      'Air density (gm/cm3)' ...
      'Sound velocity (cm/s)' ...
      'wall resistance (gm/s/cm2)' ...
      'wall mass (gm/cm2)' ...
      'wall compliance' ...
      'Close'};

h = menuvt(menustr,'pccb',0);
pccb(0);
