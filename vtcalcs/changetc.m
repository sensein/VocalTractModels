function changetc
%	CHANGETC Changes the tract definitions
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
      'Radiation load' ...
      'Wall' ...
      'Nasal tract' ...
      'Glottis' ...
      'Close'};

h = menuvt(menustr,'tccb',0);
tccb(0);
