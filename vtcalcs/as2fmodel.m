function as2fmodel
%	AS2FMODEL Starts the Area functions model
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
      'Number of sliders' ...
      'Tube length (cm)' ...
      'Tube sections' ...
      'Nasal coupling (cm2)' ...
      'Close'};

set(obj1,'Userdata',data);
h = menuas2f(menustr,'as2fcb');
as2fcb(0);
