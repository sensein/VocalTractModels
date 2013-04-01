function h0 = menuvt(menustr,callbackfn,disable)
% MENUVT Creates a vertical menu
%		MENUVT creates a vertical menu of pushbuttons. The
%		first parameter defines the labels of the pushbuttons
%		and has to be a cell structure. The second parameter
%		is a string which specifies the callback routine when 
%		the buttons are pressed. This callback routine must
%		accept one parameter which is the index of the item.
%		The Tag of each menuitem is 'mi%d' where %d again
%		represents the index.

% Copyright (c) 1999 Satrajit S. Ghosh (satra@bu.edu)
% Department of Cognitive and Neural Systems, Boston University

% $Revision: 1.00$ $Date:Fri Oct  1 17:02:05 EDT 1999$

% Bug fixes

% Modifications

% Find any active menu and close it
obj = getmenuobj;
if (~(isempty(obj))),
   close(obj);
   disp('closing obj');
end;

% Disable keep. Not yet implemented
if (nargin<3),
   disable = 1;
else,
   disable = 0;
end;

% Position menu to the right of existing figure
obj = findobj('Tag','Vtcalcs');
pos = get(obj,'position');

h0 = figure(...
   'units','normalized',...
   'menubar','none',...
   'Tag','menuvt',...
   'Resize','off',...
   'numbertitle','off',...
   'position',[pos(1)+pos(3)+0.01 0.4 0.2 0.4]);
%   'Windowstyle','modal',...

wd = 0.95/length(menustr);

% Kludge to remove calculate option
if (strcmp(callbackfn,'pccb') | strcmp(callbackfn,'tccb'))
   idx = 1;
else
   idx = 2;
end;

for i=idx:length(menustr),
   h1 = uicontrol(h0,...
      'units','normalized',...
      'position',[0.01 0.9-(i-1)*wd 0.98 0.8*wd],...
      'string',menustr(i),...
      'HorizontalAlignment','left',...
      'Busyaction','cancel',...
      'Tag',sprintf('mi%d',i),...
      'Userdata',menustr(i),...
      'callback',[callbackfn sprintf('(%d)',i)]);
   if (i==3 & disable),
      set(h1,'Enable','off');
   end;
end;
