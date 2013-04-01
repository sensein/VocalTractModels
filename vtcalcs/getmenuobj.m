function obj = getmenuobj
% GETMENUOBJ Finds an active menu and returns it. The function
%		returns empty otherwise

% Copyright (c) 1999 Satrajit S. Ghosh (satra@bu.edu)
% Department of Cognitive and Neural Systems, Boston University

% $Revision: 1.00$ $Date:Fri Oct  1 17:02:05 EDT 1999$

% Bug fixes

% Modifications


obj = findobj('Tag','menuvt');
if isempty(obj),
	obj = findobj('Tag','menuas2f');
end;
if isempty(obj),
	obj = findobj('Tag','menulam');
end;
