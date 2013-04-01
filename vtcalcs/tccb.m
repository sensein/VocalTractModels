function tccb(idx)
%	TCCB Callback function for tract config change
%		This function is called whenever the menuitems
%		associated with the tract config change are activated.
%		The identification of the menu item is provided by the
%		input variable idx.

% Copyright (c) 1999 Satrajit S. Ghosh (satra@bu.edu)
% Department of Cognitive and Neural Systems, Boston University

% $Revision: 1.00$ $Date:Fri Oct  1 17:02:05 EDT 1999$

% Bug fixes

% Modifications

% Get the main data structure from the parent figure.
obj1 = findobj('Tag','Vtcalcs');
data = get(obj1,'Userdata');

rad = { 'OFF' 'ON'};
wall = { 'RIGID' 'YIELDING'};
nasal = { 'OFF' 'ON'};
glottis = {'CLOSE' 'OPEN'};

switch (idx),
case 0,
   for i=1:4,
      str = sprintf('mi%d',i);
      obj = findobj('Tag',str);
      str = get(obj,'Userdata');
      switch (i),
      case 1,
         vn = rad(data.TC.rad_boundary+1); 
         set(obj,'string',[str{:} ' = ' vn{:}]);
      case 2,
         vn = wall(data.TC.wall+1); 
         set(obj,'string',[str{:} ' = ' vn{:}]);
      case 3,
         vn = nasal(data.TC.nasal_tract+1); 
         set(obj,'string',[str{:} ' = ' vn{:}]);
      case 4,
         vn = glottis(data.TC.glt_boundary+1); 
         set(obj,'string',[str{:} ' = ' vn{:}]);
      end;
   end;
case 1,
   data.TC.rad_boundary = ~data.TC.rad_boundary;
   set(obj1,'Userdata',data);
   tccb(0);
case 2,
   data.TC.wall = ~data.TC.wall;
   set(obj1,'Userdata',data);
   tccb(0);
case 3,
   data.TC.nasal_tract = ~data.TC.nasal_tract;
   set(obj1,'Userdata',data);
   tccb(0);
case 4,
   data.TC.glt_boundary = data.TC.glt_boundary;
   set(obj1,'Userdata',data);
   tccb(0);
case 5, 
   obj=findobj('Tag','menuvt'); 
   close(obj);
end;

   