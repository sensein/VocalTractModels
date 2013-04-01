function pccb(idx)
%	PCCB Callback function for phsycial constants change
%		This function is called whenever the menuitems
%		associated with the phsycial constants change are activated.
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

switch (idx),
case 0,
   for i=1:5,
      str = sprintf('mi%d',i);
      obj = findobj('Tag',str);
      str = get(obj,'Userdata');
      switch (i),
      case 1,
         set(obj,'string',[str{:} ' = ' sprintf('%7.5f',data.PC.ro)]);
      case 2,
         set(obj,'string',[str{:} ' = ' sprintf('%7.1f',data.PC.c)]);
      case 3,
         set(obj,'string',[str{:} ' = ' sprintf('%7.1f',data.PC.wall_resi)]);
      case 4,
         set(obj,'string',[str{:} ' = ' sprintf('%7.1f',data.PC.wall_mass)]);
      case 5,
         set(obj,'string',[str{:} ' = ' sprintf('%7.1f',data.PC.wall_comp)]);
      end;
   end;
case 1,
   newval = getnval(1.5e-3,0.9e-3,data.PC.ro);
   data.PC.ro = newval;
   set(obj1,'Userdata',data);
   pccb(0);
case 2,
   newval = getnval(3.6e+4,3.3e+4,data.PC.c);
   data.PC.c = newval;
   set(obj1,'Userdata',data);
   pccb(0);
case 3,
   newval = getnval(2000,1000,data.PC.wall_resi);
   data.PC.wall_resi = newval;
   set(obj1,'Userdata',data);
   pccb(0);
case 4,
   newval = getnval(2,1,data.PC.wall_mass);
   data.PC.wall_mass = newval;
   set(obj1,'Userdata',data);
   pccb(0);
case 5,
   newval = getnval(4e+5,2e+5,data.PC.wall_comp);
   data.PC.wall_comp = newval;
   set(obj1,'Userdata',data);
   pccb(0);
case 6, 
   obj=findobj('Tag','menuvt'); 
   close(obj);
end;

   