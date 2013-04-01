function t2cb(idx)
%	T2CB Callback function for two tube model
%		This function is called whenever the menuitems
%		associated with the two tube model are activated.
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
   for i=4:7,
      str = sprintf('mi%d',i);
      obj = findobj('Tag',str);
      str = get(obj,'Userdata');
      switch (i),
      case 4,
         set(obj,'string',[str{:} ' = ' sprintf('%4.1f',data.T2mpar.A1)]);
      case 5,
         set(obj,'string',[str{:} ' = ' sprintf('%4.1f',data.T2mpar.x1)]);
      case 6,
         set(obj,'string',[str{:} ' = ' sprintf('%4.1f',data.T2mpar.A2)]);
      case 7,
         set(obj,'string',[str{:} ' = ' sprintf('%4.1f',data.T2mpar.x2)]);
      case 8,
         set(obj,'string',[str{:} ' = ' sprintf('%5.2f',data.UTmpar.anc)]);
      end;
   end;
   t2cb(1);
case 1,
   [Af,Tf,F,B,A] = T2getdata(...
      gettc(data.TC),...
      getpc(data.PC),...
      [data.T2mpar.A1,data.T2mpar.x1,...
         data.T2mpar.A2,data.T2mpar.x2,data.T2mpar.anc]);
   plot_af(Af);
   plot_FBA(F,B,A);
   plot_tf(Tf);
   data.F = [];data.B = [];data.A=[];
   data.F = F;data.B = B;data.A=A;
   set(obj1,'Userdata',data);
case 2,
   if (length(data.F)>4),
      synth_vowel(data);
   else,
      disp('Not enough formants for synthesis');
   end;
case 4,
   newval = getnval(10,0,data.T2mpar.A1);
   data.T2mpar.A1 = newval;
   set(obj1,'Userdata',data);
   t2cb(0);
case 5,
   newval = getnval(20-data.T2mpar.x2,0,data.T2mpar.x1);
   data.T2mpar.x1 = newval;
   set(obj1,'Userdata',data);
   t2cb(0);
case 6,
   newval = getnval(10,0,data.T2mpar.A2);
   data.T2mpar.A2 = newval;
   set(obj1,'Userdata',data);
   t2cb(0);
case 7,
   newval = getnval(20-data.T2mpar.x1,0,data.T2mpar.x2);
   data.T2mpar.x2 = newval;
   set(obj1,'Userdata',data);
   t2cb(0);
case 8,
   newval = getnval(3,0,data.T2mpar.anc);
   if (data.TC.nasal_tract == 0),
      newval = 0;
   end;
   data.T2mpar.anc = newval;
   set(obj1,'Userdata',data);
   t2cb(0);
case 9, 
   obj=findobj('Tag','menuvt'); 
   close(obj);
end;

   