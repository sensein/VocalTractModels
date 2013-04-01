function p3cb(idx)
%	P3CB Callback function for 3 parameter model
%		This function is called whenever the menuitems
%		associated with the 3 parameter model are activated.
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
         set(obj,'string',[str{:} ' = ' sprintf('%6.2f',data.P3mpar.At)]);
      case 5,
         set(obj,'string',[str{:} ' = '  sprintf('%6d',data.P3mpar.Xt)]);
      case 6,
         set(obj,'string',[str{:} ' = ' sprintf('%6.2f',data.P3mpar.Al)]);
      case 7,
         set(obj,'string',[str{:} ' = ' sprintf('%6.2f',data.P3mpar.anc)]);
      end;
   end;
   p3cb(1);
case 1,
   [Af,Tf,F,B,A] = P3getdata(...
      gettc(data.TC),...
      getpc(data.PC),...
      [data.P3mpar.At, data.P3mpar.Xt,...
         data.P3mpar.Al, data.P3mpar.anc]);
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
   newval = getnval(8,0,data.P3mpar.At);
   data.P3mpar.At = newval;
   set(obj1,'Userdata',data);
   p3cb(0);
case 5,
   newval = getnval(13,3,data.P3mpar.Xt);
   data.P3mpar.Xt = newval;
   set(obj1,'Userdata',data);
   p3cb(0);
case 6,
   newval = getnval(4,0,data.P3mpar.Al);
   data.P3mpar.Al = newval;
   set(obj1,'Userdata',data);
   p3cb(0);
case 7,
   newval = getnval(3,0,data.P3mpar.anc);
   if (data.TC.nasal_tract == 0),
      newval = 0;
   end;
   data.P3mpar.anc = newval;
   set(obj1,'Userdata',data);
   p3cb(0);
case 8, 
   obj=findobj('Tag','menuvt'); 
   close(obj);
end;

   