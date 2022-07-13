function utcb(idx)
%	UTCB Callback function for unifrom tube model
%		This function is called whenever the menuitems
%		associated with the uniform tube model are activated.
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
         set(obj,'string',[str{:} ' = ' sprintf('%6.2f',data.UTmpar.area)]);
      case 5,
         set(obj,'string',[str{:} ' = '  sprintf('%6.1f',data.UTmpar.length)]);
      case 6,
         set(obj,'string',[str{:} ' = ' sprintf('%4d',data.UTmpar.nss)]);
      case 7,
         set(obj,'string',[str{:} ' = ' sprintf('%5.1f',data.UTmpar.anc)]);
      end;
   end;
   utcb(1);
case 1,
   [Af,Tfm,Tff,F,B,A] = UTgetdata(...
      gettc(data.TC),...
      getpc(data.PC),...
      [data.UTmpar.area,data.UTmpar.length,...
         data.UTmpar.nss,data.UTmpar.anc]);
   plot_af(Af);
   plot_tf(Tff,Tfm);
   plot_FBA(F,B,A);
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
   newval = getnval(15,0,data.UTmpar.area);
   data.UTmpar.area = newval;
   set(obj1,'Userdata',data);
   utcb(0);
case 5,
   newval = getnval(20,5,data.UTmpar.length);
   data.UTmpar.length = newval;
   set(obj1,'Userdata',data);
   utcb(0);
case 6,
   newval = getnval(60,5,data.UTmpar.nss);
   data.UTmpar.nss = fix(newval);
   set(obj1,'Userdata',data);
   utcb(0);
case 7,
   newval = getnval(3,0,data.UTmpar.anc);
   if (data.TC.nasal_tract == 0),
      newval = 0;
   end;
   data.UTmpar.anc = newval;
   set(obj1,'Userdata',data);
   utcb(0);
case 8, 
   obj=findobj('Tag','menuvt'); 
   close(obj);
end;

   