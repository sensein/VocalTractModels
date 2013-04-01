function as2fcb(idx)
%	AS2FCB Callback function for Areafunction model
%		This function is called whenever the menuitems
%		associated with the Area function model are activated.
%		The identification of the menu item is provided by the
%		input variable idx.

% Copyright (c) 1999 Satrajit S. Ghosh (satra@bu.edu)
% Department of Cognitive and Neural Systems, Boston University

% $Revision: 1.00$ $Date:Fri Oct  1 17:02:05 EDT 1999$

% Bug fixes

% Modifications
% 10.2.99 Automatically calculates and updates

% Get the main data structure from the parent figure.
obj1 = findobj('Tag','Vtcalcs');
data = get(obj1,'Userdata');

slst = 8; % slider start
switch (idx),
case 0,
   for i=4:slst-1,
      str = sprintf('mi%d',i);
      obj2 = findobj('Tag',str);
      str = get(obj2,'Userdata');
      switch (i),
      case 4,
         set(obj2,'string',[str{:} ' = ' sprintf('%2d',data.AS2Fpar.nss)]);
      case 5,
         set(obj2,'string',[str{:} ' = ' sprintf('%5.1f',data.AS2Fpar.length)]);
      case 6,
         set(obj2,'string',[str{:} ' = ' sprintf('%2d',data.AS2Fpar.nts)]);
      case 7,
         set(obj2,'string',[str{:} ' = ' sprintf('%5.1f',data.AS2Fpar.anc)]);
      end;
   end;
   for i=(slst+1):(slst+data.AS2Fpar.nss),
      str1 = sprintf('mistr%d',i-slst);
      obj3 = findobj('Tag',str1);
      str2 = sprintf('misl%d',i-slst);
      obj4 = findobj('Tag',str2);
      set(obj3,'string',sprintf('%d[%2.2f]cm2',i-slst,data.AS2Fpar.af(1,i-slst)));
      set(obj4,'value',data.AS2Fpar.af(1,i-slst)*100);
   end;
   as2fcb(1);
case 1,
   [Af,Tf,F,B,A] = AS2F(...
      gettc(data.TC),...
      getpc(data.PC),...
      data.AS2Fpar.af,...
      [data.AS2Fpar.nts data.AS2Fpar.anc]);
   plot_af(Af);
   plot_af(data.AS2Fpar.af,1);
   plot_tf(Tf);
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
   newval = getnval(30,5,data.AS2Fpar.nss);
   if (newval ~= data.AS2Fpar.nss),
     data.AS2Fpar.nss = fix(newval);
     data.AS2Fpar.af = [];
     data.AS2Fpar.af(1,:) = 6*ones(1,newval);
     data.AS2Fpar.af(2,:) = (data.AS2Fpar.length/newval)*ones(1,newval);
     set(obj1,'Userdata',data);
     obj=findobj('Tag','menuas2f'); 
     close(obj);
     as2fmodel;
  end;
case 5,
   newval = getnval(20,5,data.AS2Fpar.length);
   data.AS2Fpar.af(2,:) = (newval/data.AS2Fpar.nss)*ones(1,data.AS2Fpar.nss);
   data.AS2Fpar.length = newval;
   set(obj1,'Userdata',data);
   as2fcb(0);
case 6,
   newval = getnval(60,5,data.AS2Fpar.nts);
   data.AS2Fpar.nts = fix(newval);
   set(obj1,'Userdata',data);
   as2fcb(0);
case 7,
   newval = getnval(3,0,data.AS2Fpar.anc);
   if (data.TC.nasal_tract == 0),
      newval = 0;
   end;
   data.AS2Fpar.anc = newval;
   set(obj1,'Userdata',data);
   as2fcb(0);
case 8, 
   obj=findobj('Tag','menuas2f'); 
   close(obj);
end;

if (idx>slst) & (idx<=(slst+data.AS2Fpar.nss))
      str2 = sprintf('misl%d',idx-slst);
      obj2 = findobj('Tag',str2);
      newval = get(obj2,'value')/100;
      data.AS2Fpar.af(1,idx-slst) = newval;
	   set(obj1,'Userdata',data);
      as2fcb(0);
end;

