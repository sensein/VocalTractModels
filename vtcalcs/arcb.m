function arcb(idx)
%	ARCB Callback function for Articulator model
%		This function is called whenever the menuitems
%		associated with the Articulatory model are activated.
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

% Names of the vowels
% At present they can only be accessed by numbers
vowelidx = {...
        'None' 'iy' 'ey' 'eh' 'ah' 'aa' 'ao' ...
        'oh' 'uw' 'iw' 'ew' 'oe'};

% Vowel definitions in terms of LAM parameters
vowelpar = [...
        0, 0, 0, 0, 0, 0, 0;...
        0.5, -2.0, 1.0, -2.0,  1.0, -1.0, 0.0;...
        0.0, -1.0, 1.0, -2.0,  1.0, -1.0, 0.0;...
        -1.0,  0.0, 1.0, -2.0,  1.0, -0.5, 0.0;...
        -1.5,  0.5, 0.0, -0.5,  0.5, -0.5, 0.0;...
        -1.5,  2.0, 0.0, -0.5,  0.5, -0.5, 0.0;...
        -0.4,  3.0, 1.5,  0.0, -0.3,  0.0, 0.0;...
        -.7,  3.0, 1.5,  0.0, -0.6,  0.0, 0.0;...
        0.5,  2.0, 1.5, -2.0, -1.0,  1.5, 0.0;...
        0.5, -1.0, 1.0, -2.0, -0.5,  1.0, 0.0;...
        0.0, -0.2, 1.0, -1.5, -0.25, 0.5, 0.0;...
        -1.0, -0.5, 0.5, -2.0,  0.2, -0.5, 0.0];

slst = 6; % Start of slider
% This is where all the actions are carried out
switch (idx),
case 0,
    for i=4:13,
        str = sprintf('mi%d',i);
        obj = findobj('Tag',str);
        switch (i),
        case 4,
            str = get(obj,'Userdata');
            vname = vowelidx(data.AMpar.vowel+1);
            set(obj,'string',[str{:} ' = ' vname{:}]);
        case 5,
            str = get(obj,'Userdata');
            set(obj,'string',[str{:} ' = ' sprintf('%5.2f',data.AMpar.anc)]);
        case {7,8,9,10,11,12,13},
            str1 = sprintf('mistr%d',i);
            obj3 = findobj('Tag',str1);
            str3 = get(obj3,'Userdata');
            str2 = sprintf('misl%d',i);
            obj4 = findobj('Tag',str2);
            set(obj3,'string',[str3{:} sprintf('[%1.2f]',data.AMpar.ampar(i-slst))]);
            set(obj4,'value',data.AMpar.ampar(1,i-slst)*100+301);
        end;
    end;
    arcb(1);
case 1,	%	Calculate 
    [Af,Tfm,Tff,F,B,A,P1,P2] = AMgetdata(...
        gettc(data.TC),...
        getpc(data.PC),...
        [data.AMpar.ampar,data.AMpar.anc]);
    plot_af(Af);
    plot_FBA(F,B,A);
    plot_tf(Tff,Tfm);
    vt = d_opvt;
    obj = findobj('Tag','AM');
    axes(obj);
    plot(vt,data.AMpar.ampar(:),1);
    %plot_lam(P1,P2);
    data.F = [];data.B = [];data.A=[];
    data.F = F;data.B = B;data.A=A;
    set(obj1,'Userdata',data);
case 2,	% Synthesize
    if (length(data.F)>4),
        %synth_vowel(data);
        vt = d_opvt;
        synth1(vt,data.AMpar.ampar(:));
    else,
        disp('Not enough formants for synthesis');
    end;
case 4,	% Vowel
    newval = getnval(11,0,data.AMpar.vowel);
    data.AMpar.vowel = fix(newval);
    data.AMpar.ampar = vowelpar(newval+1,:);
    set(obj1,'Userdata',data);
    arcb(0);
case 5,	% Nasal coupling
    newval = getnval(3,0,data.AMpar.anc);
    if (data.TC.nasal_tract == 0),
        newval = 0;
    end;
    data.AMpar.anc = newval;
    set(obj1,'Userdata',data);
    arcb(0);
case 6, 	% close
    obj=findobj('Tag','menulam'); 
    close(obj);
case {7,8,9,10,11,12,13}, % All the sliders
    str2 = sprintf('misl%d',idx);
    obj2 = findobj('Tag',str2);
    newval = ((get(obj2,'value')-1)/100)-3;
    data.AMpar.ampar(1,idx-slst) = newval;
    set(obj1,'Userdata',data);
    arcb(0);
end;

