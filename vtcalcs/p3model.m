function p3model

obj = findobj('Tag','Vtcalcs');
data = get(obj,'Userdata');
obj = findobj('Tag','AM');
axes(obj);
cla
set(obj,'visible', 'off');

menustr = {...
      'Calculate' 'Synthesize' 'Keep (save)' ...
      'Cross area, At (cm2)', ...
      'Position, Xt (cm)', ...
      'Lip aperture, Al (cm)', ...
      'Nasal coupling (cm2)', ...
      'Close'};

h = menuvt(menustr,'p3cb');
p3cb(0);
%uiwait(h);
