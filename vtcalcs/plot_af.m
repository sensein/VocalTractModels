function plot_af(Af,bhold)
if (nargin == 1),
   bhold = 0;
end;

obj = findobj('Tag','AreaFn');
axes(obj);

cl = {'r-'};
if (~bhold),   
   cla;
   axis([0 20 0 16]);
   cl = {'b-'};
end;

hold on;
tl = sum(Af(2,:)); % tube length
wd = Af(2,:)*(tl/sum(Af(2,:)));
start = 0;

plot([start wd(1)],[Af(1,1) Af(1,1)],cl{:});
start = start +wd(1);
for i=2:length(Af),
   plot([start start],[Af(1,i-1) Af(1,i)],cl{:});
   plot([start start+wd(i)],[Af(1,i) Af(1,i)],cl{:});
	start = start +wd(i);
end;

set(obj,'XTick',[0:2:20]);
set(obj,'XTickLabel',[0:2:20])
set(obj,'YTick',[0:2:16]);
set(obj,'YTickLabel',[0:2:16])
title('Area Function');
xlabel('Distance from glottis (cm)');
ylabel('Area (cm^2)');
hold off;

