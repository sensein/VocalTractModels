function plot_lam(P1,P2)

obj = findobj('Tag','AM');
axes(obj);
cla;
hold on
ivtx = P1(1,:);
ivty = P1(2,:);
evtx = P1(3,:);
evty = P1(4,:);
lip_w = P2(1);
lip_h = P2(2);

axis([0 20 0 20]);
plot(ivtx,ivty);
plot(evtx,evty);

lip_x0 = .25*20;
lip_y0 = .5*20;
theta = 2*pi/20;
lipx = lip_w*cos([1:22]*theta)+lip_x0;
lipy = lip_h*sin([1:22]*theta)+lip_y0;

plot(lipx,lipy);
hold off
