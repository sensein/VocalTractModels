function plot_tf(Tf)

obj = findobj('Tag','Spectrum');
axes(obj);
cla;
hold on
f = linspace(0,5,length(Tf));
plot(f,Tf);
axis([0 5 -40 40]);
title('(U_m+U_n)/U_g Transfer ratio');
xlabel('Frequency (KHz)');
ylabel('Relative Mag. (dB)');
hold off
