function plot_tf(Tff,Tfm)

obj = findobj('Tag','Spectrum');
axes(obj);
cla;
hold on
plot(Tff,Tfm);
axis([0 5000 -40 40]);
title('(U_m+U_n)/U_g Transfer ratio');
xlabel('Frequency (Hz)');
ylabel('Relative Mag. (dB)');
hold off
