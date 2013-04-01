function plot_FBA(F,B,A)

obj = findobj('Tag','Formant');
axes(obj);
cla;
nF = length(F);

axis ij;
axis([0 40 0 60]);
hold on;
text(12,2,'Freq');
text(22,2,'BW');
text(32,2,'Amp');

k = (60-2)/8;
for i=1:nF,
   text(2,8+k*(i-1),sprintf('F%d',i));
   text(12,8+k*(i-1),num2str(round(F(i))));
   text(22,8+k*(i-1),num2str(round(B(i))));
   text(32,8+k*(i-1),num2str(round(A(i))));
end;


hold off;