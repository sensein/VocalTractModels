% SPLINT: Cubic spline interpolation
% allows for "decimation" as well as interpolation
% Usage: yi=splint(y,ni);
function yi=splint(y,ni);
n=length(y);
if ni<n, k=n; m=ni*n; else, k=1; m=ni; end;
yi=spline(1:n,y,1:(n-1)/(m-1):n);
yi=yi(1:k:m);
%clg; hold off; 
%axis([0 length(yi) min(yi) max(yi)]);
%plot(yi); hold on;
%plot(yi,'.'); 
%axis([0 length(y) min(y) max(y)]);
%plot(y,'o');
%hold off;
