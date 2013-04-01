function Makefile(platform)

if (strcmp(lower(platform),'unix'))
   disp('Compiling on Unix');
   mex -I. UTgetdata.c vtf_lib.c vsyn_lib.c vtt_lib.c complex.c ...
   plot_lib.c
   
   mex -I. T2getdata.c vtf_lib.c vsyn_lib.c vtt_lib.c complex.c ...
   plot_lib.c
   
   mex -I. P3getdata.c vtf_lib.c vsyn_lib.c vtt_lib.c complex.c ...
   plot_lib.c
   
   mex -I. AMgetdata.c vtf_lib.c vsyn_lib.c vtt_lib.c complex.c ...
   plot_lib.c lam_lib.c
   
   mex -I. AS2F.c vtf_lib.c vsyn_lib.c vtt_lib.c complex.c ...
   plot_lib.c lam_lib.c
elseif (strcmp(lower(platform),'windows'))
   disp('Compiling on Windows');
   mex -I. -DWINDOWS_COMPILE UTgetdata.c vtf_lib.c vsyn_lib.c vtt_lib.c complex.c ...
   plot_lib.c
   
   mex -I. -DWINDOWS_COMPILE T2getdata.c vtf_lib.c vsyn_lib.c vtt_lib.c complex.c ...
   plot_lib.c
   
   mex -I. -DWINDOWS_COMPILE P3getdata.c vtf_lib.c vsyn_lib.c vtt_lib.c complex.c ...
   plot_lib.c
   
   mex -I. -DWINDOWS_COMPILE AMgetdata.c vtf_lib.c vsyn_lib.c vtt_lib.c complex.c ...
   plot_lib.c lam_lib.c
   
   mex -I. -DWINDOWS_COMPILE AS2F.c vtf_lib.c vsyn_lib.c vtt_lib.c complex.c ...
   plot_lib.c lam_lib.c
else,
   error('Incorrect platform specified');
end;
