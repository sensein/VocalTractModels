function synth_vowel(data)
% SYNTH_VOWEL Synthesizes a vowel from formants and bandwidths
%		This function utilizes the track_synth function from the
%		Track Draw program to synthesize a vowel.

o = ones(40,1);
for i =1:5,
   F(:,i) = [o*data.F(i)];
   B(:,i) = [o*data.B(i)];
end;
AV=[o*60];
F0=[o*100];
y=track_synth(...
   200,5,16000,...
   F(:,1),F(:,2),F(:,3),F(:,4),F(:,5),...
   AV,F0,120,...
   B(:,1),B(:,2),B(:,3),B(:,4),B(:,5));
soundsc(y,16000,16);
