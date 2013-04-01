% The following file sets up acceptable ranges for:
%   - Formants, their Bandwidths and Amplitudes

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/private/diva_data.m 2     10/18/01 2:45p Satra $

% $NoKeywords: $

% Setup globals
global RELEASE

minF1 = 150;
maxF1 = 700;
minF2 = 600;
maxF2 = 2400;
minF3 = 1500;
maxF3 = 3500;

minB1 = 0;
maxB1 = 300;
minB2 = 0;
maxB2 = 450;
minB3 = 0;
maxB3 = 950;

minA1 = -10;
maxA1 = 35;
minA2 = -10;
maxA2 = 35;
minA3 = -10;
maxA3 = 30;


%RF = [minF1 maxF1;minF2 maxF2; minF3 maxF3];

minmaxF = [minF1 maxF1;minF2 maxF2; minF3 maxF3];
rangeF = minmaxF(:,2)-minmaxF(:,1);

minmaxB = [minB1 maxB1;minB2 maxB2; minB3 maxB3];
rangeB = minmaxB(:,2)-minmaxB(:,1);

minmaxA = [minA1 maxA1;minA2 maxA2; minA3 maxA3];
rangeA = minmaxA(:,2)-minmaxA(:,1);

minF0 = 80;
maxF0 = 220;

%minratio = diva_f2ratio([minF1;minF2;minF3],maxF0);
%minR1 = minratio(1);
%minR2 = minratio(2);
%minR3 = minratio(3);
%maxratio = diva_f2ratio([maxF1;maxF2;maxF3],minF0);
%maxR1 = maxratio(1);
%maxR2 = maxratio(2);
%maxR3 = maxratio(3);

minR1 = 0.04; minR2 = 0.03; minR3 = 0.02;
maxR1 = 0.9; maxR2 = 1.3; maxR3 = 0.7;
%minR1 = -1; minR2 = -1; minR3 = -1;
%maxR1 = 1.5; maxR2 = 1.5; maxR3 = 1.5;
