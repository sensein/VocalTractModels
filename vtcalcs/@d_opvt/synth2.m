function [Fmts,LPCcoeff,idxValid] = synth2(vtm,varargin)
% SYNTH2    Synthesizes the output from the vocal tract
%   This function behaves very similar to synth1 except for its 
%   output. The function returns Fmts, LPC coefficients and a list of
%   valid articulator positions.
%   Please see synth1 for additional help

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/synth2.m 4     10/24/01 11:40a Satra $

% $NoKeywords: $

% Setup globals
global RELEASE

if nargin == 2,
    if isstruct(varargin{1}),
        y = diva_synth2(data);
    end;
end;

X = varargin{1};

if (nargin>2),
    data.fs = varargin{2};
else,
    data.fs = 10000;
end;

if nargin>3,
    data.duration = varargin{3};
else,
    data.duration = 300;
end;
tf = data.duration;

if nargin >4,
    data.Ag0 = varargin{4};
    data.AgP = varargin{5};
    data.F0 = varargin{6};
else,
    data.Ag0 = [0.  0.0   0;tf  0.0   1];
    data.AgP = [0.  0.1  0;tf  0.2  1];
    data.F0 = [0.  100.  0;tf  100.  1];
end

% file containing formant ranges
diva_data;

% Setup preemphasis parameters
preemp = -0.95;

% Iterate through each vocal tract configuration
for i=1:size(X,2),   
    [F,B,A,Af] = doAM(X(:,i));
    if length(F)<3,
        F(length(F)+1:3) = NaN;
        F = F(:);
    end;
    
    data.Af = [];data.dx = [];
    data.Af(:,1) = Af(1,:)';
    data.Af(:,2) = Af(1,:)';
    data.dx = mean(Af(2,:));
    
    % Determine the acoustic signal
    out = diva_synth2(data);
    
    % Preemphasize
    y = filter([1 preemp],1,out);    
    % Get the LPC coefficients
    [A,G]=lpc(hamming(length(y)).*y(:),14);
    A = real(A);
    
    % find the formants
    if length(find(isnan(A)))==0,
        F1 = frmnts1(A(:),data.fs);
    else,
        F1 = zeros(3,1);
    end;
    
    F1 = F1(1:3);
    F1 = F1(:);
    Fmts(:,i) = F1;
    LPCcoeff(:,i) = A(:);
    
    % Determine validity of formants
    idxValid(i) = 0;
    if ( ...
            sum(abs(F1-F(1:3))<[25;50;100])==3 & ...
            (F(1)>minF1) & (F(1)<maxF1) & ...
            (F(2)>minF2) & (F(2)<maxF2) & ...
            (F(3)>minF3) & (F(3)<maxF3)),
        idxValid(i) = 1;
    end;
end