function y = synth1(vtm,varargin)
% SYNTH1    Synthesizes the output from the vocal tract
%   The function has two different forms:
%       1. Y = synth1(vtm,data);
%           In this case the input data contains the following
%           subfields: Ag0, AgP, F0, Af, dx, fs, duration
%       2. Y = synth1(vtm,X,fs,tf,Ag0,AgP,F0);
%           In this case, the fields of data are passed individually
%           to the funtion
%   Please consult Maeda's paper for a better understanding of the inputs
%
%   The typical form of the parameters Ag0, AgP and F0 are of the form:
%   [timestamp1 value1 interpolation_scheme1;timestamp2 value2 interpolation_scheme2; ...]
%   where interpolation_scheme takes values: 
%           0 - set
%           1 - linear interpolation
%           2 - exponential interpolation

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/synth1.m 5     10/19/01 11:19a Satra $

% $NoKeywords: $

% Setup globals
global RELEASE

% If there are only two inputs, assume that the second input is
% the struct containing information about the production
if nargin == 2,
    if isstruct(varargin{1}),
        data = varargin{1}
        y = diva_synth2(data);
        if nargout==0,
            soundsc(y,data.fs);
        end;
        return;
    end;
end;

X = varargin{1};

% Setup the data structure
% a. sampling rate
if (nargin>2),
    data.fs = varargin{2};
else,
    data.fs = 10000;
end;

% b. duration of utterance
if nargin>3,
    data.duration = varargin{3};
else,
    data.duration = 300;
end;
tf = data.duration;

% c. Fields, ag0,agP,F0
if nargin >4,
    data.Ag0 = varargin{4};
    data.AgP = varargin{5};
    data.F0 = varargin{6};
else,
    data.Ag0 = [0.  0.0   0;tf  0.0   1];
    data.AgP = [0.  0.1  0;tf  0.2  1];
    data.F0 = [0.  100.  0;tf  100.  1];
end

% Iterate through inputs.
for i=1:size(X,2),   
    [F,B,A,Af] = doAM(X(:,i));
    data.Af = [];data.dx = [];
    data.Af(:,1) = Af(1,:)';
    data.Af(:,2) = Af(1,:)';
    data.dx = mean(Af(2,:));
    
    out = diva_synth2(data);
    
    % set output to zero if area function is too small
    if min(data.Af(:,1))<1e-1,
        y(:,i) = zeros(size(out));
    else,
        y(:,i) = out;
    end;    
    if nargout==0,
        soundsc(y,data.fs);
    end;
end