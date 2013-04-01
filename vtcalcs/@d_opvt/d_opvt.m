function op = d_opvt(varargin)
%D_OPVT Constructor for diva operator for vocal tract object
%   The operator does two things regarding a vocal tract.
%   a. Plots the shape of the vocal tract
%   b. Synthesizes the sound corresponding to the vocal tract
%      shape

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/d_opvt.m 3     10/18/01 2:45p Satra $

% $NoKeywords: $

% Setup globals
global RELEASE

switch nargin
case 0
    % if no input arguments, create a default object  
    % based on the Maeda vocal tract that has 7 dimensions
    op.idim = 7;
    op = class(op,'d_opvt');
case 1
    % if single argument of class operator, return it
    if (isa(varargin{1},'d_opvt'))
        op = varargin{1};
   else
        error('Wrong argument type')
    end 
otherwise
    error('Wrong number of input arguments')
end
