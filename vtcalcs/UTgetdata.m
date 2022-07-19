function [af,tfm,tff,f,b,a] = UTgetdata(TC,PC,UTpar)
%	UTGETDATA This function computes the uniform tube transfer 
%				 function (Implemented as a mex function)
%
%		UTGETDATA calculates the area functions, the transfer
%		function, the formants, bandwidths and amplitudes for
%		the given set of uniform tube parameters.
%
%		Input Parameters:
%		TC		: tract configuration from gettc
%		PC		: physical constants from getpc
%		UTpar	: Uniform tube model parameters 
%				  vector size(1,4)
%		Note	: data is defined in vtsimul
%
%		Output Parameters:
%		af	:	area functions
%		tfm	:	transfer function magnitude
%       tff :   transfer function frequency
%		f	:	formants
%		b	:	bandwidths
%		a	:	amplitudes
%
%		Usage:
%		[Af,Tfm,Tff,F,B,A] = UTgetdata(...
%	      gettc(data.TC),...
%	      getpc(data.PC),...
%			[data.UTmpar.area,data.UTmpar.length,...
%        data.UTmpar.nss,data.UTmpar.anc]);

% Copyright (c) 1999 Satrajit S. Ghosh (satra@bu.edu)
% Department of Cognitive and Neural Systems, Boston University

% $Revision: 1.00$ $Date:Fri Oct  1 17:02:05 EDT 1999$

% Bug fixes

% Modifications
