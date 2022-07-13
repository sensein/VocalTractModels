function [af,tfm,tff,f,b,a] = P3getdata(TC,PC,P3par)
%	P3GETDATA This function computes the 3 parameter transfer 
%				 function (Implemented as a mex function)
%
%		P3GETDATA calculates the area functions, the transfer
%		function, the formants, bandwidths and amplitudes for
%		the given set of 3 parameters.
%
%		Input Parameters:
%		TC		: tract configuration from gettc
%		PC		: physical constants from getpc
%		P3par	: Three parameter model parameters 
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
%		[Af,Tfm,Tff,F,B,A] = P3getdata(...
%				gettc(data.TC),...
%				getpc(data.PC),...
%				[data.P3mpar.At, data.P3mpar.Xt,...
%				data.P3mpar.Al, data.P3mpar.anc]);

% Copyright (c) 1999 Satrajit S. Ghosh (satra@bu.edu)
% Department of Cognitive and Neural Systems, Boston University

% $Revision: 1.00$ $Date:Fri Oct  1 17:02:05 EDT 1999$

% Bug fixes

% Modifications
