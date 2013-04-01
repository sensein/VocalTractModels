function [af,tf,f,b,a] = T2getdata(TC,PC,T2par)
%	T2GETDATA This function computes the two tube transfer 
%				 function (Implemented as a mex function)
%
%		T2GETDATA calculates the area functions, the transfer
%		function, the formants, bandwidths and amplitudes for
%		the given set of two tube parameters.
%
%		Input Parameters:
%		TC		: tract configuration from gettc
%		PC		: physical constants from getpc
%		T2par	: Uniform tube model parameters 
%				  vector size(1,5)
%		Note	: data is defined in vtsimul
%
%		Output Parameters:
%		af	:	area functions
%		tf	:	transfer function
%		f	:	formants
%		b	:	bandwidths
%		a	:	amplitudes
%
%		Usage:
%		[Af,Tf,F,B,A] = T2getdata(...
%			gettc(data.TC),...
%			getpc(data.PC),...
%			[data.T2mpar.A1,data.T2mpar.x1,...
%        data.T2mpar.A2,data.T2mpar.x2,data.T2mpar.anc]);

% Copyright (c) 1999 Satrajit S. Ghosh (satra@bu.edu)
% Department of Cognitive and Neural Systems, Boston University

% $Revision: 1.00$ $Date:Fri Oct  1 17:02:05 EDT 1999$

% Bug fixes

% Modifications
