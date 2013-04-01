function [af,tf,f,b,a,p1,p2] = AMgetdata(TC,PC,AMpar)
%	AMGETDATA This function computes the LAM transfer function
%		(Implemented as a mex function)
%
%		AMGETDATA calculates the area functions, the transfer
%		function, the formants, bandwidths and amplitudes for
%		the given set of LAM parameters. In addition it also 
%		returns two varaibles p1 and p2 which can be used to 
%		plot the outline of the vocal tract using the function
%		plot_lam.
%
%		Input Parameters:
%		TC		: tract configuration from gettc
%		PC		: physical constants from getpc
%		AMpar	: LAM parameters vector size(1,8)
%		Note	: data is defined in vtsimul
%
%		Output Parameters:
%		af	:	area functions
%		tf	:	transfer function
%		f	:	formants
%		b	:	bandwidths
%		a	:	amplitudes
%		p1	:	plotting parameters for plot_lam
%		p2 :	plotting parameters for plot_lam
%
%		Usage:
%		[Af,Tf,F,B,A,P1,P2] = AMgetdata(...
%		      gettc(data.TC),...
%		      getpc(data.PC),...
%		      [data.AMpar.ampar,data.AMpar.anc]);

% Copyright (c) 1999 Satrajit S. Ghosh (satra@bu.edu)
% Department of Cognitive and Neural Systems, Boston University

% $Revision: 2 $ $Date: 10/18/01 2:45p $

% Bug fixes

% Modifications
