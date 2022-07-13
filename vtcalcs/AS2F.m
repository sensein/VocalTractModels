function [af,tfm,tff,f,b,a] = AS2F(TC,PC,AF,misc)
%	AS2F This function computes the transfer function from
%		  a given set of area functions.
%		  (Implemented as a mex function)
%
%		AS2F calculates the area functions, the transfer
%		function, the formants, bandwidths and amplitudes for
%		the given set of area functions with individual lengths.
%
%		Input Parameters:
%		TC		: tract configuration from gettc
%		PC		: physical constants from getpc
%		AF		: Area functions vector size(2,n)
%				  The first row contains the areas and second row
%				  contains the lengths. n is the number of areas.
%		misc	: a vector of two elements. The first specifies
%				  the number of tube segments of equal length
%				  and the second is the nasal coupling area.
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
%		[Af,Tfm,Tff,F,B,A] = AS2F(...
%			gettc(data.TC),...
%			getpc(data.PC),...
%			data.AS2Fpar.af,...
%			[data.AS2Fpar.nts data.AS2Fpar.anc]);

% Copyright (c) 1999 Satrajit S. Ghosh (satra@bu.edu)
% Department of Cognitive and Neural Systems, Boston University

% $Revision: 1.00$ $Date:Fri Oct  1 17:02:05 EDT 1999$

% Bug fixes

% Modifications
