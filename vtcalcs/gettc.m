function tccfg = gettc(TCcfg)
%	GETPC Converts the data structure for tract configuration
%		into a row vector

tccfg = [...
      TCcfg.rad_boundary TCcfg.wall...
      TCcfg.nasal_tract TCcfg.glt_boundary];
