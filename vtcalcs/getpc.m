function pccfg = getpc(PCcfg)
%	GETPC Converts the data structure for physical constants
%		into a row vector

pccfg = [...
      PCcfg.ro PCcfg.c PCcfg.wall_resi ...
      PCcfg.wall_mass PCcfg.wall_comp];
