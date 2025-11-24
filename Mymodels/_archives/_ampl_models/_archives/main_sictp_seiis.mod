model defConstraints_sict_seiis.mod;

s.t. c30: Prevalence_H >=bnd_sup_0;
s.t. c31: Prevalence_X >=bnd_sup_0;


# Objective function to be minimized.
minimize objective: 1;

