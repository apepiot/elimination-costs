
model defConstraints_sictp_seiis_hX.mod;

#s.t. cobj: Prevalence_HX*N >= bnd_sup_0*N;

s.t. cobj1: Prevalence_H*N >= bnd_sup_0*N;
s.t. cobj2: Prevalence_X*N >= bnd_sup_0*N;


# Objective function to be minimized.
maximize objective: rho_hX;
