
model defConstraints_rho_scg.mod;

s.t. cobj: Prevalence_SCG*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_scg;