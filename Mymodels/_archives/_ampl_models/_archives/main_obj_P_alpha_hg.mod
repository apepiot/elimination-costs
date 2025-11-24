
model defConstraints_rho_hg.mod;

s.t. cobj: Prevalence_HG*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_hg;