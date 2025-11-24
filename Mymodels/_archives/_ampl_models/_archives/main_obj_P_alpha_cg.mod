
model defConstraints_rho_cg.mod;

s.t. cobj: Prevalence_CG*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_cg;