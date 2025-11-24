
model defConstraints_rho_sc.mod;

s.t. cobj: Prevalence_SC*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_sc;