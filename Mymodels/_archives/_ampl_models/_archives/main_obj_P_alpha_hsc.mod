
model defConstraints_rho_hsc.mod;

s.t. cobj: Prevalence_HSC*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_hsc;