
model defConstraints_rho_hscg.mod;

s.t. cobj: Prevalence_HSCG*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_hscg;