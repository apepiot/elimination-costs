
model defConstraints_rho_hsg.mod;

s.t. cobj: Prevalence_HSG*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_hsg;