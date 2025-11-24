
model defConstraints_rho_hc.mod;

s.t. cobj: Prevalence_HC*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_hc;
