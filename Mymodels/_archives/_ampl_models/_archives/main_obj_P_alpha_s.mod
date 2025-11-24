
model defConstraints_rho_s.mod;

s.t. cobj: Prevalence_S*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_s;