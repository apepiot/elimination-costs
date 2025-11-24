

model defConstraints_rho_c.mod;

s.t. cobj: Prevalence_C*N >= bnd_sup_0*N;
# Objective function to be minimized.
maximize objective: rho_c;