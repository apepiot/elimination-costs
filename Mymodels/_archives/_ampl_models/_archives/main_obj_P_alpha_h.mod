

model defConstraints_rho_h.mod;

s.t. cobj: Prevalence_H*N >= bnd_sup_0*N;
# Objective function to be minimized.
maximize objective: rho_h;