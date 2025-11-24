

model defConstraints_rho_g.mod;

s.t. cobj: Prevalence_G*N >= bnd_sup_0*N;
# Objective function to be minimized.
maximize objective: rho_g;