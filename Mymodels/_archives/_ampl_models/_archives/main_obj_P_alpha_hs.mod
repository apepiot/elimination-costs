
model defConstraints_rho_hs.mod;

s.t. cobj: Prevalence_HS*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_hs;