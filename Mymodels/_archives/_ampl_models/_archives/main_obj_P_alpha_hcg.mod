
model defConstraints_rho_hcg.mod;

s.t. cobj: Prevalence_HCG*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_hcg;