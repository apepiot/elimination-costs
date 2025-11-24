
model defConstraints_rho_sg.mod;

s.t. cobj: Prevalence_SG*N >= bnd_sup_0*N;

# Objective function to be minimized.
maximize objective: rho_sg;