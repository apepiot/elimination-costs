model declare_my_var.mod;
model defConstraints_rho_hscg.mod;

# Objective function to be minimized.
minimize Cost: rho_hscg*(c - Prevalence_HSCG);

