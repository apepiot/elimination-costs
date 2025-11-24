# Voluntary testing rates #

param rho_s >=0;
param rho_h >=0;
param rho_g >=0;

param rho_hs >=0;
param rho_hg >=0;
param rho_sc >=0;
param rho_sg >=0;
param rho_cg >=0;
param rho_hc >=0;

param rho_scg >=0;
param rho_hsc >=0;
param rho_hcg >=0;
param rho_hsg >=0;

param rho_hscg >= 0;

param bnd_sup_0 >= 0;

param up_bnd_alpha >=0;
param inf_bnd_alpha >=0;

var rho_c >= inf_bnd_alpha, <= up_bnd_alpha;

# Main parameters and constraints
model setParamVar_global_v7.mod; 

# Equivalent model
s.t. c_mod_2 : sum{k in G_infected} (Y[k]) = 0 ;
s.t. c_mod_3 : sum{k in S_infected} (Y[k]) = 0 ;
s.t. c_mod_4 : sum{k in H_infected} (Y[k]) = 0 ;

s.t. c_mod_5 : Prevalence_C*N_equ*p_c_0 >= bnd_sup_0*N_equ*p_c_0;


# Objective function to be minimized.
minimize Cost: rho_c*(c - Prevalence_C);

