# Voluntary testing rates #
param rho_h >=0;
param rho_s >=0;
param rho_c >=0;
param rho_g >=0;

param rho_hg >=0;
param rho_hc >=0;
param rho_sc >=0;
param rho_sg >=0;
param rho_cg >=0;

param rho_scg >=0;
param rho_hsc >=0;
param rho_hcg >=0;
param rho_hsg >=0;

param rho_hscg >= 0;

param bnd_sup_0 >= 0;

param up_bnd_alpha >=0;
var rho_hs>=0, <= up_bnd_alpha;

# Main parameters and constraints
model setParamVar_global_v7.mod; 


# Equivalent model
s.t. c_mod_1 : sum{k in G_infected} (Y[k]) = 0 ;
s.t. c_mod_2 : sum{k in C_infected} (Y[k]) = 0 ;

# Prevalence of asymptomatics
set infected_hs := H_infectious union S_infected; 
set infected_asympt_hs := infected_hs;
var Prevalence_HS = sum{k in infected_asympt_hs} (Y[k])/N_equ;


# Objective function to be minimized.

#s.t. cobj: Prevalence_HSCG*N >= bnd_sup_0*N;
s.t. cobj1 : Prevalence_H*N_equ >= bnd_sup_0*N_equ;
s.t. cobj2 : Prevalence_S*N_equ >= bnd_sup_0*N_equ;

maximize objective: rho_hs;