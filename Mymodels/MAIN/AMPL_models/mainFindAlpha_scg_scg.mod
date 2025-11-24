# Voluntary testing rates #

param rho_h >=0;
param rho_s >=0;
param rho_c >=0;
param rho_g >=0;

param rho_hs >=0;
param rho_hc >=0;
param rho_hg >=0;
param rho_sc >=0;
param rho_sg >=0;
param rho_cg >=0;

param rho_hsc >=0;
param rho_hsg >=0;
param rho_hcg >=0;

param rho_hscg >=0;

param bnd_sup_0>=0;

param up_bnd_alpha >=0;
var rho_scg >=0, <= up_bnd_alpha;


# Prevalence of asymptomatics
model setParamVar_global_v7.mod; 
set infected_scg := S_infected union C_infected union G_infected; 
set infected_asympt_scg := infected_scg diff {C_IS union G_IS};
var Prevalence_SCG = sum{k in infected_asympt_scg} (Y[k])/N_equ;

s.t. c_mod1 : sum{k in H_infectious} (Y[k]) = 0;

# Objective function to be minimized.

s.t. cobj2 : Prevalence_S*N_equ >= bnd_sup_0*N_equ;
s.t. cobj3 : Prevalence_C*N_equ >= bnd_sup_0*N_equ;
s.t. cobj4 : Prevalence_G*N_equ >= bnd_sup_0*N_equ;

maximize objective: rho_scg;