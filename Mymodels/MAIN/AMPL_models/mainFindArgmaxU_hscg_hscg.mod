# Voluntary testing rates #

param rho_h >=0;
param rho_s >=0;
param rho_c >=0;
param rho_g >=0;

param rho_hc >=0;
param rho_hg >=0;
param rho_sc >=0;
param rho_sg >=0;
param rho_cg >=0;
param rho_hs >=0;

param rho_scg >=0;
param rho_hsg >=0;
param rho_hsc >=0;
param rho_hcg >=0;

param bnd_sup_0 >= 0;

param up_bnd_alpha >=0;
param inf_bnd_alpha >=0;

param f>=0;

var rho_hscg >= inf_bnd_alpha, <= up_bnd_alpha;

# Main parameters and constraints
model setParamVar_global_v7.mod; 

s.t. c_mod_4 : Prevalence_H*N_equ*p_h_0 >= bnd_sup_0*N_equ*p_h_0;
s.t. c_mod_5 : Prevalence_S*N_equ*p_s_0 >= bnd_sup_0*N_equ*p_s_0;
s.t. c_mod_6 : Prevalence_C*N_equ*p_c_0 >= bnd_sup_0*N_equ*p_c_0;
s.t. c_mod_7 : Prevalence_G*N_equ*p_g_0 >= bnd_sup_0*N_equ*p_g_0;

# Prevalence of asymptomatics
set infected_hscg := H_infectious union S_infected union C_infected union G_infected; 
set infected_asympt_hscg := infected_hscg diff G_IS diff C_IS;
var Prevalence_HSCG = sum{k in infected_asympt_hscg} (Y[k])/N_equ;

var Prevalence_H_HSCG  = sum{k in (infected_asympt_hscg diff G_infectious diff C_infectious diff S_infected)} (Y[k])/N_equ;
#var Prevalence_HSCG_b  = max(min(f*Prevalence_H_HSCG + (Prevalence_HSCG-Prevalence_H_HSCG),1),0);
var Prevalence_HSCG_b  = f*Prevalence_H_HSCG + (Prevalence_HSCG-Prevalence_H_HSCG);


# Objective function to be minimized.
minimize Cost: rho_hscg*(c - Prevalence_HSCG_b);

