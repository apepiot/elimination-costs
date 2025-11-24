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

param rho_hscg >= 0;

param bnd_sup_0 >= 0;

param up_bnd_alpha >=0;
param inf_bnd_alpha >=0;

param f>=0;

var rho_hcg >= inf_bnd_alpha, <= up_bnd_alpha;

# Main parameters and constraints
model setParamVar_global_v7.mod; 

# Equivalent model
s.t. c_mod_2 : sum{k in S_infected} (Y[k]) = 0 ;
s.t. cobj1 : Prevalence_H*N_equ*p_h_0 >= bnd_sup_0*N_equ*p_h_0;
s.t. cobj2 : Prevalence_C*N_equ*p_c_0 >= bnd_sup_0*N_equ*p_c_0;
s.t. cobj3 : Prevalence_G*N_equ*p_g_0 >= bnd_sup_0*N_equ*p_g_0;


# Prevalence of asymptomatics
set infected_hcg := H_infectious union C_infected union G_infected; 
set infected_asympt_hcg := infected_hcg diff G_IS diff C_IS;
var Prevalence_HCG = sum{k in infected_asympt_hcg} (Y[k])/N_equ;
var Prevalence_H_HCG  = sum{k in (infected_asympt_hcg diff G_infectious diff C_infectious)} (Y[k])/N_equ;
#var Prevalence_HCG_b  = max(min(f*Prevalence_H_HCG + (Prevalence_HCG-Prevalence_H_HCG),1),0);
var Prevalence_HCG_b  = f*Prevalence_H_HCG + (Prevalence_HCG-Prevalence_H_HCG);
# Prevalence_HCG_b >= 0, <= 1;
#s.t. Prevalence_HCG_b_con: Prevalence_HCG_b <= f*Prevalence_H_HCG + (Prevalence_HCG-Prevalence_H_HCG);
#f*Prevalence_H_HCG + (Prevalence_HCG-Prevalence_H_HCG); 

# Objective function to be minimized.
minimize Cost: rho_hcg*(c - Prevalence_HCG_b);

