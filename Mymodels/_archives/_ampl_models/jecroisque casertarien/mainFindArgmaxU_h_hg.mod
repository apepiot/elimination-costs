# Voluntary testing rates #

param rho_h >=0;
param rho_s >=0;
param rho_c >=0;
param rho_g >=0;

param rho_hs >=0;
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
param inf_bnd_alpha >=0;

var rho_hg >= inf_bnd_alpha, <= up_bnd_alpha;

# Main parameters and constraints
model setParamVar_global_v7.mod; 

# Equivalent model
s.t. c_mod_1 : sum{k in S_infected} (Y[k]) = 0 ;
s.t. c_mod_2 : sum{k in C_infected} (Y[k]) = 0 ;
s.t. cobj1 : Prevalence_H*N_equ*p_h_0 >= bnd_sup_0*N_equ*p_h_0;
s.t. cobj2 : Prevalence_G*N_equ*p_g_0 >= bnd_sup_0*N_equ*p_g_0;

# Prevalence of asymptomatics
set infected_hg := H_infectious union G_infected; 
set infected_asympt_hg := infected_hg diff G_IS;
var Prevalence_HG = sum{k in infected_asympt_hg} (Y[k])/N_equ;
 
#Prevalence of asymptomatic, untreated of HIV
set asympt_h := H_infectious diff G_IS;
var Prevalence_H_HG = sum{k in asympt_h}(Y[k])/N_equ;

# Objective function to be minimized.
minimize Cost: rho_hg*(c - Prevalence_H_HG);


