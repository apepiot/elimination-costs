# Voluntary testing rates #

param rho_h >=0;
param rho_s >=0;
param rho_c >=0;
param rho_g >=0;

param rho_hs >=0;
param rho_hg >=0;
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

var rho_hc>=inf_bnd_alpha, <= up_bnd_alpha;

# Main parameters and constraints
model setParamVar_global_v7.mod; 

# Equivalent model
s.t. c_mod_1 : sum{k in S_infected} (Y[k]) = 0 ;
s.t. c_mod_2 : sum{k in G_infected} (Y[k]) = 0 ;

# Prevalence of asymptomatics
set infected_hc := H_infectious union C_infected; 
set infected_asympt_hc := infected_hc diff C_IS;
var Prevalence_HC = sum{k in infected_asympt_hc} (Y[k])/N_equ;

s.t. cobj1 : Prevalence_H*N_equ*p_h_0 >= bnd_sup_0*N_equ*p_h_0;
s.t. cobj2 : Prevalence_C*N_equ*p_c_0 >= bnd_sup_0*N_equ*p_c_0;

#Prevalence of asymptomatic, untreated of HIV
set asympt_h := H_infectious diff C_IS;
var Prevalence_H_HC = sum{k in asympt_h}(Y[k])/N_equ;
 
# Objective function to be minimized.
minimize Cost: rho_hc*(c - Prevalence_H_HC);

