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
param rho_scg >=0;

param bnd_sup_0>=0;

param up_bnd_alpha >=0;
var rho_hscg >=0, <= up_bnd_alpha;


# Prevalence of asymptomatics
model setParamVar_global_v7.mod; 
set infected_tot := H_infectious union S_infected union C_infected union G_infected; 
set infected_asympt_tot := infected_tot diff {C_IS union G_IS};
var Prevalence_HSCG = sum{k in infected_asympt_tot} (Y[k])/N_equ;


# Objective function to be minimized.

#s.t. cobj: Prevalence_HSCG*N >= bnd_sup_0*N;
s.t. cobj1 : Prevalence_H*N_equ*p_h_0 >= bnd_sup_0*N_equ*p_h_0;
s.t. cobj2 : Prevalence_S*N_equ*p_s_0 >= bnd_sup_0*N_equ*p_s_0;
s.t. cobj3 : Prevalence_C*N_equ*p_c_0 >= bnd_sup_0*N_equ*p_c_0;
s.t. cobj4 : Prevalence_G*N_equ*p_g_0 >= bnd_sup_0*N_equ*p_g_0;

maximize objective: rho_hscg;