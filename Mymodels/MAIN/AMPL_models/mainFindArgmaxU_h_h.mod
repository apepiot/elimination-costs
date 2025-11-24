# Voluntary testing rates #

param rho_s >=0;
param rho_c >=0;
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

param f >=0;          #bias parameter for HIV prevalence

param up_bnd_alpha >=0;
param inf_bnd_alpha >=0;

var rho_h>=inf_bnd_alpha, <= up_bnd_alpha;

# Main parameters and constraints
model setParamVar_global_v7.mod; 

# Equivalent model
s.t. c_mod_1 : sum{k in S_infected} (Y[k]) = 0 ;
s.t. c_mod_2 : sum{k in G_infected} (Y[k]) = 0 ;
s.t. c_mod_3 : sum{k in C_infected} (Y[k]) = 0 ;


# Prevalence of asymptomatics
#Prevalence_H already defined
var Prevalence_H_b = f*Prevalence_H;
#var Prevalence_H_b = max(min(f*Prevalence_H,1),0);


s.t. cobj1 : Prevalence_H*N_equ >= bnd_sup_0*N_equ;
#s.t. cobj2 : Prevalence_C*N_equ >= bnd_sup_0*N_equ;

#---- idem mainFindAlpha_hc_hc.mod

#s.t. cmod: Prevalence_HC*N>=bnd_sup_0*N;
 
# Objective function to be minimized.
minimize Cost: rho_h*(c - Prevalence_H_b);

