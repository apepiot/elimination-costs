model setParamVar_global_P_v7.mod;

set infected_tot := H_infectious union S_infected union C_infected union G_infected; 
set infected_asympt_tot := infected_tot diff {C_IS union G_IS};
var Prevalence_HSCG = sum{k in infected_asympt_tot} (Y[k])/N_equ;


# Objective function 
#minimize objective: 1;
maximize OBJ: Prevalence_H + Prevalence_S + Prevalence_C + Prevalence_G;