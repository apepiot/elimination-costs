model defConstraints_sictp_seiis_U.mod;

#s.t. c30: Prevalence_H>*N=bnd_sup_0*N;
#s.t. c31: Prevalence_X*N>=bnd_sup_0*N;

s.t. c32: Prevalence_HX*N>=bnd_sup_0*N;
 
# Objective function to be minimized.
minimize Cost: rho_hX*(c - Prevalence_HX);

