

#var rho_cg >= 0, <= 15; #up_bnd_alpha;
var rho_g >=0 , <= 15;

model defConstraints.mod;


# Objective function to be minimized.
#minimize Cost: rho_cg*c - rho_cg*Prevalence_C_and_G;
minimize Cost: rho_g*c - rho_g*Prevalence_G;

#data;
#let {i in 1..560} Y[i] := 1;
# option solver knitroampl;
# option knitro_options 'ms_enable=0 feastol=1e-6 maxtime_real = 1';
#option knitro_options 'ms_enable=1 ms_maxsolves=5 feastol=1e-5 maxtime_real=3';
