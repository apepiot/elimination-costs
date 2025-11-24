# Define data
param betaIh >=0;
param betaCh >=0;
param sigmah >=0;
param thetah >=0;
param zetah >=0;
param eta_h_prep >=0;
param ph >=0;

param betaX >=0;
param gammaX >=0;
param nuX >=0;
param epsX >=0;
param sigmaX >=0;


param mu >=0;
param b >=0;
param up_bnd_alpha >=0;

param rho_hX >=0;
param rho_h >=0;
param rho_X >=0;
param eta_X_prep >=0;
param eta_X_art >=0;
param VTunderART >=0;
param bnd_sup_0 >=0;# 0.001*(b/mu);

param rho_Xg=0;

#--------------------------------------------------------------------------------------

set acute_H := {2,5,9,12,16,19,23,26};
set chronic_H := {3,6,10,13,17,20,24,27};
set infected_H = acute_H union chronic_H;

set infectious_X := {15..28};
set infected_X := {8..28};


#--------------------------------------------------------------------------------------

param N := b/mu;

# Define variables and enforce that they be non-negative.
var Y{j in 1..28} >= 0;
var Lambdah;
var LambdaX;

var Prevalence_H = (sum{k in infected_H}(Y[k]))/N;
var Prevalence_X = (sum{k in infected_X}(Y[k]))/N;

var Prevalence_HX = (sum{k in {2,3,5,6}union {8..21}}(Y[k]))/N; #enlever les gens sous prep

s.t. clambdah: Lambdah = (betaIh*(sum{k in acute_H} (Y[k]))/N + betaCh*(sum{i in chronic_H} (Y[i]))/N);
s.t. clambdaX: LambdaX = betaX*(sum{k in infectious_X} (Y[k]))/N;

model ODE_SICTPSEIIS_v7.mod;

s.t. c29: sum{i in 1..28} (Y[i]) = b/mu;
