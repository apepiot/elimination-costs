
#---------- Parameters ----------#
param p_h_0 >=0, <=1; #infecton H activee : 1
param p_s_0 >=0, <=1; #infecton S activee
param p_c_0 >=0, <=1; #infecton C activee
param p_g_0 >=0, <=1; #infecton G activee

param betaIh >=0;
param betaCh >=0;
param sigmah >=0;
param thetah >=0;
param zetah >=0;
param eta_h_prep >=0;
param ph_mod >=0;
param ph = ph_mod*p_h_0;

param betas >=0;
param sigmas >=0;
param gamma3s >=0;
param taus >=0;
param thetas >=0;

param beta_c >=0;
param gamma_c >=0;
param nu_c >=0;
param eps_c >=0;
param sigma_c >=0;

param beta_g >=0;
param gamma_g >=0;
param nu_g >=0;
param eps_g >=0;
param sigma_g >=0;

param mu >=0;
param b >=0;
param N_equ := b/mu;

param VTunderART >=0;

param eta_s_prep >=0;
param eta_c_prep >=0;
param eta_g_prep >=0;
param eta_s_art >=0;
param eta_c_art >=0;
param eta_g_art >=0;

param c ;

#---------- Populations ----------#

#HIV
set H_S := 1 .. 554 by 7;
set H_P := 4 .. 557 by 7;
set H_I  := 2 .. 555 by 7;
set H_C  := 3 .. 556 by 7;
set H_Ip := 5 .. 558 by 7;
set H_Cp := 6 .. 559 by 7;
set H_acute := H_I union H_Ip;
set H_chronic := H_C union H_Cp;
set H_infectious := H_acute union H_chronic; #and untreated
set H_infected := {1..560} diff {H_S union H_P}; #with treated


#Syphilis
set S_S1 := 1 .. 526 by 35;
set S_S2 := 2 .. 527 by 35;
set S_S3 := 3 .. 528 by 35;
set S_S4 := 4 .. 529 by 35;
set S_S5 := 5 .. 530 by 35;
set S_S6 := 6 .. 531 by 35;
set S_S7 := 7 .. 532 by 35;
set S_susceptible := S_S1 union S_S2 union S_S3 union S_S4 union S_S5 union S_S6 union S_S7; 

set S_E1 := 8 .. 533 by 35;
set S_E2 := 9 .. 534 by 35;
set S_E3 := 10 .. 535 by 35;
set S_E4 := 11 .. 536 by 35;
set S_E5 := 12 .. 537 by 35;
set S_E6 := 13 .. 538 by 35;
set S_E7 := 14 .. 539 by 35;

set S_E := S_E1 union S_E2 union S_E3 union S_E4 union S_E5 union S_E6 union S_E7;

set S_infected := {1..560} diff S_susceptible;

set S_infectious := S_infected diff S_E;

#Ct
set C_IA1 := 71 .. 105;
set C_IA2 := 211 .. 245;
set C_IA3 := 351 .. 385;
set C_IA4 := 491 .. 525;
set C_IA := C_IA1 union C_IA2 union C_IA3 union C_IA4;

set C_IS1 := 106 .. 140;
set C_IS2 := 246 .. 280;
set C_IS3 := 386 .. 420;
set C_IS4 := 526 .. 560;
set C_IS := C_IS1 union C_IS2 union C_IS3 union C_IS4;

set C_infectious := C_IA union C_IS;

set C_infected := {36..140} union {176..280} union {316..420} union {456..560};

#Ng
set G_infectious := {281 .. 560};
set G_infected := {141..560};
set G_IS := {421 .. 560};
#set CG_infected := {36..560};#infected_C union infected_G; #{36..560};

#---------- Variables ----------#

var Y{j in 1..560} >= 0;
var Lambda_h;
var Lambda_s;
var Lambda_c;
var Lambda_g;

# Prevalences #
var Prevalence_H = sum{k in H_infectious} (Y[k])/N_equ;
var Prevalence_S = sum{k in S_infected} (Y[k])/N_equ;
var Prevalence_C = sum{k in C_infected} (Y[k])/N_equ;
var Prevalence_G = sum{k in G_infected} (Y[k])/N_equ;

s.t. clambdah: Lambda_h = (betaIh*(sum{k in H_acute} (Y[k]))/N_equ + betaCh*(sum{i in H_chronic} (Y[i]))/N_equ);
s.t. clambdas: Lambda_s = betas*(sum{k in S_infectious} (Y[k]))/N_equ;
s.t. clambdaX: Lambda_c = beta_c*(sum{k in C_infectious} (Y[k]))/N_equ;
s.t. clambdaY: Lambda_g = beta_g*(sum{k in G_infectious} (Y[k]))/N_equ;


s.t. c_mod_h : (1-p_h_0)*Prevalence_H = 0 ;
s.t. c_mod_s : (1-p_s_0)*Prevalence_S = 0 ;
s.t. c_mod_c : (1-p_c_0)*Prevalence_C = 0 ;
s.t. c_mod_g : (1-p_g_0)*Prevalence_G = 0 ;

#---------- ODE/non linear equations system ----------#

model ODE_SICTPSEIIISSEIIS2_v7.mod;
s.t. c561: sum{i in 1..560} (Y[i]) = N_equ;

