param betaIh := 13.194799683714759;
param betaCh := 1.406761238295781;
param sigmah := 6.890538634865351;
param thetah := 0.204544709146265;
param zetah := 0.631739291427463;
param eta_h_prep := 4;
param ph := 0.6000;
#param ph := 0.;

param mu := 1/35;
param b := 1;

param rho_h := 0.4998;
param N := b/mu;

var Y{j in 1..7} >= 0;
var Lambdah;

s.t. conslambda: Lambdah = betaIh*(Y[2]+Y[5])/N + betaCh*(Y[3]+Y[6])/N;

#Objective
maximize obj: 1;#Y[1];

s.t. c1: (1-ph)*b - Lambdah*Y[1]-mu*Y[1]=0;
s.t. c2: Lambdah*Y[1] - (sigmah+rho_h+mu)*Y[2]=0;
s.t. c3: sigmah*Y[2] - (thetah+mu)*Y[3]=0;
s.t. c4: ph*b - (1-zetah)*Lambdah*Y[4] - mu*Y[4]=0;
s.t. c5: (1-zetah)*Lambdah*Y[4] - (sigmah+eta_h_prep+mu)*Y[5]=0;
s.t. c6: sigmah*Y[5] - (thetah+eta_h_prep+mu)*Y[6]=0;
s.t. c7: rho_h*Y[2]+thetah*Y[3] + eta_h_prep*Y[5] + (thetah+eta_h_prep)*Y[6] - mu*Y[7]=0;
s.t. c8: sum{i in 1..7} (Y[i]) = b/mu;
#s.t. c9: (Y[2]+Y[3]+Y[5]+Y[6]+Y[7])/N - 0.718529607876459 = 0;

data;

let {i in 1..7} Y[i] := Uniform01() * N;