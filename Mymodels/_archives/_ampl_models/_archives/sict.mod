param betaIh := 8.91;
param betaCh := 1.19;
param sigmah := 0.094;
param thetah := 0.99;

param mu := 1/35;
param b := 1;

param rho_h := 0.4998;
param N := b/mu;

var Y{j in 1..4} >= 0;
var Lambdah;
s.t. conslambda: Lambdah = betaIh*(Y[2])/N + betaCh*(Y[3])/N;
#let Lambdah := 4;

#Objective
maximize obj: 1;

s.t. c1: b - Lambdah*Y[1] - mu*Y[1]=0;
s.t. c2: Lambdah*Y[1] - (sigmah+rho_h+mu)*Y[2] = 0;
s.t. c3: sigmah*Y[2] - (thetah+mu)*Y[3] = 0;
s.t. c4: thetah*Y[3] + rho_h*Y[2] - (sigmah+mu)*Y[4] = 0;
s.t. c5: sum{i in 1..4} (Y[i]) = b/mu;

data;

let {i in 1..4} Y[i] := Uniform01() * N;