% SEIS (gono, chlam)
clear all
syms b lambda beta sigma0 gamma mu S E I rho

N = b/mu;
lambda = (beta*I)/N;
sigma = sigma0+rho;

dS = b - lambda*S + gamma*I + rho*E - mu*S;
dE = lambda*S - sigma*E - mu*E;
dI = sigma0*E  - gamma*I - mu*I;

eqn1 = dS==0;
eqn2 = dE==0;
eqn3 = dI==0;


sol = solve([eqn1,eqn2,eqn3], [S,E,I])
