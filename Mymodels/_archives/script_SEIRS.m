% SEIRS (syphilis)
clear all
syms b lambda beta sigma0 gamma ksi mu S E I R% rho

N = b/mu;
lambda = (beta*I)/N;
sigma = sigma0;%+rho;

dS = b - lambda*S + ksi*R - mu*S;
dE = lambda*S - sigma*E - mu*E;
dI = sigma0*E  - gamma*I - mu*I;
dR = gamma*I-ksi*R-mu*R;

eqn1 = dS==0;
eqn2 = dE==0;
eqn3 = dI==0;
eqn4 = dR==0;

sol = solve([eqn1,eqn2,eqn3,eqn4], [S,E,I,R])

%R
-(b*gamma^2*mu + b*sigma0*gamma^2 + b*gamma*mu^2 + b*sigma0*gamma*mu - b*beta*sigma0*gamma)/(beta*mu^3 + beta*mu^2*sigma0 + beta*gamma*mu^2 + beta*ksi*mu^2 + beta*gamma*ksi*mu + beta*gamma*mu*sigma0 + beta*ksi*mu*sigma0)
