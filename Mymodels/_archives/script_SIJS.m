% SIJS (STI)
clear all
syms b lambda beta p sigma gamma0 rho mu nu S I J

%beta=0.5;p=0.3;sigma=1/12;gamma0=6/7;rho=0.5;mu=1/35;nu=1/5; b=2; 

N = b/mu;
lambda = beta*(I+J)/N;

dS  = b - lambda*S + ((1-p)*nu+rho)*I + (gamma0+nu)*J - mu*S;
dI  = lambda*S - ((1-p)*nu + p*sigma + rho)*I - mu*I;
dJ = p*sigma*I - (gamma0+nu)*J - mu*J;

eqn1 = dS ==0;
eqn2 = dI==0;
eqn3 = dJ==0;

sol = solve([eqn1,eqn2,eqn3], [S,I,J]);

%S
Ssol = (b*mu^2 + b*nu^2 + b*gamma0*mu + b*gamma0*nu + b*gamma0*rho + 2*b*mu*nu + b*mu*rho + b*nu*rho - b*nu^2*p - b*gamma0*nu*p + b*gamma0*p*sigma - b*mu*nu*p + b*mu*p*sigma + b*nu*p*sigma)/(mu*(beta*gamma0 + beta*mu + beta*nu + beta*p*sigma));
Ssol - b*(gamma0 + mu + nu)*(mu + nu + rho - nu*p + p*sigma)/(beta*mu*(gamma0 + mu + nu + p*sigma));
%I
Inum = -((gamma0 + mu + nu)*(b*mu^2 + b*nu^2 - b*beta*gamma0 - b*beta*mu - b*beta*nu + b*gamma0*mu + b*gamma0*nu + b*gamma0*rho + 2*b*mu*nu + b*mu*rho + b*nu*rho - b*nu^2*p - b*beta*p*sigma - b*gamma0*nu*p + b*gamma0*p*sigma - b*mu*nu*p + b*mu*p*sigma + b*nu*p*sigma))/(beta*gamma0^2*mu + 2*beta*gamma0*mu^2 + 2*beta*gamma0*mu*nu + 2*beta*gamma0*mu*p*sigma + beta*mu^3 + 2*beta*mu^2*nu + 2*beta*mu^2*p*sigma + beta*mu*nu^2 + 2*beta*mu*nu*p*sigma + beta*mu*p^2*sigma^2);
%J
Jnum = -(b*mu^2*p*sigma + b*nu^2*p*sigma - b*beta*p^2*sigma^2 + b*gamma0*p^2*sigma^2 + b*mu*p^2*sigma^2 + b*nu*p^2*sigma^2 - b*nu^2*p^2*sigma - b*beta*gamma0*p*sigma - b*beta*mu*p*sigma - b*beta*nu*p*sigma + b*gamma0*mu*p*sigma + b*gamma0*nu*p*sigma + b*gamma0*p*rho*sigma + 2*b*mu*nu*p*sigma + b*mu*p*rho*sigma + b*nu*p*rho*sigma - b*gamma0*nu*p^2*sigma - b*mu*nu*p^2*sigma)/(beta*gamma0^2*mu + 2*beta*gamma0*mu^2 + 2*beta*gamma0*mu*nu + 2*beta*gamma0*mu*p*sigma + beta*mu^3 + 2*beta*mu^2*nu + 2*beta*mu^2*p*sigma + beta*mu*nu^2 + 2*beta*mu*nu*p*sigma + beta*mu*p^2*sigma^2);

Rp = beta*(gamma0+nu+mu+p*sigma)/((gamma0+nu+mu)*(p*sigma+(1-p)*nu+rho+mu));
S = b/(mu*Rp); %ok
I = b/mu*(1-1/Rp)*(gamma0+nu+mu)/(gamma0+nu+p*sigma+mu); %ok
J = b/mu*(1-1/Rp)*p*sigma/(gamma0+nu+p*sigma+mu); %ok

%% utility function
syms c
P = 1-1/Rp;
U = rho*(P-c);
dU = diff(U,rho)
rhohat = solve(dU==0,rho)
R0 = beta*(gamma0+nu+mu+p*sigma)/((gamma0+nu+mu)*(p*sigma+(1-p)*nu+mu));
c1 = 1-1/R0;c2=1/R0-1;
rhohat = (p*sigma+(1-p)*nu+mu)/2*((1-c)*R0-1);%ok
alpha = (p*sigma+(1-p)*nu+mu)*(R0-1); %ok

%% plot U
clear all
p = 0.8; sigma=365/30;nu=1;mu=1/35;
P0 = 0.1;
R0 =2;
%R0 = beta*(gamma0+nu+mu+p*sigma)/((gamma0+nu+mu)*(p*sigma+(1-p)*nu+mu));
%beta = R0/((gamma0+nu+mu+p*sigma)/((gamma0+nu+mu)*(p*sigma+(1-p)*nu+mu)));
vecC = -1:0.1:1;
rhohat = (p*sigma+(1-p)*nu+mu)./2*((1-vecC)*R0-1);%ok
alpha = (p*sigma+(1-p)*nu+mu)*(R0-1); %ok
c0 = 1-1/R0;c1=1/R0-1;
plot_paper1_procedure2(rhohat, c1, c0, alpha, vecC)

