% SIJS_v2 (STI)
%diff from SIJS : we add a rho rate from J to S
clear all
syms b lambda beta p sigma gamma0 rho mu nu S I J

beta=0.5;p=0.3;sigma=1/12;gamma0=6/7;rho=5;mu=1/35;nu=1/5; b=2; 
N = b/mu;
lambda = beta*(I+J)/N;

dS  = b - lambda*S + ((1-p)*nu+rho)*I + (gamma0+nu+rho)*J - mu*S;
dI  = lambda*S - ((1-p)*nu + p*sigma + rho)*I - mu*I;
dJ = p*sigma*I - (gamma0+nu+rho)*J - mu*J;

eqn1 = dS ==0;
eqn2 = dI==0;
eqn3 = dJ==0;

sol = solve([eqn1,eqn2,eqn3], [S,I,J]);

Rp = beta*(gamma0+nu+mu+p*sigma+rho)/((gamma0+nu+mu+rho)*(p*sigma+(1-p)*nu+rho+mu));

%S
Ssol = (b*mu^2 + b*nu^2 + b*rho^2 + b*gamma0*mu + b*gamma0*nu + b*gamma0*rho + 2*b*mu*nu + 2*b*mu*rho + 2*b*nu*rho - b*nu^2*p - b*gamma0*nu*p + b*gamma0*p*sigma - b*mu*nu*p - b*nu*p*rho + b*mu*p*sigma + b*nu*p*sigma + b*p*rho*sigma)/(mu*(beta*gamma0 + beta*mu + beta*nu + beta*rho + beta*p*sigma))
%I
Inum = -((gamma0 + mu + nu + rho)*(b*mu^2 + b*nu^2 + b*rho^2 - b*beta*gamma0 - b*beta*mu - b*beta*nu - b*beta*rho + b*gamma0*mu + b*gamma0*nu + b*gamma0*rho + 2*b*mu*nu + 2*b*mu*rho + 2*b*nu*rho - b*nu^2*p - b*beta*p*sigma - b*gamma0*nu*p + b*gamma0*p*sigma - b*mu*nu*p - b*nu*p*rho + b*mu*p*sigma + b*nu*p*sigma + b*p*rho*sigma))/(beta*gamma0^2*mu + 2*beta*gamma0*mu^2 + 2*beta*gamma0*mu*nu + 2*beta*gamma0*mu*p*sigma + 2*beta*gamma0*mu*rho + beta*mu^3 + 2*beta*mu^2*nu + 2*beta*mu^2*p*sigma + 2*beta*mu^2*rho + beta*mu*nu^2 + 2*beta*mu*nu*p*sigma + 2*beta*mu*nu*rho + beta*mu*p^2*sigma^2 + 2*beta*mu*p*rho*sigma + beta*mu*rho^2)
%J
Jnum = -(b*mu^2*p*sigma + b*nu^2*p*sigma + b*p*rho^2*sigma - b*beta*p^2*sigma^2 + b*gamma0*p^2*sigma^2 + b*mu*p^2*sigma^2 + b*nu*p^2*sigma^2 - b*nu^2*p^2*sigma + b*p^2*rho*sigma^2 - b*beta*gamma0*p*sigma - b*beta*mu*p*sigma - b*beta*nu*p*sigma - b*beta*p*rho*sigma + b*gamma0*mu*p*sigma + b*gamma0*nu*p*sigma + b*gamma0*p*rho*sigma + 2*b*mu*nu*p*sigma + 2*b*mu*p*rho*sigma + 2*b*nu*p*rho*sigma - b*gamma0*nu*p^2*sigma - b*mu*nu*p^2*sigma - b*nu*p^2*rho*sigma)/(beta*gamma0^2*mu + 2*beta*gamma0*mu^2 + 2*beta*gamma0*mu*nu + 2*beta*gamma0*mu*p*sigma + 2*beta*gamma0*mu*rho + beta*mu^3 + 2*beta*mu^2*nu + 2*beta*mu^2*p*sigma + 2*beta*mu^2*rho + beta*mu*nu^2 + 2*beta*mu*nu*p*sigma + 2*beta*mu*nu*rho + beta*mu*p^2*sigma^2 + 2*beta*mu*p*rho*sigma + beta*mu*rho^2)

S = b/(mu*Rp); %ok
I = b/mu*(1-1/Rp)*(gamma0+rho+nu+mu)/(gamma0+rho+nu+p*sigma+mu); %ok
J = b/mu*(1-1/Rp)*p*sigma/(gamma0+nu+rho+p*sigma+mu); %ok

%% utility function
syms c
P = 1-1/Rp;
U = rho*(P-c);
dU = diff(U,rho)
%rhohat = solve(dU==0,rho) %no solution

R0 = beta*(gamma0+nu+mu+p*sigma)/((gamma0+nu+mu)*(p*sigma+(1-p)*nu+mu));
%reprendre ici

c0 = 1-1/R0;%P(0)
sousracine = gamma0^2 + 2*gamma0*nu*p - 2*gamma0*p*sigma + 2*beta*gamma0 + nu^2*p^2 - 2*nu*p^2*sigma + 2*beta*nu*p + p^2*sigma^2 + 2*beta*p*sigma + beta^2;
sousracine2 = (beta+p*sigma+gamma0+nu*p)^2 - 4*p*sigma*(gamma0+nu*p);
%sousracine2 = (beta-(gamma0+(1-p)*nu+nu+p*sigma))^2;
alpha1 = beta/2 - gamma0/2 - mu - nu + (sousracine)^(1/2)/2 + (nu*p)/2 - (p*sigma)/2
alpha2 = (beta-(gamma0+2*mu+(1-p)*nu+nu+p*sigma))/2 + sqrt((gamma0+nu*p+beta+p*sigma)^2-4*p*sigma*(gamma0+nu*p))/2

%% plot U
clear all;
beta=0.5;p=0.3;sigma=1/12;gamma0=6/7;rho=5;mu=1/35;nu=1/5; b=2; 

alpha = (beta-(gamma0+2*mu+(1-p)*nu+nu+p*sigma))/2 + sqrt((gamma0+nu*p+beta+p*sigma)^2-4*p*sigma*(gamma0+nu*p))/2;

%[U0,c0] = U_SIJSv2(p, beta, sigma, gamma0, mu, nu, 0, 0);
R0 = beta*(gamma0+nu+mu+p*sigma)/((gamma0+nu+mu)*(p*sigma+(1-p)*nu+mu));
c0=1-1/R0;
[U1,c1] = U_SIJSv2(p, beta, sigma, gamma0, mu, nu, alpha, 0);

options = optimset('Display','off'); %options for minsearch
i=1; %vecC =(c1-(c2-c1)/2):(c2-c1)/1000:(c2+(c2-c1)/2);
vecC = -1.22:0.001:1.5;
vecRhomax1 = []; vecRhomax2 = []; vecRhomax3 = [];
for c=vecC    
    fun = @(rho) -U_SIJSv2(p, beta, sigma, gamma0, mu, nu, rho, c);
    vecRhomax1(i) = min(max(fminsearch(fun,0,options),0),alpha);
    i=i+1; %c
end

plot(vecC,vecRhomax1)

plot_paper1_procedure2(vecRhomax1, c1, c0, alpha, vecC)

