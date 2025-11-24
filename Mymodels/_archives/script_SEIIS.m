% SEIIS (STI) (voir script_SEIIS_v2)
clear all
syms b lambda beta p sigma gamma0 rho mu nu S E IS IA

%beta=0.5;p=0.3;sigma=1/12;gamma0=6/7;rho=0.5;mu=1/35;nu=1/5; b=2; 

N = b/mu;
lambda = beta*(IS+IA)/N;

dS  = b - lambda*S + (gamma0+nu)*IS + (rho+nu)*IA - mu*S;
dE  = lambda*S - sigma*E - mu*E;
dIS = p*sigma*E - (gamma0+nu)*IS - mu*IS;
dIA = (1-p)*sigma*E - (rho+nu)*IA - mu*IA;

eqn1 = dS ==0;
eqn2 = dE==0;
eqn3 = dIS==0;
eqn4 = dIA==0;
sol = solve([eqn1,eqn2,eqn3,eqn4], [S,E,IS,IA]);

%S 
(b*mu^3 + b*mu^2*rho + b*mu^2*sigma + b*nu^2*sigma + b*gamma0*mu^2 + b*mu*nu^2 + 2*b*mu^2*nu + b*gamma0*mu*nu + b*gamma0*mu*rho + b*gamma0*mu*sigma + b*gamma0*nu*sigma + b*gamma0*rho*sigma + b*mu*nu*rho + 2*b*mu*nu*sigma + b*mu*rho*sigma + b*nu*rho*sigma)/(mu*sigma*(beta*gamma0 + beta*mu + beta*nu - beta*gamma0*p + beta*p*rho))
b*(mu + sigma)*(gamma0 + mu + nu)*(mu + nu + rho)/(beta*mu*sigma*(gamma0 + mu + nu - gamma0*p + p*rho))
%E
-((mu + nu + rho)*(b*mu^4 + b*mu^3*rho + b*mu^3*sigma + b*nu^3*sigma + b*gamma0^2*mu^2 + 3*b*mu^2*nu^2 + 2*b*gamma0*mu^3 + b*mu*nu^3 + 3*b*mu^3*nu - b*beta*gamma0^2*sigma - b*beta*mu^2*sigma + 2*b*gamma0*mu*nu^2 + 4*b*gamma0*mu^2*nu + b*gamma0^2*mu*nu - b*beta*nu^2*sigma + 2*b*gamma0*mu^2*rho + b*gamma0^2*mu*rho + 2*b*gamma0*mu^2*sigma + b*gamma0^2*mu*sigma + 2*b*gamma0*nu^2*sigma + b*gamma0^2*nu*sigma + b*gamma0^2*rho*sigma + b*mu*nu^2*rho + 2*b*mu^2*nu*rho + 3*b*mu*nu^2*sigma + 3*b*mu^2*nu*sigma + b*mu^2*rho*sigma + b*nu^2*rho*sigma - 2*b*beta*gamma0*mu*sigma - 2*b*beta*gamma0*nu*sigma - 2*b*beta*mu*nu*sigma + 2*b*gamma0*mu*nu*rho + 4*b*gamma0*mu*nu*sigma + 2*b*gamma0*mu*rho*sigma + 2*b*gamma0*nu*rho*sigma + 2*b*mu*nu*rho*sigma + b*beta*gamma0^2*p*sigma + b*beta*gamma0*mu*p*sigma + b*beta*gamma0*nu*p*sigma - b*beta*gamma0*p*rho*sigma - b*beta*mu*p*rho*sigma - b*beta*nu*p*rho*sigma))/(sigma*(- beta*gamma0^2*mu^2*p + beta*gamma0^2*mu^2 - beta*gamma0^2*mu*nu*p + beta*gamma0^2*mu*nu + sigma*beta*gamma0^2*mu*p^2 - beta*gamma0^2*mu*p*rho - 2*sigma*beta*gamma0^2*mu*p + beta*gamma0^2*mu*rho + sigma*beta*gamma0^2*mu - beta*gamma0*mu^3*p + 2*beta*gamma0*mu^3 - 2*beta*gamma0*mu^2*nu*p + 4*beta*gamma0*mu^2*nu - 2*sigma*beta*gamma0*mu^2*p + 2*beta*gamma0*mu^2*rho + 2*sigma*beta*gamma0*mu^2 - beta*gamma0*mu*nu^2*p + 2*beta*gamma0*mu*nu^2 - 2*sigma*beta*gamma0*mu*nu*p + 2*beta*gamma0*mu*nu*rho + 2*sigma*beta*gamma0*mu*nu - 2*sigma*beta*gamma0*mu*p^2*rho + beta*gamma0*mu*p*rho^2 + 2*sigma*beta*gamma0*mu*p*rho + beta*mu^4 + 3*beta*mu^3*nu + beta*mu^3*p*rho + beta*mu^3*rho + sigma*beta*mu^3 + 3*beta*mu^2*nu^2 + 2*beta*mu^2*nu*p*rho + 2*beta*mu^2*nu*rho + 2*sigma*beta*mu^2*nu + beta*mu^2*p*rho^2 + 2*sigma*beta*mu^2*p*rho + beta*mu*nu^3 + beta*mu*nu^2*p*rho + beta*mu*nu^2*rho + sigma*beta*mu*nu^2 + beta*mu*nu*p*rho^2 + 2*sigma*beta*mu*nu*p*rho + sigma*beta*mu*p^2*rho^2))
%IS
-(p*(mu + nu + rho)*(b*mu^3 + b*mu^2*rho + b*mu^2*sigma + b*nu^2*sigma + b*gamma0*mu^2 + b*mu*nu^2 + 2*b*mu^2*nu - b*beta*mu*sigma + b*gamma0*mu*nu - b*beta*nu*sigma + b*gamma0*mu*rho + b*gamma0*mu*sigma + b*gamma0*nu*sigma + b*gamma0*rho*sigma + b*mu*nu*rho + 2*b*mu*nu*sigma + b*mu*rho*sigma + b*nu*rho*sigma - b*beta*gamma0*sigma + b*beta*gamma0*p*sigma - b*beta*p*rho*sigma))/(- beta*gamma0^2*mu^2*p + beta*gamma0^2*mu^2 - beta*gamma0^2*mu*nu*p + beta*gamma0^2*mu*nu + sigma*beta*gamma0^2*mu*p^2 - beta*gamma0^2*mu*p*rho - 2*sigma*beta*gamma0^2*mu*p + beta*gamma0^2*mu*rho + sigma*beta*gamma0^2*mu - beta*gamma0*mu^3*p + 2*beta*gamma0*mu^3 - 2*beta*gamma0*mu^2*nu*p + 4*beta*gamma0*mu^2*nu - 2*sigma*beta*gamma0*mu^2*p + 2*beta*gamma0*mu^2*rho + 2*sigma*beta*gamma0*mu^2 - beta*gamma0*mu*nu^2*p + 2*beta*gamma0*mu*nu^2 - 2*sigma*beta*gamma0*mu*nu*p + 2*beta*gamma0*mu*nu*rho + 2*sigma*beta*gamma0*mu*nu - 2*sigma*beta*gamma0*mu*p^2*rho + beta*gamma0*mu*p*rho^2 + 2*sigma*beta*gamma0*mu*p*rho + beta*mu^4 + 3*beta*mu^3*nu + beta*mu^3*p*rho + beta*mu^3*rho + sigma*beta*mu^3 + 3*beta*mu^2*nu^2 + 2*beta*mu^2*nu*p*rho + 2*beta*mu^2*nu*rho + 2*sigma*beta*mu^2*nu + beta*mu^2*p*rho^2 + 2*sigma*beta*mu^2*p*rho + beta*mu*nu^3 + beta*mu*nu^2*p*rho + beta*mu*nu^2*rho + sigma*beta*mu*nu^2 + beta*mu*nu*p*rho^2 + 2*sigma*beta*mu*nu*p*rho + sigma*beta*mu*p^2*rho^2)
%IA
(b*mu^4*p - b*mu^3*rho - b*mu^3*sigma - b*nu^3*sigma - b*gamma0^2*mu^2 - 3*b*mu^2*nu^2 - 2*b*gamma0*mu^3 - b*mu*nu^3 - 3*b*mu^3*nu - b*mu^4 + b*beta*gamma0^2*sigma + b*beta*mu^2*sigma - 2*b*gamma0*mu*nu^2 - 4*b*gamma0*mu^2*nu - b*gamma0^2*mu*nu + b*beta*nu^2*sigma + 2*b*gamma0*mu^3*p - 2*b*gamma0*mu^2*rho - b*gamma0^2*mu*rho - 2*b*gamma0*mu^2*sigma - b*gamma0^2*mu*sigma - 2*b*gamma0*nu^2*sigma - b*gamma0^2*nu*sigma + b*mu*nu^3*p + 3*b*mu^3*nu*p - b*gamma0^2*rho*sigma - b*mu*nu^2*rho - 2*b*mu^2*nu*rho - 3*b*mu*nu^2*sigma - 3*b*mu^2*nu*sigma + b*mu^3*p*rho + b*mu^3*p*sigma + b*nu^3*p*sigma - b*mu^2*rho*sigma - b*nu^2*rho*sigma + b*gamma0^2*mu^2*p + 3*b*mu^2*nu^2*p + b*beta*gamma0^2*p^2*sigma + 2*b*beta*gamma0*mu*sigma + 2*b*beta*gamma0*nu*sigma + 2*b*beta*mu*nu*sigma - 2*b*gamma0*mu*nu*rho - 4*b*gamma0*mu*nu*sigma - 2*b*gamma0*mu*rho*sigma - 2*b*gamma0*nu*rho*sigma - 2*b*mu*nu*rho*sigma - 2*b*beta*gamma0^2*p*sigma - b*beta*mu^2*p*sigma + 2*b*gamma0*mu*nu^2*p + 4*b*gamma0*mu^2*nu*p + b*gamma0^2*mu*nu*p - b*beta*nu^2*p*sigma + 2*b*gamma0*mu^2*p*rho + b*gamma0^2*mu*p*rho + 2*b*gamma0*mu^2*p*sigma + b*gamma0^2*mu*p*sigma + 2*b*gamma0*nu^2*p*sigma + b*gamma0^2*nu*p*sigma + b*gamma0^2*p*rho*sigma + b*mu*nu^2*p*rho + 2*b*mu^2*nu*p*rho + 3*b*mu*nu^2*p*sigma + 3*b*mu^2*nu*p*sigma + b*mu^2*p*rho*sigma + b*nu^2*p*rho*sigma + b*beta*gamma0*mu*p^2*sigma + b*beta*gamma0*nu*p^2*sigma - b*beta*gamma0*p^2*rho*sigma - b*beta*mu*p^2*rho*sigma - b*beta*nu*p^2*rho*sigma - 3*b*beta*gamma0*mu*p*sigma - 3*b*beta*gamma0*nu*p*sigma + b*beta*gamma0*p*rho*sigma - 2*b*beta*mu*nu*p*sigma + b*beta*mu*p*rho*sigma + 2*b*gamma0*mu*nu*p*rho + b*beta*nu*p*rho*sigma + 4*b*gamma0*mu*nu*p*sigma + 2*b*gamma0*mu*p*rho*sigma + 2*b*gamma0*nu*p*rho*sigma + 2*b*mu*nu*p*rho*sigma)/(- beta*gamma0^2*mu^2*p + beta*gamma0^2*mu^2 - beta*gamma0^2*mu*nu*p + beta*gamma0^2*mu*nu + sigma*beta*gamma0^2*mu*p^2 - beta*gamma0^2*mu*p*rho - 2*sigma*beta*gamma0^2*mu*p + beta*gamma0^2*mu*rho + sigma*beta*gamma0^2*mu - beta*gamma0*mu^3*p + 2*beta*gamma0*mu^3 - 2*beta*gamma0*mu^2*nu*p + 4*beta*gamma0*mu^2*nu - 2*sigma*beta*gamma0*mu^2*p + 2*beta*gamma0*mu^2*rho + 2*sigma*beta*gamma0*mu^2 - beta*gamma0*mu*nu^2*p + 2*beta*gamma0*mu*nu^2 - 2*sigma*beta*gamma0*mu*nu*p + 2*beta*gamma0*mu*nu*rho + 2*sigma*beta*gamma0*mu*nu - 2*sigma*beta*gamma0*mu*p^2*rho + beta*gamma0*mu*p*rho^2 + 2*sigma*beta*gamma0*mu*p*rho + beta*mu^4 + 3*beta*mu^3*nu + beta*mu^3*p*rho + beta*mu^3*rho + sigma*beta*mu^3 + 3*beta*mu^2*nu^2 + 2*beta*mu^2*nu*p*rho + 2*beta*mu^2*nu*rho + 2*sigma*beta*mu^2*nu + beta*mu^2*p*rho^2 + 2*sigma*beta*mu^2*p*rho + beta*mu*nu^3 + beta*mu*nu^2*p*rho + beta*mu*nu^2*rho + sigma*beta*mu*nu^2 + beta*mu*nu*p*rho^2 + 2*sigma*beta*mu*nu*p*rho + sigma*beta*mu*p^2*rho^2)

sigma=1/10;beta=0.5;gamma0=3/10; mu=1/35;nu=1;p=0.7;rho=1;b=3;
lambdamat = beta*simplify(sol.IA+sol.IS)*mu/b;
lambdamat = -(gamma0*mu^2 + mu*nu^2 + 2*mu^2*nu + mu^2*rho + mu^2*sigma + nu^2*sigma + mu^3 - beta*gamma0*sigma - beta*mu*sigma + gamma0*mu*nu - beta*nu*sigma + gamma0*mu*rho + gamma0*mu*sigma + gamma0*nu*sigma + gamma0*rho*sigma + mu*nu*rho + 2*mu*nu*sigma + mu*rho*sigma + nu*rho*sigma + beta*gamma0*p*sigma - beta*p*rho*sigma)/(gamma0*mu + gamma0*nu + gamma0*rho + gamma0*sigma + 2*mu*nu + mu*rho + nu*rho + mu*sigma + nu*sigma + mu^2 + nu^2 - gamma0*p*sigma + p*rho*sigma)

% 
Rp = (beta*sigma*(gamma0 + mu + nu - gamma0*p + p*rho))/((mu + sigma)*(gamma0 + mu + nu)*(mu + nu + rho));
S  = b/(mu*Rp);
E = lambdamat*S/(sigma+mu);
lambda = (Rp*(sigma+mu)*beta*(1-1/Rp))/(Rp*(sigma+mu)+beta);

syms c
P=1-1./Rp
U=rho*(P-c)

dU=diff(U,rho)



