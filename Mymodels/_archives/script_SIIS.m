% SIIS (STI)
clear all
syms b lambda beta p sigma gamma0 rho mu nu S IA IS 

%beta=0.5;p=0.3;sigma=1/12;gamma0=6/7;rho=0.5;mu=1/35;nu=1/5; b=2; 
gamma = gamma0+rho;

N = b/mu;
lambda = beta*(IA+IS)/N;

dS  = b - lambda*S + (1-p)*(rho+nu)*IA + (gamma+nu)*IS - mu*S;
dIA = lambda*S - ((1-p)*(rho+nu) + p*sigma)*IA - mu*IA;
dIS = p*sigma*IA - (gamma+nu)*IS - mu*IS;

eqn1 = dS ==0;
eqn2 = dIA==0;
eqn3 = dIS==0;

sol = solve([eqn1,eqn2,eqn3], [S,IA,IS]);

%sol.S:
(b*mu^2 + b*nu^2 + b*rho^2 - b*p*rho^2 + b*gamma0*mu + b*gamma0*nu + b*gamma0*rho + 2*b*mu*nu + 2*b*mu*rho + 2*b*nu*rho - b*nu^2*p - b*gamma0*nu*p - b*gamma0*p*rho + b*gamma0*p*sigma - b*mu*nu*p - b*mu*p*rho - 2*b*nu*p*rho + b*mu*p*sigma + b*nu*p*sigma + b*p*rho*sigma)/(mu*(beta*gamma0 + beta*mu + beta*nu + beta*rho + beta*p*sigma))
R = beta*(p*sigma + gamma+mu+nu)/((gamma+mu+nu)*(p*sigma+mu+(1-p)*(nu+rho)))
solve(R==1,rho)
alpha = - nu - (beta - gamma0 - 2*mu + gamma0*p + mu*p - p*sigma + (gamma0^2*p^2 - 2*gamma0^2*p + gamma0^2 + 2*gamma0*mu*p^2 - 2*gamma0*mu*p + 2*gamma0*p^2*sigma - 2*gamma0*p*sigma - 2*beta*gamma0*p + 2*beta*gamma0 + mu^2*p^2 + 2*mu*p^2*sigma - 2*beta*mu*p + p^2*sigma^2 - 4*beta*p^2*sigma + 2*beta*p*sigma + beta^2)^(1/2))/(2*(p - 1))
alpha = (gamma0 - beta + 2*mu - gamma0*p - mu*p + p*sigma + (gamma0^2*p^2 - 2*gamma0^2*p + gamma0^2 + 2*gamma0*mu*p^2 - 2*gamma0*mu*p + 2*gamma0*p^2*sigma - 2*gamma0*p*sigma - 2*beta*gamma0*p + 2*beta*gamma0 + mu^2*p^2 + 2*mu*p^2*sigma - 2*beta*mu*p + p^2*sigma^2 - 4*beta*p^2*sigma + 2*beta*p*sigma + beta^2)^(1/2))/(2*(p - 1)) - nu

b/mu/R

%sol.IA:
-((gamma0 + mu + nu + rho)*(b*mu^2 + b*nu^2 + b*rho^2 - b*p*rho^2 - b*beta*gamma0 - b*beta*mu - b*beta*nu - b*beta*rho + b*gamma0*mu + b*gamma0*nu + b*gamma0*rho + 2*b*mu*nu + 2*b*mu*rho + 2*b*nu*rho - b*nu^2*p - b*beta*p*sigma - b*gamma0*nu*p - b*gamma0*p*rho + b*gamma0*p*sigma - b*mu*nu*p - b*mu*p*rho - 2*b*nu*p*rho + b*mu*p*sigma + b*nu*p*sigma + b*p*rho*sigma))/(beta*gamma0^2*mu + 2*beta*gamma0*mu^2 + 2*beta*gamma0*mu*nu + 2*beta*gamma0*mu*p*sigma + 2*beta*gamma0*mu*rho + beta*mu^3 + 2*beta*mu^2*nu + 2*beta*mu^2*p*sigma + 2*beta*mu^2*rho + beta*mu*nu^2 + 2*beta*mu*nu*p*sigma + 2*beta*mu*nu*rho + beta*mu*p^2*sigma^2 + 2*beta*mu*p*rho*sigma + beta*mu*rho^2)

%sol.IS
-(b*mu^2*p*sigma + b*nu^2*p*sigma + b*p*rho^2*sigma - b*beta*p^2*sigma^2 + b*gamma0*p^2*sigma^2 + b*mu*p^2*sigma^2 + b*nu*p^2*sigma^2 - b*nu^2*p^2*sigma + b*p^2*rho*sigma^2 - b*p^2*rho^2*sigma - b*beta*gamma0*p*sigma - b*beta*mu*p*sigma - b*beta*nu*p*sigma - b*beta*p*rho*sigma + b*gamma0*mu*p*sigma + b*gamma0*nu*p*sigma + b*gamma0*p*rho*sigma + 2*b*mu*nu*p*sigma + 2*b*mu*p*rho*sigma + 2*b*nu*p*rho*sigma - b*gamma0*nu*p^2*sigma - b*gamma0*p^2*rho*sigma - b*mu*nu*p^2*sigma - b*mu*p^2*rho*sigma - 2*b*nu*p^2*rho*sigma)/(beta*gamma0^2*mu + 2*beta*gamma0*mu^2 + 2*beta*gamma0*mu*nu + 2*beta*gamma0*mu*p*sigma + 2*beta*gamma0*mu*rho + beta*mu^3 + 2*beta*mu^2*nu + 2*beta*mu^2*p*sigma + 2*beta*mu^2*rho + beta*mu*nu^2 + 2*beta*mu*nu*p*sigma + 2*beta*mu*nu*rho + beta*mu*p^2*sigma^2 + 2*beta*mu*p*rho*sigma + beta*mu*rho^2)



%%
% clear all
% syms b lambda beta p sigma gamma0 rho mu nu S IA IS R
% gamma=gamma0;
% 
% %R = beta*(p*sigma + gamma+mu+nu)/((gamma+mu+nu)*(p*sigma+mu+(1-p)*(nu+rho)))
% 
% S = b/(mu*R);
% N = b/mu;
% lambda = beta*(IA+IS)/N;
% 
% dS  = b - lambda*S + (1-p)*(rho+nu)*IA + (gamma+nu)*IS - mu*S;
% dIA = lambda*S - ((1-p)*(rho+nu) + p*sigma)*IA - mu*IA;
% dIS = p*sigma*IA - (gamma+nu)*IS - mu*IS;
% eqn1 = dS ==0;
% eqn2 = dIA==0;
% eqn3 = dIS==0;
% 
% sol = solve([eqn3,eqn1], [IA,IS]);
% 
% sol.IA
% -(b*(gamma0 + mu + nu)*(gamma0*mu - beta*mu - beta*nu - beta*gamma0 + gamma0*nu + gamma0*rho + 2*mu*nu + mu*rho + nu*rho - nu^2*p + mu^2 + nu^2 - beta*p*sigma - gamma0*nu*p - gamma0*p*rho + gamma0*p*sigma - mu*nu*p - mu*p*rho - nu*p*rho + mu*p*sigma + nu*p*sigma))/(beta*(gamma0*mu + mu*nu + mu^2 + mu*p*sigma)*(gamma0 + mu + nu + p*sigma))
% (b*(R - 1)*(gamma0 + mu + nu))/(beta*gamma0 + beta*mu + beta*nu - R*nu^2 - R*gamma0*nu - R*gamma0*rho - R*mu*nu - R*mu*rho - R*nu*rho + beta*p*sigma + R*nu^2*p + R*mu*p*rho + R*nu*p*rho - R*nu*p*sigma + R*gamma0*nu*p + R*gamma0*p*rho - R*gamma0*p*sigma + R*mu*nu*p)
% 
% sol.IS
% -(b*p*sigma*(gamma0*mu - beta*mu - beta*nu - beta*gamma0 + gamma0*nu + gamma0*rho + 2*mu*nu + mu*rho + nu*rho - nu^2*p + mu^2 + nu^2 - beta*p*sigma - gamma0*nu*p - gamma0*p*rho + gamma0*p*sigma - mu*nu*p - mu*p*rho - nu*p*rho + mu*p*sigma + nu*p*sigma))/(beta*(gamma0*mu + mu*nu + mu^2 + mu*p*sigma)*(gamma0 + mu + nu + p*sigma))
% (b*p*sigma*(R - 1))/(beta*gamma0 + beta*mu + beta*nu - R*nu^2 - R*gamma0*nu - R*gamma0*rho - R*mu*nu - R*mu*rho - R*nu*rho + beta*p*sigma + R*nu^2*p + R*mu*p*rho + R*nu*p*rho - R*nu*p*sigma + R*gamma0*nu*p + R*gamma0*p*rho - R*gamma0*p*sigma + R*mu*nu*p)
% 


%% Utility function
clear all
syms lambda beta p sigma gamma0 rho mu nu c
gamma = gamma0+rho;
R = beta*(p*sigma + gamma+mu+nu)/((gamma+mu+nu)*(p*sigma+mu+(1-p)*(nu+rho)));
P = 1-1/R;
U = rho*(P-c);
dU=diff(U,rho);
solve(dU==0,rho)

%if p=1
sousracine = sigma*(mu + sigma)*(gamma0 + mu + nu + sigma)*(mu - beta + sigma + beta*c);
sol1 = -(gamma0*mu - beta*mu - beta*nu - beta*gamma0 - beta*sigma + gamma0*sigma + mu*nu + 2*mu*sigma + nu*sigma + mu^2 + sigma^2 + beta*c*gamma0 + beta*c*mu + beta*c*nu + beta*c*sigma - sousracine^(1/2) )/(mu - beta + sigma + beta*c)
sol2 = -(gamma0*mu - beta*mu - beta*nu - beta*gamma0 - beta*sigma + gamma0*sigma + mu*nu + 2*mu*sigma + nu*sigma + mu^2 + sigma^2 + beta*c*gamma0 + beta*c*mu + beta*c*nu + beta*c*sigma + sousracine^(1/2))/(mu - beta + sigma + beta*c)
sol1simp = -((gamma0 + mu + nu + sigma)*(mu - beta + sigma + beta*c) - sousracine^(1/2))/(mu - beta + sigma + beta*c)
sol1simp2 = -(gamma0 + mu + nu + sigma) + sousracine^(1/2)/(mu - beta + sigma + beta*c) %ok


simplify(gamma0*mu - beta*mu - beta*nu - beta*gamma0 - beta*sigma + gamma0*sigma + mu*nu + 2*mu*sigma + nu*sigma + mu^2 + sigma^2 + beta*c*gamma0 + beta*c*mu + beta*c*nu + beta*c*sigma)


%%
clear all
syms b lambda beta p sigma gamma0 rho mu nu S IA IS P
gamma=gamma0+rho;
R = beta*(p*sigma + gamma+mu+nu)/((gamma+mu+nu)*(p*sigma+mu+(1-p)*(nu+rho)));
eqn = P==1-1/R;
sol = solve(eqn,beta)

%% chlam
P=0.074;p=11/100;sigma=1/(11/365);gamma0=1/((30-11)/365);nu=1/(497/365); rho=0;mu=1/35;
gamma=gamma0+rho;
beta = -((mu + p*sigma - (nu + rho)*(p - 1))*(gamma0 + mu + nu + rho))/((P - 1)*(gamma0 + mu + nu + rho + p*sigma));
R = beta*(p*sigma + gamma+mu+nu)/((gamma+mu+nu)*(p*sigma+mu+(1-p)*(nu+rho)));
R

%% gono
P=0.07;p=5/100;sigma=1/(5/365);gamma0=1/(5/365);nu=1/(5/12); rho=0;mu=1/35;
gamma=gamma0+rho;
beta = -((mu + p*sigma - (nu + rho)*(p - 1))*(gamma0 + mu + nu + rho))/((P - 1)*(gamma0 + mu + nu + rho + p*sigma));
R = beta*(p*sigma + gamma+mu+nu)/((gamma+mu+nu)*(p*sigma+mu+(1-p)*(nu+rho)));
R


%% utility function
rho = 0:0.01:1;
Rp = beta*(p*sigma + gamma+mu+nu)./((gamma+mu+nu).*(p*sigma+mu+(1-p).*(nu+rho)));
Pp = 1-1./Rp;
c = 0;
U = rho.*(Pp-c);
plot(rho,max(U,0))


%% argmaxU
alpha = - nu - (beta - gamma0 - 2*mu + gamma0*p + mu*p - p*sigma + (gamma0^2*p^2 - 2*gamma0^2*p + gamma0^2 + 2*gamma0*mu*p^2 - 2*gamma0*mu*p + 2*gamma0*p^2*sigma - 2*gamma0*p*sigma - 2*beta*gamma0*p + 2*beta*gamma0 + mu^2*p^2 + 2*mu*p^2*sigma - 2*beta*mu*p + p^2*sigma^2 - 4*beta*p^2*sigma + 2*beta*p*sigma + beta^2)^(1/2))/(2*(p - 1));

[U0,c0] = U_SIIS(p, beta, sigma, gamma0, mu, nu, 0, 0);
[U1,c1] = U_SIIS(p, beta, sigma, gamma0, mu, nu, alpha, 0);

options = optimset('Display','off'); %options for minsearch
i=1; %vecC =(c1-(c2-c1)/2):(c2-c1)/1000:(c2+(c2-c1)/2);
vecC = -1:0.001:1;
vecRhomax1 = []; vecRhomax2 = []; vecRhomax3 = [];
for c=vecC    
    fun = @(rho) -U_SIIS(p, beta, sigma, gamma0, mu, nu, rho, c);
    vecRhomax1(i) = min(max(fminsearch(fun,0,options),0),alpha);
    i=i+1; %c
end

plot(vecC,vecRhomax1)

plot_paper1_procedure2(vecRhomax1, c1, c0, alpha, vecC)


