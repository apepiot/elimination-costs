% SCIR (HCV model)
clear all
syms b lambda betaI betaC sigma gamma0 theta0 mu S I C R rho

N = b/mu;
lambda = (betaI*I+betaC*C)/N;
theta = theta0 + rho;
gamma = gamma0 + rho;

dS = b - lambda*S-mu*S;
dI = lambda*S - (sigma+gamma)*I - mu*I;
dC = sigma*I - theta*C - mu*C;
dR = theta*C + gamma*I - mu*R;

eqn1 = dS==0;
eqn2 = dI==0;
eqn3 = dC==0;
eqn4 = dR==0;

sol = solve([eqn1,eqn2,eqn3,eqn4], [S,I,C,R])

S = (b*mu^2 + b*gamma*mu + b*gamma*theta + b*mu*sigma + b*mu*theta + b*sigma*theta)/(mu*(betaI*mu + betaC*sigma + betaI*theta))
Rp = (betaI*(theta+mu) + betaC*sigma)/((mu+theta)*(sigma+gamma+mu))
S = b/mu/Rp;
I = -(b*mu^2 - b*betaI*mu + b*gamma*mu - b*betaC*sigma - b*betaI*theta + b*gamma*theta + b*mu*sigma + b*mu*theta + b*sigma*theta)/(betaI*mu^2 + betaC*sigma^2 + betaI*gamma*mu + betaC*gamma*sigma + betaI*gamma*theta + betaC*mu*sigma + betaI*mu*sigma + betaI*mu*theta + betaI*sigma*theta)
I = b/(sigma+mu+gamma)*(1-1/Rp)

C = -(sigma*(b*mu^2 - b*betaI*mu + b*gamma*mu - b*betaC*sigma - b*betaI*theta + b*gamma*theta + b*mu*sigma + b*mu*theta + b*sigma*theta))/(betaI*mu^3 + betaC*mu*sigma^2 + betaC*mu^2*sigma + betaI*mu^2*sigma + betaI*mu*theta^2 + 2*betaI*mu^2*theta + betaC*sigma^2*theta + betaI*sigma*theta^2 + betaI*gamma*mu^2 + betaI*gamma*theta^2 + betaC*gamma*mu*sigma + 2*betaI*gamma*mu*theta + betaC*gamma*sigma*theta + betaC*mu*sigma*theta + 2*betaI*mu*sigma*theta)
C = b*sigma/(theta+mu)/(sigma+gamma+mu)*(1-1/Rp)

R = -(b*gamma^2*mu^2 + 2*b*gamma^2*mu*theta + b*gamma^2*theta^2 + b*gamma*mu^3 + b*gamma*mu^2*sigma + 2*b*gamma*mu^2*theta - b*betaI*gamma*mu^2 + 3*b*gamma*mu*sigma*theta - b*betaC*gamma*mu*sigma + b*gamma*mu*theta^2 - 2*b*betaI*gamma*mu*theta + 2*b*gamma*sigma*theta^2 - b*betaC*gamma*sigma*theta - b*betaI*gamma*theta^2 + b*mu^2*sigma*theta + b*mu*sigma^2*theta + b*mu*sigma*theta^2 - b*betaI*mu*sigma*theta + b*sigma^2*theta^2 - b*betaC*sigma^2*theta - b*betaI*sigma*theta^2)/(betaI*mu^4 + betaC*mu^3*sigma + betaI*mu^3*sigma + 2*betaI*mu^3*theta + betaC*mu^2*sigma^2 + betaI*mu^2*theta^2 + betaI*gamma*mu^3 + betaC*gamma*mu^2*sigma + betaI*gamma*mu*theta^2 + 2*betaI*gamma*mu^2*theta + betaC*mu*sigma^2*theta + betaC*mu^2*sigma*theta + betaI*mu*sigma*theta^2 + 2*betaI*mu^2*sigma*theta + betaC*gamma*mu*sigma*theta)
R = b/mu*(theta*sigma+gamma*(theta+mu))/((theta+mu)*(sigma+gamma+mu))*(1-1/Rp)

%Rp=1 iff rho = :
%betaI/2 - gamma0/2 - mu - sigma/2 - theta0/2 - (betaI^2 - 2*betaI*gamma0 - 2*betaI*sigma + 2*betaI*theta0 + gamma0^2 + 2*gamma0*sigma - 2*gamma0*theta0 + sigma^2 - 2*sigma*theta0 + 4*betaC*sigma + theta0^2)^(1/2)/2
betaI/2 - gamma0/2 - mu - sigma/2 - theta0/2 + (betaI^2 - 2*betaI*gamma0 - 2*betaI*sigma + 2*betaI*theta0 + gamma0^2 + 2*gamma0*sigma - 2*gamma0*theta0 + sigma^2 - 2*sigma*theta0 + 4*betaC*sigma + theta0^2)^(1/2)/2

% ratio tests/infections
lambda = (betaI.*I+betaC*C)*mu/b;
rho*(S+I+C)./(lambda*S);


%% utility function maximization
clear all
syms b lambda betaI betaC sigma gamma0 theta0 mu S I C R rho epsilon
betaI=betaC;
%theta0=gamma0;

theta = theta0 + rho;
gamma = gamma0 + rho;
Rp = (betaI*(theta+mu) + betaC*sigma)/((mu+theta)*(sigma+gamma+mu));

%prevalence
%P = mu*((theta+mu+sigma)/((theta+mu)*(sigma+gamma+mu)))*(1-1/Rp)
P = mu*((theta+mu+sigma)/(betaI*(theta+mu) + betaC*sigma))*(Rp-1)

syms c
U = rho*(P-c);
dU=diff(U,rho);
solve(dU==0,rho)

%if betaI=betaC and theta0 = gamma0
-(gamma0*mu - (betaC*mu*(gamma0 + mu)*(mu + betaC*c))^(1/2) + mu^2 + betaC*c*gamma0 + betaC*c*mu)/(mu + betaC*c)
-(gamma0*mu + (betaC*mu*(gamma0 + mu)*(mu + betaC*c))^(1/2) + mu^2 + betaC*c*gamma0 + betaC*c*mu)/(mu + betaC*c)


%% numerical simulations
clear all;
% Definition des parametres 
[betaI,betaC,gamma,sigma,theta,s2,b,mu,rho] = random_parameters(true, true);
tspan = 0:1:10000;
Y0 = [1,1,1,1];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_SICR(t,Y,b, betaI, betaC, sigma, theta, gamma, mu,0,'frequency'),tspan,Y0, options);
T = Ys(end,:)
Nnum = sum(T) %b/mu

Rp = (betaI*(theta+mu) + betaC*sigma)/((mu+theta)*(sigma+gamma+mu))
S = b/mu/Rp
I= b/(sigma+mu+gamma)*(1-1/Rp)
C = b*sigma/(theta+mu)/(sigma+gamma+mu)*(1-1/Rp)
prevalence = mu*((theta+mu+sigma)/((theta+mu)*(sigma+gamma+mu)))*(1-1/Rp)

%%
betaI = 0.71*7.4; betaC = 0.08*7.4; theta = 1/3.1; sigma = 52/8.2; mu=1/35; gamma=0;
Rp = (betaI*(theta+mu) + betaC*sigma)/((mu+theta)*(sigma+gamma+mu))

%% argmaxU
clear all;close all;
%betaI = 0.71*7.4; betaC = 0.08*7.4; theta0 = 1/3.1; sigma = 52/8.2; mu=1/35; gamma0=0.05;
R     = 1.19;
mu    = 1/35;
sigma = 365/(8.2*7); %8.2 semaines
theta0 = 9.8;
gamma0 = 365/(8.2*7); %proportion de diag au stade precoce : 40
ratioBeta = 9;
betaC = R*(sigma+gamma0+mu)*(theta0+mu)/(ratioBeta*(theta0+mu)+1);
betaI = ratioBeta*betaC;
alpha = betaI/2 - gamma0/2 - mu - sigma/2 - theta0/2 + (betaI^2 - 2*betaI*gamma0 - 2*betaI*sigma + 2*betaI*theta0 +...
    gamma0^2 + 2*gamma0*sigma - 2*gamma0*theta0 + sigma^2 - 2*sigma*theta0 + 4*betaC*sigma + theta0^2)^(1/2)/2;

%alphasyms = betaI/2 - gamma0/2 - mu - sigma/2 - theta0/2 + (betaI^2 - 2*betaI*gamma0 - 2*betaI*sigma + 2*betaI*theta0 + gamma0^2 + 2*gamma0*sigma - 2*gamma0*theta0 + sigma^2 - 2*sigma*theta0 + 4*betaC*sigma + theta0^2)^(1/2)/2
%souslaracine : ok
alpha = (betaI-gamma0-2*mu-sigma-theta0)/2+((betaI-gamma0-sigma+theta0)^2+4*betaC*sigma)^(1/2)/2;

%[U0,c2] = U_SICR(betaI, betaC, sigma, gamma0, theta0, mu, 0, 0, alpha); %c2 = c0
%[U1,c1] = U_SICR(betaI, betaC, sigma, gamma0, theta0, mu, alpha, 0,alpha); %c1=c1
param.betaI = betaI; param.betaC=betaC;param.sigma=sigma;param.gamma=gamma0;param.theta=theta0;param.alpha=alpha; 
[U0,c2] = U_SICR_v3(param, mu, 0, 0, 1);
[U0,c1] = U_SICR_v3(param, mu, alpha, 0, 1);

options = optimset('Display','off'); %options for minsearch
i=1; 
vecC =(c1-(c2-c1)/2):(c2-c1)/100:(c2+(c2-c1)/2);
vecRhomax1 = [];
for c=vecC    
    fun = @(rho) -U_SICR_v3(param, mu, rho, c, 1); %U_SICR(betaI, betaC, sigma, gamma0, theta0, mu, rho, c,alpha);
    vecRhomax1(i) = min(max(fminsearch(fun,0,options),0),alpha);
    i=i+1; %c
end

%plot(vecC,vecRhomax1)
plot_paper1_procedure2(vecRhomax1, c1, c2, alpha, vecC,'')

%%
%ratio test/infections
b=0.5;
theta = theta0 + vecRhomax1;
gamma = gamma0 + vecRhomax1;
Rp = (betaI.*(theta+mu) + betaC*sigma)./((mu+theta).*(sigma+gamma+mu));
S = b./(mu.*Rp);
I = b./(sigma+mu+gamma).*(1-1./Rp);
C = b*sigma./((theta+mu).*(sigma+gamma+mu)).*(1-1./Rp);
lambda = (betaI.*I+betaC*C)*mu/b;
ratio = vecRhomax1.*(S+I+C)./(lambda.*S).*(vecC>=c1);
plot_paper1_procedure2(log(ratio), c1, c2, max(log(ratio)), vecC)

figure(2)
plot(vecC.*(vecC>=c1),log(ratio))

figure(4)
plot(vecC,vecRhomax1.*(S+I+C))
plot_paper1_procedure2(vecRhomax1.*(S+I+C), c1, c2, max(vecRhomax1.*(S+I+C)), vecC)

figure(5)
plot_paper1_procedure2((S+I+C)./(I+C), c1, c2, max((S+I+C)./(I+C)), vecC)


%xlim(c1,vecC(end))

%%
