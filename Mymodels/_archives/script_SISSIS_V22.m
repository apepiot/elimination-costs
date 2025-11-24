%% Script mymodel SISxSIS (V2.2) %KIT TESTING
clear all
close all
%% Parametres
% Populations initiales
N0 = 150;
S0 = 95;
I10 = 2;
I20 = 3;
I120 = N0-S0-I10-I20;

% Definition des parametres 
beta1 = 0.8;
beta2 = 0.6/4;
gamma1 = 0.02;
gamma2 = 0.02;
b = 0.05;
mu = 0.01;
s1 = 0.9;
s2 = 0.8;
rho = 0.2;

%[beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');


% Parametres du systeme d'ODE 
tspan = 0:1:5000;
Y0 = [S0; I10; I20; I120];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) SISxSIS_V22(t,Y,b,beta1, beta2,gamma1,gamma2, s1,s2, rho, mu, 'frequency'),tspan,Y0, options);

Ys(end,:)
sum(Ys(end,:))

% Plot
% figure(1)
% plot(ts,Ys)%./(sum(Ys')'));
% title('Evolution');
% xlabel('Time t');
% ylabel('Y');

%% Endemic equilibrium
gamma2p = gamma2 + s2*rho;
gamma1p = gamma1 + s1*rho;
R1p = beta1/(gamma1p + mu);
R2p = beta2/(gamma2p + mu);


Nequ = b/mu;
E0 = [b/mu,0,0,0];
E1 = [Nequ/R1p, Nequ - Nequ/R1p,0,0]
E2 = [Nequ/R2p, 0, Nequ - Nequ/R2p,0]

% num = gamma2p/R1p + gamma1p/R2p + mu + (s1*s2*rho)*(1-1/R1p-1/R2p)
% denom = (beta1+beta2-mu-s1*s2*rho)
% S12 = b/mu*num/denom -Ys(end,1)
% S12 = b/mu*((gamma1+s1*rho)/R2p + (gamma2+s2*rho)/R1p + mu + s1*s2*rho*I12*mu/b)/(beta1+beta2-mu)%OK
% 
% I12 = S12 + Nequ*(1-1/R1p-1/R2p)
% I12 = b/mu*(1-1/R2p)*(1-1/R1p)*(beta1+beta2)/(beta1+beta2-mu-s1*s2*rho)-Ys(end,4) %OK

% E12 = [S12, Nequ/R2p-S12, Nequ/R1p - S12, S12 + Nequ*(1-1/R1p - 1/R2p)]
% Ys(end,:)

%Coinfection endemic state 
gamma1t = gamma1p-s1*s2*rho;
gamma2t = gamma2p-s1*s2*rho;
num     = gamma1t/R2p+gamma2t/R1p+mu+s1*s2*rho;
denum   = beta1+beta2-mu-s1*s2*rho;
S12     = Nequ*num/denum;
[S12, Nequ/R2p-S12, Nequ/R1p-S12, S12 + Nequ*(1-1/R1p - 1/R2p)]
Ys(end,:)

% %% Utility function %(case s1=0)
% clear all
% syms beta1 beta2 gamma1 gamma2 mu s1 s2 rho b s c
% %[beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');
% syms rho
% %mu=0;
% s1=1
% s2=1
% 
% R1p = beta1/(gamma1 + s1*rho + mu);
% R2p = beta2/(gamma2 + s2*rho + mu);
% gamma2p = gamma2 + s2*rho;
% gamma1p = gamma1 + s1*rho;
% 
% num = gamma2p/R1p + gamma1p/R2p + mu + (s1*s2*rho)*(1-1/R1p-1/R2p); %??
% denom = (beta1+beta2-mu-s1*s2*rho);
% s12 = num/denom;
% 
% U = rho*((1-s12)-c);
% dU = diff(U,rho)
% sol = solve(dU==0,rho) 
% 
% solve(simplify(dU)==0, 'maxdegree', 3)
% %cas general, premiere racine qui n'est pas celle qu'on cherche a priori
% % ((((2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)^3/(216*beta2^3*rho^9*s2^6) - (beta1^2*gamma1*gamma2 - beta1^3*beta2 - 2*beta1^2*beta2^2 - beta1*beta2^3 + beta2^2*gamma1*gamma2 - 2*beta1*beta2*mu^2 + 3*beta1*beta2^2*mu + 3*beta1^2*beta2*mu - beta1*gamma1*mu^2 + beta1^2*gamma1*mu - beta2*gamma2*mu^2 + beta2^2*gamma2*mu + 2*beta1^2*gamma1*rho*s2 + 2*beta2^2*gamma1*rho*s2 - 2*beta2*mu^2*rho*s2 + 2*beta2^2*mu*rho*s2 + 2*beta1*beta2*gamma1*gamma2 + beta1*beta2*gamma1*mu + beta1*beta2*gamma2*mu - beta1*gamma1*gamma2*mu - beta2*gamma1*gamma2*mu + 4*beta1*beta2*gamma1*rho*s2 + 2*beta1*beta2*mu*rho*s2 - 2*beta1*gamma1*mu*rho*s2 - 2*beta2*gamma1*mu*rho*s2)/(4*beta2*rho^3*s2^2) + ((2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)*(2*beta1*mu^2*rho - 2*beta1^2*mu*rho - 3*beta1^2*rho^2*s2 - 3*beta2^2*rho^2*s2 + 3*beta1^2*rho^2*s2^2 - 2*beta1^2*gamma2*rho - 2*beta2^2*gamma2*rho - 2*beta1*beta2*mu*rho + 2*beta1*gamma2*mu*rho + 2*beta2*gamma2*mu*rho - 6*beta1*beta2*rho^2*s2 - 4*beta1*beta2^2*rho*s2 - 4*beta1^2*beta2*rho*s2 + 2*beta1^2*gamma2*rho*s2 + 2*beta2^2*gamma1*rho*s2 + 3*beta1*mu*rho^2*s2 - 2*beta1*mu^2*rho*s2 + 2*beta1^2*mu*rho*s2 + 3*beta2*mu*rho^2*s2 - 2*beta2*mu^2*rho*s2 + 2*beta2^2*mu*rho*s2 + 3*beta1*beta2*rho^2*s2^2 + beta1*gamma1*rho^2*s2^2 + beta2*gamma1*rho^2*s2^2 - 3*beta1*mu*rho^2*s2^2 + beta2*mu*rho^2*s2^2 - 4*beta1*beta2*gamma2*rho + 2*beta1*beta2*gamma1*rho*s2 + 2*beta1*beta2*gamma2*rho*s2 + 8*beta1*beta2*mu*rho*s2 - 2*beta1*gamma2*mu*rho*s2 - 2*beta2*gamma1*mu*rho*s2))/(24*beta2^2*rho^6*s2^4))^2 - ((2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)^2/(36*beta2^2*rho^6*s2^4) + (2*beta1*mu^2*rho - 2*beta1^2*mu*rho - 3*beta1^2*rho^2*s2 - 3*beta2^2*rho^2*s2 + 3*beta1^2*rho^2*s2^2 - 2*beta1^2*gamma2*rho - 2*beta2^2*gamma2*rho - 2*beta1*beta2*mu*rho + 2*beta1*gamma2*mu*rho + 2*beta2*gamma2*mu*rho - 6*beta1*beta2*rho^2*s2 - 4*beta1*beta2^2*rho*s2 - 4*beta1^2*beta2*rho*s2 + 2*beta1^2*gamma2*rho*s2 + 2*beta2^2*gamma1*rho*s2 + 3*beta1*mu*rho^2*s2 - 2*beta1*mu^2*rho*s2 + 2*beta1^2*mu*rho*s2 + 3*beta2*mu*rho^2*s2 - 2*beta2*mu^2*rho*s2 + 2*beta2^2*mu*rho*s2 + 3*beta1*beta2*rho^2*s2^2 + beta1*gamma1*rho^2*s2^2 + beta2*gamma1*rho^2*s2^2 - 3*beta1*mu*rho^2*s2^2 + beta2*mu*rho^2*s2^2 - 4*beta1*beta2*gamma2*rho + 2*beta1*beta2*gamma1*rho*s2 + 2*beta1*beta2*gamma2*rho*s2 + 8*beta1*beta2*mu*rho*s2 - 2*beta1*gamma2*mu*rho*s2 - 2*beta2*gamma1*mu*rho*s2)/(6*beta2*rho^3*s2^2))^3)^(1/2) - (beta1^2*gamma1*gamma2 - beta1^3*beta2 - 2*beta1^2*beta2^2 - beta1*beta2^3 + beta2^2*gamma1*gamma2 - 2*beta1*beta2*mu^2 + 3*beta1*beta2^2*mu + 3*beta1^2*beta2*mu - beta1*gamma1*mu^2 + beta1^2*gamma1*mu - beta2*gamma2*mu^2 + beta2^2*gamma2*mu + 2*beta1^2*gamma1*rho*s2 + 2*beta2^2*gamma1*rho*s2 - 2*beta2*mu^2*rho*s2 + 2*beta2^2*mu*rho*s2 + 2*beta1*beta2*gamma1*gamma2 + beta1*beta2*gamma1*mu + beta1*beta2*gamma2*mu - beta1*gamma1*gamma2*mu - beta2*gamma1*gamma2*mu + 4*beta1*beta2*gamma1*rho*s2 + 2*beta1*beta2*mu*rho*s2 - 2*beta1*gamma1*mu*rho*s2 - 2*beta2*gamma1*mu*rho*s2)/(4*beta2*rho^3*s2^2) + (2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)^3/(216*beta2^3*rho^9*s2^6) + ((2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)*(2*beta1*mu^2*rho - 2*beta1^2*mu*rho - 3*beta1^2*rho^2*s2 - 3*beta2^2*rho^2*s2 + 3*beta1^2*rho^2*s2^2 - 2*beta1^2*gamma2*rho - 2*beta2^2*gamma2*rho - 2*beta1*beta2*mu*rho + 2*beta1*gamma2*mu*rho + 2*beta2*gamma2*mu*rho - 6*beta1*beta2*rho^2*s2 - 4*beta1*beta2^2*rho*s2 - 4*beta1^2*beta2*rho*s2 + 2*beta1^2*gamma2*rho*s2 + 2*beta2^2*gamma1*rho*s2 + 3*beta1*mu*rho^2*s2 - 2*beta1*mu^2*rho*s2 + 2*beta1^2*mu*rho*s2 + 3*beta2*mu*rho^2*s2 - 2*beta2*mu^2*rho*s2 + 2*beta2^2*mu*rho*s2 + 3*beta1*beta2*rho^2*s2^2 + beta1*gamma1*rho^2*s2^2 + beta2*gamma1*rho^2*s2^2 - 3*beta1*mu*rho^2*s2^2 + beta2*mu*rho^2*s2^2 - 4*beta1*beta2*gamma2*rho + 2*beta1*beta2*gamma1*rho*s2 + 2*beta1*beta2*gamma2*rho*s2 + 8*beta1*beta2*mu*rho*s2 - 2*beta1*gamma2*mu*rho*s2 - 2*beta2*gamma1*mu*rho*s2))/(24*beta2^2*rho^6*s2^4))^(1/3) + ((2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)^2/(36*beta2^2*rho^6*s2^4) + (2*beta1*mu^2*rho - 2*beta1^2*mu*rho - 3*beta1^2*rho^2*s2 - 3*beta2^2*rho^2*s2 + 3*beta1^2*rho^2*s2^2 - 2*beta1^2*gamma2*rho - 2*beta2^2*gamma2*rho - 2*beta1*beta2*mu*rho + 2*beta1*gamma2*mu*rho + 2*beta2*gamma2*mu*rho - 6*beta1*beta2*rho^2*s2 - 4*beta1*beta2^2*rho*s2 - 4*beta1^2*beta2*rho*s2 + 2*beta1^2*gamma2*rho*s2 + 2*beta2^2*gamma1*rho*s2 + 3*beta1*mu*rho^2*s2 - 2*beta1*mu^2*rho*s2 + 2*beta1^2*mu*rho*s2 + 3*beta2*mu*rho^2*s2 - 2*beta2*mu^2*rho*s2 + 2*beta2^2*mu*rho*s2 + 3*beta1*beta2*rho^2*s2^2 + beta1*gamma1*rho^2*s2^2 + beta2*gamma1*rho^2*s2^2 - 3*beta1*mu*rho^2*s2^2 + beta2*mu*rho^2*s2^2 - 4*beta1*beta2*gamma2*rho + 2*beta1*beta2*gamma1*rho*s2 + 2*beta1*beta2*gamma2*rho*s2 + 8*beta1*beta2*mu*rho*s2 - 2*beta1*gamma2*mu*rho*s2 - 2*beta2*gamma1*mu*rho*s2)/(6*beta2*rho^3*s2^2))/((((2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)^3/(216*beta2^3*rho^9*s2^6) - (beta1^2*gamma1*gamma2 - beta1^3*beta2 - 2*beta1^2*beta2^2 - beta1*beta2^3 + beta2^2*gamma1*gamma2 - 2*beta1*beta2*mu^2 + 3*beta1*beta2^2*mu + 3*beta1^2*beta2*mu - beta1*gamma1*mu^2 + beta1^2*gamma1*mu - beta2*gamma2*mu^2 + beta2^2*gamma2*mu + 2*beta1^2*gamma1*rho*s2 + 2*beta2^2*gamma1*rho*s2 - 2*beta2*mu^2*rho*s2 + 2*beta2^2*mu*rho*s2 + 2*beta1*beta2*gamma1*gamma2 + beta1*beta2*gamma1*mu + beta1*beta2*gamma2*mu - beta1*gamma1*gamma2*mu - beta2*gamma1*gamma2*mu + 4*beta1*beta2*gamma1*rho*s2 + 2*beta1*beta2*mu*rho*s2 - 2*beta1*gamma1*mu*rho*s2 - 2*beta2*gamma1*mu*rho*s2)/(4*beta2*rho^3*s2^2) + ((2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)*(2*beta1*mu^2*rho - 2*beta1^2*mu*rho - 3*beta1^2*rho^2*s2 - 3*beta2^2*rho^2*s2 + 3*beta1^2*rho^2*s2^2 - 2*beta1^2*gamma2*rho - 2*beta2^2*gamma2*rho - 2*beta1*beta2*mu*rho + 2*beta1*gamma2*mu*rho + 2*beta2*gamma2*mu*rho - 6*beta1*beta2*rho^2*s2 - 4*beta1*beta2^2*rho*s2 - 4*beta1^2*beta2*rho*s2 + 2*beta1^2*gamma2*rho*s2 + 2*beta2^2*gamma1*rho*s2 + 3*beta1*mu*rho^2*s2 - 2*beta1*mu^2*rho*s2 + 2*beta1^2*mu*rho*s2 + 3*beta2*mu*rho^2*s2 - 2*beta2*mu^2*rho*s2 + 2*beta2^2*mu*rho*s2 + 3*beta1*beta2*rho^2*s2^2 + beta1*gamma1*rho^2*s2^2 + beta2*gamma1*rho^2*s2^2 - 3*beta1*mu*rho^2*s2^2 + beta2*mu*rho^2*s2^2 - 4*beta1*beta2*gamma2*rho + 2*beta1*beta2*gamma1*rho*s2 + 2*beta1*beta2*gamma2*rho*s2 + 8*beta1*beta2*mu*rho*s2 - 2*beta1*gamma2*mu*rho*s2 - 2*beta2*gamma1*mu*rho*s2))/(24*beta2^2*rho^6*s2^4))^2 - ((2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)^2/(36*beta2^2*rho^6*s2^4) + (2*beta1*mu^2*rho - 2*beta1^2*mu*rho - 3*beta1^2*rho^2*s2 - 3*beta2^2*rho^2*s2 + 3*beta1^2*rho^2*s2^2 - 2*beta1^2*gamma2*rho - 2*beta2^2*gamma2*rho - 2*beta1*beta2*mu*rho + 2*beta1*gamma2*mu*rho + 2*beta2*gamma2*mu*rho - 6*beta1*beta2*rho^2*s2 - 4*beta1*beta2^2*rho*s2 - 4*beta1^2*beta2*rho*s2 + 2*beta1^2*gamma2*rho*s2 + 2*beta2^2*gamma1*rho*s2 + 3*beta1*mu*rho^2*s2 - 2*beta1*mu^2*rho*s2 + 2*beta1^2*mu*rho*s2 + 3*beta2*mu*rho^2*s2 - 2*beta2*mu^2*rho*s2 + 2*beta2^2*mu*rho*s2 + 3*beta1*beta2*rho^2*s2^2 + beta1*gamma1*rho^2*s2^2 + beta2*gamma1*rho^2*s2^2 - 3*beta1*mu*rho^2*s2^2 + beta2*mu*rho^2*s2^2 - 4*beta1*beta2*gamma2*rho + 2*beta1*beta2*gamma1*rho*s2 + 2*beta1*beta2*gamma2*rho*s2 + 8*beta1*beta2*mu*rho*s2 - 2*beta1*gamma2*mu*rho*s2 - 2*beta2*gamma1*mu*rho*s2)/(6*beta2*rho^3*s2^2))^3)^(1/2) - (beta1^2*gamma1*gamma2 - beta1^3*beta2 - 2*beta1^2*beta2^2 - beta1*beta2^3 + beta2^2*gamma1*gamma2 - 2*beta1*beta2*mu^2 + 3*beta1*beta2^2*mu + 3*beta1^2*beta2*mu - beta1*gamma1*mu^2 + beta1^2*gamma1*mu - beta2*gamma2*mu^2 + beta2^2*gamma2*mu + 2*beta1^2*gamma1*rho*s2 + 2*beta2^2*gamma1*rho*s2 - 2*beta2*mu^2*rho*s2 + 2*beta2^2*mu*rho*s2 + 2*beta1*beta2*gamma1*gamma2 + beta1*beta2*gamma1*mu + beta1*beta2*gamma2*mu - beta1*gamma1*gamma2*mu - beta2*gamma1*gamma2*mu + 4*beta1*beta2*gamma1*rho*s2 + 2*beta1*beta2*mu*rho*s2 - 2*beta1*gamma1*mu*rho*s2 - 2*beta2*gamma1*mu*rho*s2)/(4*beta2*rho^3*s2^2) + (2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)^3/(216*beta2^3*rho^9*s2^6) + ((2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)*(2*beta1*mu^2*rho - 2*beta1^2*mu*rho - 3*beta1^2*rho^2*s2 - 3*beta2^2*rho^2*s2 + 3*beta1^2*rho^2*s2^2 - 2*beta1^2*gamma2*rho - 2*beta2^2*gamma2*rho - 2*beta1*beta2*mu*rho + 2*beta1*gamma2*mu*rho + 2*beta2*gamma2*mu*rho - 6*beta1*beta2*rho^2*s2 - 4*beta1*beta2^2*rho*s2 - 4*beta1^2*beta2*rho*s2 + 2*beta1^2*gamma2*rho*s2 + 2*beta2^2*gamma1*rho*s2 + 3*beta1*mu*rho^2*s2 - 2*beta1*mu^2*rho*s2 + 2*beta1^2*mu*rho*s2 + 3*beta2*mu*rho^2*s2 - 2*beta2*mu^2*rho*s2 + 2*beta2^2*mu*rho*s2 + 3*beta1*beta2*rho^2*s2^2 + beta1*gamma1*rho^2*s2^2 + beta2*gamma1*rho^2*s2^2 - 3*beta1*mu*rho^2*s2^2 + beta2*mu*rho^2*s2^2 - 4*beta1*beta2*gamma2*rho + 2*beta1*beta2*gamma1*rho*s2 + 2*beta1*beta2*gamma2*rho*s2 + 8*beta1*beta2*mu*rho*s2 - 2*beta1*gamma2*mu*rho*s2 - 2*beta2*gamma1*mu*rho*s2))/(24*beta2^2*rho^6*s2^4))^(1/3) + (2*beta1*rho^3*s2^2 + 3*beta2^2*rho^2*s2 - 2*beta1*rho^3*s2^3 + 2*beta2*rho^3*s2^2 + 3*beta1*beta2*rho^2*s2 + beta1*gamma2*rho^2*s2 + beta2*gamma2*rho^2*s2 + beta1*mu*rho^2*s2 - 3*beta2*mu*rho^2*s2 + 2*beta1*beta2*rho^2*s2^2 - beta1*gamma2*rho^2*s2^2 - beta2*gamma1*rho^2*s2^2 - beta1*mu*rho^2*s2^2 - beta2*mu*rho^2*s2^2)/(6*beta2*rho^3*s2^2)
% 
% 
% %s1= 0
% A = (mu + (gamma2*(gamma1 + mu))/beta1 + (gamma1*(gamma2 + mu))/beta2 - (beta1+beta2-mu));
% sol = -A/(2*(s2*(gamma1 + mu)/beta1 + (gamma1*s2)/beta2)); %OK
% 
% %sol pour s1=1=s2 mu=0
% R1 = beta1/(gamma1+mu);
% R2 = beta2/(gamma2+mu);
% r = R1*R2; % a voir...
% 
% sol = (beta1+beta2)*(1-sqrt(1-(r-1)/(2*r+beta1/gamma2+beta2/gamma1)))
% 
% %sol pour mu diff de 0, s1=s2=1
% (2*beta1*beta2^2 + 2*beta1^2*beta2 + beta1^2*gamma1 + beta2^2*gamma2 - ((beta1 + beta2)*(beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)*(beta1 + beta2 - mu))^(1/2) + beta1*beta2*gamma1 + beta1*beta2*gamma2 - 2*beta1*beta2*mu - beta1*gamma1*mu - beta2*gamma2*mu)/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)
% (2*beta1*beta2^2 + 2*beta1^2*beta2 + beta1^2*gamma1 + beta2^2*gamma2 + ((beta1 + beta2)*(beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)*(beta1 + beta2 - mu))^(1/2) + beta1*beta2*gamma1 + beta1*beta2*gamma2 - 2*beta1*beta2*mu - beta1*gamma1*mu - beta2*gamma2*mu)/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)
% 
% r=beta1/gamma1*beta2/gamma2;
% (beta1+beta2-mu)-sqrt((beta1+beta2)*(beta1+beta2-mu)*(1-(-1+r)/(2*r+beta1/gamma2+beta2/gamma1)))
% 
% %sol pour s1=s2=1 cdiff de 0
% % (2*beta1*beta2^2 + 2*beta1^2*beta2 -...
% %     ((beta1 + beta2)*(beta1 + gamma2)*(beta2 + gamma1)*(beta1 + beta2 - mu)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2 - beta1*beta2*c))^(1/2) +...
% %     beta1^2*gamma1 + beta2^2*gamma2 + beta1*beta2*gamma1 + beta1*beta2*gamma2 - 2*beta1*beta2*mu - beta1*gamma1*mu - beta2*gamma2*mu - beta1*beta2^2*c - beta1^2*beta2*c + beta1*beta2*c*mu)/...
% %     (beta1*gamma1 + beta2*gamma2 - beta1*beta2*(c - 2))
% % (2*beta1*beta2^2 + 2*beta1^2*beta2 +...
% %     ((beta1 + beta2)*(beta1 + gamma2)*(beta2 + gamma1)*(beta1 + beta2 - mu)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2 - beta1*beta2*c))^(1/2) +...
% %     beta1^2*gamma1 + beta2^2*gamma2 + beta1*beta2*gamma1 + beta1*beta2*gamma2 - 2*beta1*beta2*mu - beta1*gamma1*mu - beta2*gamma2*mu - beta1*beta2^2*c - beta1^2*beta2*c + beta1*beta2*c*mu)/...
% %     (beta1*gamma1 + beta2*gamma2 - beta1*beta2*(c - 2))
% 
% %verif
% r=beta1*beta2/(gamma1*gamma2);
% rhohat = (beta1+beta2-mu) - sqrt((beta1+beta2)*(beta1+beta2-mu)*((beta1/gamma2+beta2/gamma1+r+1)/(beta1/gamma2+beta2/gamma1+(2-c)*r)))
% 
% 
% %% Utility function (case s1=s2=1,mu=0)
% clear all
% syms beta1 beta2 gamma1 gamma2 mu s1 s2 rho b
% s1=1; s2=1; mu=0;
% 
% R1 = beta1./(gamma1 + mu);
% R2 = beta2./(gamma2 + mu);
% 
% gamma2p = gamma2 + s2.*rho;
% gamma1p = gamma1 + s1.*rho;
% R1p = beta1./(gamma1p + mu);
% R2p = beta2./(gamma2p + mu);
% gamma1t = gamma1p - s1*s2*rho;
% gamma2t = gamma2p - s1*s2*rho;
% 
% num = gamma2t./R1p + gamma1t./R2p + mu + s1*s2.*rho;
% denom = (beta1+beta2-mu-s1*s2*rho);
% s12 = num./denom;
% 
% U = rho.*(1-s12);
% 
% dU = diff(U,rho);
% sol = solve(dU==0,rho)
% 
% %%
% [beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');
% s1=1;s2=1;mu=0;
% R1 = beta1/(gamma1 + mu);
% R2 = beta2/(gamma2 + mu);
% sol1 = (2*beta1*beta2^2 + 2*beta1^2*beta2 + beta1^2*gamma1 + beta2^2*gamma2 + beta1*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta2*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta1*beta2*gamma1 + beta1*beta2*gamma2)/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)
% sol2 = (2*beta1*beta2^2 + 2*beta1^2*beta2 + beta1^2*gamma1 + beta2^2*gamma2 - beta1*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) - beta2*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta1*beta2*gamma1 + beta1*beta2*gamma2)/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)
% 
% r=R1*R2;
% sol3 = (beta1+beta2)*(1-sqrt(1- (r-1)/(2*r+beta1/gamma2+beta2/gamma1)))%idem sol2
% sol4 = (beta1+beta2)*(1+sqrt(1- (r-1)/(2*r+beta1/gamma2+beta2/gamma1)))%idem sol1
% 
% %U12(sol3)
% %((1 - ((beta1*beta2)/(gamma1*gamma2) - 1)/(beta1/gamma2 + beta2/gamma1 + (2*beta1*beta2)/(gamma1*gamma2)))^(1/2) - 1)*(beta1 + beta2)*(((gamma2*(gamma1 - ((1 - ((beta1*beta2)/(gamma1*gamma2) - 1)/(beta1/gamma2 + beta2/gamma1 + (2*beta1*beta2)/(gamma1*gamma2)))^(1/2) - 1)*(beta1 + beta2)))/beta1 - ((1 - ((beta1*beta2)/(gamma1*gamma2) - 1)/(beta1/gamma2 + beta2/gamma1 + (2*beta1*beta2)/(gamma1*gamma2)))^(1/2) - 1)*(beta1 + beta2) + (gamma1*(gamma2 - ((1 - ((beta1*beta2)/(gamma1*gamma2) - 1)/(beta1/gamma2 + beta2/gamma1 + (2*beta1*beta2)/(gamma1*gamma2)))^(1/2) - 1)*(beta1 + beta2)))/beta2)/(beta1 + beta2 + ((1 - ((beta1*beta2)/(gamma1*gamma2) - 1)/(beta1/gamma2 + beta2/gamma1 + (2*beta1*beta2)/(gamma1*gamma2)))^(1/2) - 1)*(beta1 + beta2)) - 1)
% %((beta1 + beta2)*((((beta1 + gamma2)*(beta2 + gamma1))/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) - 1)*(beta1*beta2 + beta1*gamma1 + beta2*gamma2 + gamma1*gamma2 - 2*beta1*beta2*(((beta1 + gamma2)*(beta2 + gamma1))/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) - beta1*gamma1*(((beta1 + gamma2)*(beta2 + gamma1))/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) - beta2*gamma2*(((beta1 + gamma2)*(beta2 + gamma1))/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2)))/(beta1*beta2*(((beta1 + gamma2)*(beta2 + gamma1))/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2))
% 
% 
% %%
% 
% % Verification manuelle
% % dUbis = 1 - (gamma1*gamma2*(1/beta1+1/beta2)+rho*(gamma1/beta2+gamma2/beta1+1))/(beta1+beta2-rho) - rho*((gamma1/beta2+gamma2/beta1+1)*(beta1+beta2-rho)+gamma1*gamma2*(1/beta1+1/beta2)+rho*(gamma1/beta2+gamma2/beta1+1))/(beta1+beta2-rho)^2
% % solve(dUbis==0,rho)
% % (2*beta1*beta2^2 + 2*beta1^2*beta2 + beta1^2*gamma1 + beta2^2*gamma2 + beta1*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta2*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta1*beta2*gamma1 + beta1*beta2*gamma2)/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)
% % (2*beta1*beta2^2 + 2*beta1^2*beta2 + beta1^2*gamma1 + beta2^2*gamma2 - beta1*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) - beta2*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta1*beta2*gamma1 + beta1*beta2*gamma2)/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)
%  
% %Verif avec le polynome
% % pol = (gamma1/beta2+gamma2/beta1+2)*rho^2 +...
% %     rho*(-2*(beta1+beta2) - 2*(beta1+beta2)*(gamma1/beta2+gamma2/beta1+1)) +...
% %     (beta1+beta2)^2  - (beta1+beta2)*gamma1*gamma2*(1/beta1+1/beta2)
% % A = (gamma1/beta2+gamma2/beta1+2);
% % B = (-2*(beta1+beta2) - 2*(beta1+beta2)*(gamma1/beta2+gamma2/beta1+1));
% % C = (beta1+beta2)^2  - (beta1+beta2)*gamma1*gamma2*(1/beta1+1/beta2);
% % dP = solve(pol==0,rho)
% %  (2*beta1*beta2^2 + 2*beta1^2*beta2 + beta1^2*gamma1 + beta2^2*gamma2 + beta1*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta2*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta1*beta2*gamma1 + beta1*beta2*gamma2)/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)
% %  (2*beta1*beta2^2 + 2*beta1^2*beta2 + beta1^2*gamma1 + beta2^2*gamma2 - beta1*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) - beta2*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta1*beta2*gamma1 + beta1*beta2*gamma2)/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)
% % rho = B^2 -4*A*C
% % simplify(rho)
% % rho = (4*(beta1 + beta2)^2*(beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))/(beta1^2*beta2^2)
% % (2*beta1 + 2*beta2 - ((2*beta1 + 2*beta2 + (2*beta1 + 2*beta2)*(gamma2/beta1 + gamma1/beta2 + 1))^2 - ((beta1 + beta2)^2 - gamma1*gamma2*(beta1 + beta2)*(1/beta1 + 1/beta2))*((4*gamma2)/beta1 + (4*gamma1)/beta2 + 8))^(1/2) + (2*beta1 + 2*beta2)*(gamma2/beta1 + gamma1/beta2 + 1))/((2*gamma2)/beta1 + (2*gamma1)/beta2 + 4)
% % (2*beta1 + 2*beta2 + ((2*beta1 + 2*beta2 + (2*beta1 + 2*beta2)*(gamma2/beta1 + gamma1/beta2 + 1))^2 - ((beta1 + beta2)^2 - gamma1*gamma2*(beta1 + beta2)*(1/beta1 + 1/beta2))*((4*gamma2)/beta1 + (4*gamma1)/beta2 + 8))^(1/2) + (2*beta1 + 2*beta2)*(gamma2/beta1 + gamma1/beta2 + 1))/((2*gamma2)/beta1 + (2*gamma1)/beta2 + 4)
% 
% 
% %solution s1=s2=1, mu=0
% %ecriture 1
% (beta1+beta2)*(1-sqrt(1+(gamma1*gamma2/(beta1*beta2)-1)/(2+gamma1/beta2+gamma2/beta1)))
% 
% %ecriture 2
% a = gamma1/beta2 + gamma2/beta1 +1;
% c = beta1+beta2;
% b = gamma1*gamma2*(1/beta1+1/beta2);
% c*(1-sqrt((1+a)*(a+gamma1*gamma2/(beta1*beta2)))/(1+a))
% 
% 
% 
% 
% [beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');
% 
% %% case beta1=beta2 gamma1=gamma2
% beta2 = beta1;
% gamma2 = gamma1;
% 
% 
% (2*beta1*beta2^2 + 2*beta1^2*beta2 + beta1^2*gamma1 + beta2^2*gamma2 + beta1*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta2*((beta1 + gamma2)*(beta2 + gamma1)*(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta1*beta2*gamma1 + beta1*beta2*gamma2)/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)
% 
% 2*beta1 + sqrt(2*(beta1^2+gamma1*beta1))
% 2*beta1 - sqrt(2*(beta1^2+gamma1*beta1)) 
% 
% %% PLOT of U
% clear all
% close all
% arret =false;
% 
% while(~arret)
% [beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');
% % beta2 = beta1;
% % gamma2 = gamma1;
% mu=0;
% s2=1;
% s1=s2;
% R1 = beta1/(gamma1+mu);
% R2 = beta2/(gamma2+mu);
% 
% alpha1 = beta1/s1*(1-1/R1);
% alpha2 = beta2/s2*(1-1/R2);
% 
% if(alpha1<alpha2 & beta2<beta1 & gamma2>2*gamma1)
%     arret = true;
%     maxalpha = alpha2;
% 
%     vecRho = 0:maxalpha/5000:(maxalpha*1.2);
%     [U,P,vecR10,vecR20,P1,P2,P12] = utilityFunctionSIS2_V22(beta1,beta2,gamma1,gamma2,s1,s2,vecRho,b,mu)
% 
%     %maxU
%     [maxU,imax] = max(U);
%     rhomax = vecRho(imax);
% 
%     U12 = vecRho.*P12;
% 
%     figure(4)
%     plot(vecRho,max(U12,0), '--')
%     hold on;
%     plot(vecRho, max(vecRho.*P1,0),'--')
%     plot(vecRho, max(vecRho.*P2,0),'--')
%     plot(vecRho,U)
%     ylim([0 maxU])
% 
%     title([{'Utility U of the SISxSIS model function of \rho'},...
%             {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
%             num2str(gamma2), ' s_1=', num2str(s1),' s_2=', num2str(s2), '\mu=', num2str(mu)]}, ...
%             {[' R_1(0)=' num2str(round(R1,2)), ' R_2(0)=' num2str(round(R2,2)),...
%             ' \rho_{max}=', num2str(round(rhomax,2)),...
%             ' \rho_1''=',num2str(round(alpha1,2)),' \rho_2''=',num2str(round(alpha2,2))]}])
%     xlabel("\rho","fontweight","bold")
%     ylabel("U(\rho)","fontweight","bold")
%     zlabel('U_{12}=(\rho)(1-S/N)')
% 
%     legend('U_{12}=\rho \Pi_{12}(\rho)','U_1 = \rho \Pi_1(\rho)','U_2 = \rho \Pi_2(\rho)','U')
% end
% end
% 
% %% 
%  %U=0 : cas tous parametres egaux
% (beta1-gamma1)
% 
% %U=0 : cas s1=1=s
% s=s1;
% (beta1+beta2-gamma1/R2-gamma2/R1-2*mu)/(gamma1/beta2+gamma2/beta1+2)/s
% 
% % 04/09 Trouver argmax
% %Cas s1=s2=1, mu=0;
% r=R1*R2;
% rho12_th = (beta1+beta2)*(1- sqrt(1- (r-1)/(2*r + beta1/gamma2+beta2/gamma1)))
% [maxU,imax] = max(U12);
% vecRho(imax)
% 
% %% 07/06 - U12(rho12hat)
% syms beta1 beta2 gamma1 gamma2 mu s1 s2 rho 
% mu=0;s1=1;s2=1;
% gamma1p = gamma1 + s1*rho;      gamma2p = gamma2 + s2*rho;
% gamma1t = gamma1p - s1*s2*rho;  gamma2t = gamma2p - s1*s2*rho;
% R1p = beta1/(gamma1p+mu);       R2p = beta2/(gamma2p+mu); 
% R1 = beta1/(gamma1+mu);         R2 = beta2/(gamma2+mu); 
% 
% 
% %rho that maximizes U12
% r = R1*R2;
% rho12 = (beta1+beta2)*(1 - sqrt(1 - (r-1)/(2*r+beta1/gamma2+beta2/gamma1)));
%     
% rho = rho12;
% 
% %cas general
% %U12 = rho*(1- (gamma1t/R2p + gamma2t/R2p + mu + rho*s1*s2*(1-1/R1p-1/R2))/(beta1+beta2 - mu - s1*s2*rho));
% %cas mu=0;s1=1;s2=1;
% U12 = rho*(1 - (gamma1/R2p + gamma2/R2p + rho)/(beta1+beta2 - rho));
% 
% simplify(U12)
% %(((((beta1 + gamma2)*(beta2 + gamma1))/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) - 1)*(beta1*beta2 + gamma1*gamma2 + gamma1*rho + gamma2*rho - 2*beta2^2*(((beta1 + gamma2)*(beta2 + gamma1))/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2) + beta2^2 + gamma2^2 - 2*beta1*beta2*(((beta1 + gamma2)*(beta2 + gamma1))/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2)))/(beta2*(((beta1 + gamma2)*(beta2 + gamma1))/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2))^(1/2))
% 
% 
% %% 07/06 - Ui(rhoi)
% 
% %rho that maximizes rhoi
% rhoi = beta1/(2*s1)*(1-1/R1);
% Ui = beta1/(4*s1)*(1-1/R1)^2;
% 
% solve(U12-Ui==0, gamma2)
% 
% sol = solve(rho12 - rhoi==0,gamma2)
% %under some conditions
% a1 = (beta1*(- 2*beta1^2*beta2 - 3*beta1^2*gamma1 + 2*beta1*gamma1^2 + 4*beta2^3 + 8*beta2^2*gamma1 + 6*beta2*gamma1^2 + gamma1^3))/(3*beta1^2*beta2 + 4*beta1^2*gamma1 + 4*beta1*beta2^2 + 6*beta1*beta2*gamma1 - beta2*gamma1^2)
% a2 = -(beta1*(beta2/gamma1 + ((rhoi/(beta1 + beta2) - 1)^2 - 1)*((2*beta2)/gamma1 + 1)))/((beta2*((rhoi/(beta1 + beta2) - 1)^2 - 1))/gamma1 - 1)
% 
% 
% D = rhoi^2/(beta1+beta2)^2 - rhoi/(beta1+beta2);
% a3 = beta1*(D+beta2/gamma1+2)/(-beta2*D/gamma1+1)
% 
% %% 11/09 max de U quand beta1=beta2, gamma1=gamma2
% clear all; close all;
% [beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true')
% 
% %beta2=beta1; gamma2=gamma1;
% s1=1;s2=1;mu=0
% 
% alpha1 = beta1/s1*(1-(gamma1+mu)/beta1);
% vecRho = 0:alpha1/1000:alpha1*1.2;
% [U,P,R1p,R2p,P1,P2,P12] = utilityFunctionSIS2_V22(beta1,beta2,gamma1,gamma2,s1,s2,vecRho,b,mu);
% 
% plot(vecRho, U)
% 
% [Umax,imax] = max(U);
% rhomax = vecRho(imax)
% 
% rhoth = 2*beta1*(1-sqrt((1+gamma1/beta1)/2))
% 
% rhoth*(1-P12)
% %simplify(rhoth*(1-P12))
% %-(2*beta1*((2*(gamma1 - 2*beta1*((gamma1/(2*beta1) + 1/2)^(1/2) - 1))^2)/beta1 + (2*(gamma1 - 2*beta1*((gamma1/(2*beta1) + 1/2)^(1/2) - 1))^2*((gamma1/(2*beta1) + 1/2)^(1/2) - 1)*((2*(gamma1 - 2*beta1*((gamma1/(2*beta1) + 1/2)^(1/2) - 1)))/beta1 - 1))/beta1)*((gamma1/(2*beta1) + 1/2)^(1/2) - 1))/(2*beta1 + (2*(gamma1 - 2*beta1*((gamma1/(2*beta1) + 1/2)^(1/2) - 1))^2*((gamma1/(2*beta1) + 1/2)^(1/2) - 1))/beta1)
%   %kjjhg
%   
%   
% %% 04/02 comparison between rho12ES and rho2/2 (s1=s2=1)
% syms beta1 beta2 gamma1 gamma2 mu
% mu=0;
% %[beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true')
% r = beta1*beta2/(gamma1*gamma2);
% rho12 = (beta1+beta2-mu)-sqrt((beta1+beta2-mu)*(beta1+beta2)*(1-(r-1)/(2*r+beta1/gamma2+beta2/gamma1)))
% rho2 = beta2/2*(1-(gamma2+mu)/beta2)
% 
% sol = solve(rho12-rho2==0,gamma1)
% 
% %gamma1
% %-(beta2*(beta1/gamma2 + ((beta1 + beta2 - mu + (beta2*((gamma2 + mu)/beta2 - 1))/2)^2/((beta1 + beta2)*(beta1 + beta2 - mu)) - 1)*((2*beta1)/gamma2 + 1)))/((beta1*((beta1 + beta2 - mu + (beta2*((gamma2 + mu)/beta2 - 1))/2)^2/((beta1 + beta2)*(beta1 + beta2 - mu)) - 1))/gamma2 - 1)
% (beta2*(4*beta1^3 + 8*beta1^2*gamma2 - 2*beta1*beta2^2 + 6*beta1*gamma2^2 - 3*beta2^2*gamma2 + 2*beta2*gamma2^2 + gamma2^3))/(4*beta1^2*beta2 + 3*beta1*beta2^2 + 6*beta1*beta2*gamma2 - beta1*gamma2^2 + 4*beta2^2*gamma2)
% 
% %mu
% % (2*beta1^2*beta2 + beta2*gamma2^2 - beta2^2*gamma2 + 2*beta1*(-beta2*(beta1 + beta2)*(beta2 + gamma1)*(beta1 - beta2 - gamma1 + gamma2))^(1/2) + 2*gamma2*(-beta2*(beta1 + beta2)*(beta2 + gamma1)*(beta1 - beta2 - gamma1 + gamma2))^(1/2) - beta1*beta2*gamma1 + 2*beta1*beta2*gamma2 - beta1*gamma1*gamma2 - 2*beta2*gamma1*gamma2)/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)
% % -(beta2^2*gamma2 - beta2*gamma2^2 - 2*beta1^2*beta2 + 2*beta1*(-beta2*(beta1 + beta2)*(beta2 + gamma1)*(beta1 - beta2 - gamma1 + gamma2))^(1/2) + 2*gamma2*(-beta2*(beta1 + beta2)*(beta2 + gamma1)*(beta1 - beta2 - gamma1 + gamma2))^(1/2) + beta1*beta2*gamma1 - 2*beta1*beta2*gamma2 + beta1*gamma1*gamma2 + 2*beta2*gamma1*gamma2)/(2*beta1*beta2 + beta1*gamma1 + beta2*gamma2)
% 
% %% numeriquement (seulement si s1=s2=1 ??)
% clear all
% k=0;cond=1;
% while(cond & k<= 10000)
%     [beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');
%     s1=1;s2=1;
%     alpha1 = beta1/s1*(1-(gamma1+mu)/beta1);
%     alpha2 = beta2/s2*(1-(gamma2+mu)/beta2);
% 
%     if(alpha1<alpha2  && alpha1>0)
%         k=k+1;
% 
%         r = beta1*beta2/(gamma1*gamma2);
%         rho12 = (beta1+beta2-mu)-sqrt((beta1+beta2-mu)*(beta1+beta2)*(1-(r-1)/(2*r+beta1/gamma2+beta2/gamma1)));
%         rho2 = beta2/2*(1-(gamma2+mu)/beta2);
%         rho1hat = alpha1/2;
%         cond = rho12>rho1hat;
%          
%         %
%         rho = 0:alpha1/100:alpha1;
%         R1p = beta1./(gamma1 + s1*rho + mu);
%         R2p = beta2./(gamma2 + s2*rho + mu);
%         gamma2p = gamma2 + s2*rho;
%         gamma1p = gamma1 + s1*rho;
% 
%         num = gamma2p./R1p + gamma1p./R2p + mu + (s1*s2*rho).*(1-1./R1p-1./R2p);
%         denom = (beta1+beta2-mu-s1.*s2.*rho);
%         S12 = b./mu.*num./denom;
%         U = rho.*(1-S12*mu/b);
%         
%         [maxU,imax] = max(U);
%         rho12 = rho(imax);
%         cond = rho12>=rho1hat;     
%     end
% end
% 
% 
% %% zones combined testing
% clear all;close all;
% len=200;
% gamma2=356/21; gamma1=12/1; 
% s1=1;s2=1;
% mu=1/35;
% 
% R10 = linspace(1.1,5,len);
% R20 = linspace(1.1,5,len);
% %R10 = 3.86;R20=4.49
% 
% [R1,R2] = meshgrid(R10,R20);
% 
% BETA1 = R1.*(gamma1+mu);
% BETA2 = R2.*(gamma2+mu);
% ALPHA1 = BETA1./s1.*(1-1./R1);
% ALPHA2 = BETA2./s2.*(1-1./R2);
% c11   = 1./R1-1;
% c22   = 1./R2-1;
% 
% N = 2;
% CEL1=[];
%         
% %calculer c1_12 et c2_12;
% for i=1:(len)
%     i
%     for j=1:(len)
%         c11 = 1/R1(i,j)-1; c01 = -c11;
%         c22 = 1/R2(i,j)-1; c02 = -c22;
%         mincnn = min(c11,c22);maxc0n = max(c01,c02);
%         maxcnn = max(c11,c22);
%         %interval of c
%         vecC = linspace(1.5*mincnn,maxc0n,50);
%         diff = 0;step=1;
%         while (diff<step)
%             
%             step = vecC(2) - vecC(1);
% 
%             [tab,tabco,tabcn] = findRhohat(N,0,0,[BETA1(i,j),gamma1,1;BETA2(i,j),gamma2,1],[],[],mu,5,vecC);        
%             vecAlpha = [ALPHA1(i,j),ALPHA2(i,j)]; 
%             pres = 0.1*abs(vecC(1)-vecC(end));
%             cs1dis = findThresholds(2,0,0, tab, vecAlpha, [BETA1(i,j),BETA2(i,j)], [gamma1,gamma2], [s1,s2], 5, mu,vecC);
% 
%             if(size(cs1dis,2)==0)%did not find cs1dis in the interval vecC
%                 %then cs1dis<min(cnn);
%                 cs1dis=vecC(1);
%             end
%             
%             if(ALPHA1(i,j)<ALPHA2(i,j))
%                 %alpha1=alphaj %j : first disease eliminated
%                 cjj = c11;
%             else
%                 cjj = c22;
%             end
%             diff = abs(cs1dis-cjj); %affiner vecC
%             vecC = linspace(1.1*min(cs1dis,cjj),1.1*max(cs1dis,cjj),10);
%             if(diff<step)
%                 break;
%             end
%        end
%        CEL1(i,j) = cs1dis; %store costs for which 1 disease has been elminated
% 
%     end
% end
% 
% yR2 = (gamma1+mu)/(gamma2+mu)*(R10-1)+1;
% %%
% %which disease elminates the other ?
% % e.g. if c11<c112 and alpha1<alpha2 then combined testing eliminates disease 1
% % more easily thanks to disease 2
% close all;
% n=15;%nb de subdvisions colormap
% titre = ['IST1 : ', num2str(round(1/gamma1,2)),' years and IST 2 : ', num2str(round(1/gamma2,2)),' years'];
% 
% %%%%%% !!!!!!!!!!!!!! c11 et c22 ?? !!!!!!!!!!!
% %C11 = 1./R1-1; C22 = 1./R2-1;
% dis2drives1 = (ALPHA1<=ALPHA2).*(c11-CEL1)./c11;
% dis1drives2 = (ALPHA1>=ALPHA2).*(c22-CEL1)./c22;
% tot = dis2drives1+dis1drives2; minv = min(min(tot));maxv=max(max(tot));
% dis2drives1((ALPHA1>ALPHA2) & dis2drives1==0)=NaN;
% dis1drives2((ALPHA1<ALPHA2) & dis1drives2==0)=NaN;
% 
% %coordinates to plot text
% yR10end = (gamma1+mu)/(gamma2+mu)*(R10(end)-1)+1;
% xR20end = (gamma2+mu)/(gamma1+mu)*(R20(end)-1)+1;
% if(R10(end)>xR20end) %voir les notes du 14/04 (triangle superieur)
%    xG1 = (min(R10)+min(R10)+xR20end)/3;%pas tout a fait ça avec R10(end)
%    yG1 = (min(R20)+yR10end+max(R20))/3;
%    xG2 = (min(R10)+max(R10)+max(R10)+xR20end)/4;
%    yG2 = (min(R20)+min(R20)+max(R20)+max(R20))/4;
% else
%    xG1 = (min(R10)+min(R10)+R10(end)+R10(end))/4;%quadrilatere superieur
%    yG1 = (min(R20)+max(R20)+max(R20)+yR10end)/4;
%    xG2 = (min(R10)+max(R10)+max(R10))/3;
%    yG2 = (min(R20)+min(R20)+yR10end)/3;
% end
% 
% %plot : diff between cjj and cjij
% figure(1) 
% surf1=surf(R1,R2,tot);
% surf1.EdgeColor = 'none';
% text(xG1,yG1,maxv,'\rho\prime_1<\rho\prime_2','HorizontalAlignment','center','FontSize',20)
% text(xG2,yG2,maxv,'\rho\prime_1>\rho\prime_2','HorizontalAlignment','center','FontSize',20)
% new_cb = subdivisedColormap([[0.2,0.2,0.2];[1,1,1];[0.466,0.674,0.1880]],n, 'quad'); %2^n+1
% new_cb2 = new_cb(ceil((2^n+1)*0.5*(1-abs(minv)/maxv)):end,:);
% colormap(new_cb2)
% cb = colorbar;
% cb.YTick = [-0.5 0 0.5 1 1.5];
% cb.YTickLabel = {'-50%', 'No Change', '+50%','+100%', '+150%'};
% hold on;
% plot3(R10,yR2,maxv*ones(length(R10),1),'-k');
% ylim([1,max(R20)]);
% title([{'Is combined testing better than specific testing ?'},...
%     {'$(c_{j2} - c_{j1\times2})/c_{j2}$ such that $\rho_j\prime<\rho_i\prime$'},{titre}], 'Interpreter','latex')
% xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
% ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')
% %%
% %plot : areas of alpha1 vs alpha2
% figure(2)
% surf2=surf(R1,R2,double(ALPHA1>ALPHA2));
% surf2.EdgeColor = 'none';
% text(xG1,yG1,maxv,'\rho\prime_1<\rho\prime_2','HorizontalAlignment','center','FontSize',20)
% text(xG2,yG2,maxv,'\rho\prime_1>\rho\prime_2','HorizontalAlignment','center','FontSize',20)
% colormap([[1,1,1];[0.8,0.8,0.8]]);
% 
% hold on;
% plot3(R10,yR2,ones(length(R10),1),'-k');
% ylim([1,max(R20)]);
% title([{'zones $\rho_j\prime<\rho_i\prime$ '},titre], 'Interpreter','latex')
% xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
% ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')
% 
% %%
% %plot : beta1 vs beta2
% figure(3)
% surf3 = surf(R1,R2,double(BETA1>BETA2));
% surf3.EdgeColor = 'none';
% hold on;
% plot3(R10,yR2,maxv*ones(length(R10),1),'-k');
% ylim([1,max(R20)]);
% text(xG1,yG1,maxv,'\beta_1<\beta_2','HorizontalAlignment','center','FontSize',20,'Color','w')
% text(xG2,yG2,maxv,'\beta_1>\beta_2','HorizontalAlignment','center','FontSize',20)
% title([{'$\beta_1$ vs $\beta_2$'},titre], 'Interpreter','latex')
% xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
% ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')
% 
% %%
% figure(4)
% %plot :  gradient of change
% %colorbars
% greyred = subdivisedColormap([[0.2,0.2,0.2];[1,1,1];[0.85,0.325,0.098]],n, 'quad'); %2^n+1
% %put 0 as a median of the colorbar
% amp = maxv-minv; %amplitude des valeurs
% %maxv/amplitude %pourcentage de valeurs +
% %minv/amplitude %pourcentage de valeurs -
% greyred2 = greyred(ceil((2^n+1)*0.5*(1-abs(minv)/maxv)):end,:);
% greyyellow = subdivisedColormap([[0.2,0.2,0.2];[1,1,1];[0.929,0.8,0.125]],n, 'quad'); %2^n+1
% greyyellow2 = greyyellow(ceil((2^n+1)*0.5*(1-abs(minv)/maxv)):end,:);
% 
% %surf1=surf(R1,R2,tot)
% surf1=surf(R1,R2,dis2drives1)
% surf1.EdgeColor = 'none';
% colormap(greyred2)
% cb = colorbar;
% cb.YTick = [-0.5 0 0.5 1 1.5];
% cb.YTickLabel = {'-50%', 'No Change', '+50%','+100%', '+150%'};
% freezeColors %adding to the second area where aplha1>alpha2
% title([{'Is combined testing better than specific testing ?'},...
%     {'$(c_{j2} - c_{j1\times2})/c_{j2}$ such that $\rho_j\prime<\rho_i\prime$'},{titre}], 'Interpreter','latex')
% xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
% ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')
% hold on
% %%
% surf2=surf(R1,R2,dis1drives2);
% surf2.EdgeColor = 'none';
% cb2 = colorbar('southoutside');
% colormap(greyyellow2)
% cb2.YTick = [-0.5 0 0.5 1 1.5];
% cb2.YTickLabel = {'-50%', 'No Change', '+50%','+100%', '+150%'};
% plot3(R10,yR2,maxv*ones(length(R10),1),'-k','LineWidth',1.5);
% ylim([1,max(R20)]);
% 
% %%
% figure(5)
% tot2 = (ALPHA1<ALPHA2 & (c11<CEL1))*1 +...
%     (ALPHA1<ALPHA2 & (c11>=CEL1))*2 +...
%     (ALPHA1>=ALPHA2 & (c22<CEL1))*4 +...
%     (ALPHA1>=ALPHA2 & (c22>CEL1))*3
% surf4=surf(R1,R2,tot2);
% surf4.EdgeColor = 'none';
% new_cb2 = [[0.85,0.325,0.098];[0.8,0.8,0.8];[0.8,0.8,0.8];[0.9290,0.6940,0.125]]; %2^n+1
% colormap(new_cb2)
% text(xG1,yG1,4,'\rho\prime_1<\rho\prime_2','HorizontalAlignment','center','FontSize',20)
% text(xG2,yG2,4,'\rho\prime_1>\rho\prime_2','HorizontalAlignment','center','FontSize',20)
% hold on;
% plot3(R10,yR2,4*ones(length(R10),1),'-k','LineWidth',1.5);
% ylim([1,max(R20)]);
% title([{'Is combined testing better than specific testing ?'},...
%     {'$(c_{j2} - c_{j1\times2})/c_{j2}$ such that $\rho_j\prime<\rho_i\prime$'},{titre}], 'Interpreter','latex')
% xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
% ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')
% 

%% Utility septembre 2024
clear all; close all;
%--------------------------
vecC = 0.001:0.0005:0.003
b = 1000; mu = 1/30.6;

gamma1 = 365/14*0.11;
s1 = 1;%0.85; % urine testing ct
rho1b = 0.12;

gamma2 = 0.8*365/14; 
s2 = 1;0.96; % ng
rho2b = 0.12;

P1b = 5.3/100; R1b = 1/(1-P1b);
beta1 = R1b*(gamma1+s1*rho1b+mu);
P2b = 4.7/100; R2b = 1/(1-P2b);
beta2 = R2b*(gamma2+s2*rho2b+mu);

R10 = beta1/(gamma1 + mu);
R20 = beta2/(gamma2 + mu);
alpha1 = beta1*(1-1/R10)/s1;
alpha2 = beta2*(1-1/R20)/s2;
vecRho = sort([alpha1,alpha2,linspace(0,max(alpha1,alpha2)*1.5,200)]); %1:ct, 2:ng
%--------------------------
% syms rho
% gamma2p = gamma2 + s2*rho;
% gamma1p = gamma1 + s1*rho;
% R1p = beta1/(gamma1p + mu);
% R2p = beta2/(gamma2p + mu);
% 
% Nequ = b/mu;
% E0 = [b/mu,0,0,0];
% E1 = [Nequ/R1p, Nequ - Nequ/R1p,0,0];
% E2 = [Nequ/R2p, 0, Nequ - Nequ/R2p,0];
% 
% %Coinfection endemic state 
% gamma1t = gamma1p-s1*s2*rho;
% gamma2t = gamma2p-s1*s2*rho;
% num     = gamma1t/R2p+gamma2t/R1p+mu+s1*s2*rho;
% denum   = beta1+beta2-mu-s1*s2*rho;
% S12     = Nequ*num/denum;
% [S12, Nequ/R2p-S12, Nequ/R1p-S12, S12 + Nequ*(1-1/R1p - 1/R2p)]



i=1;P12=[];P1=[];P2=[];
Nequ = b/mu;
for rho=vecRho
    gamma2p = gamma2 + s2*rho;
	gamma1p = gamma1 + s1*rho;
    R1p = beta1/(gamma1p + mu);
	R2p = beta2/(gamma2p + mu);
    gamma1t = gamma1p-s1*s2*rho;
    gamma2t = gamma2p-s1*s2*rho;
    num     = gamma1t/R2p+gamma2t/R1p+mu+s1*s2*rho;
    denum   = beta1+beta2-mu-s1*s2*rho;
    S12     = Nequ*num/denum;
    P12(i) = max(1-S12/Nequ,0);
    P1(i) = max(1-1/R1p,0);
    P2(i) = max(1-1/R2p,0);

    i=i+1;
end

% figure()
% hold on
% plot(vecRho,P1)
% plot(vecRho,P2)
% plot(vecRho,P12)
%ylim([0,1])
for c=vecC
% utility
U12 = vecRho.*(P12-c);
U1 = vecRho.*(P1-c);
U2 = vecRho.*(P2-c);
figure(1)
hold on
% plot(vecRho,U1)
% plot(vecRho,U2)
% plot(vecRho,U12)
plot(vecRho, max(U12,U2),'DisplayName',['c=',num2str(c)])
end
legend()
%legend('U_1','U_2','U_{12}')

[aa,bb] = max(U12)
vecRho(bb)

tab = [vecRho',P1',P2',P12']

%%
close all

addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\_archives\')
vecC=sort([0.135,0.132,0.133,0.137,-0.051,-0.052,-0.055,-0.25:0.01:0,0.001:0.0001:0.003,0.004:0.001:0.005,0.01:0.01:0.4]);
rho12=[]; rho2=[]; i=1; rhohat =[]; rho1=[]; rho12th=[]; up_bnd=[];
for c=vecC
    CSIS2 = @(rho) -U12_SISSIS(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    rho12(i) = min([max(fmincon(CSIS2,0),0),max(alpha1,alpha2)]);
    r=beta1*beta2/(gamma1*gamma2);
    A = (beta1/gamma2+beta2/gamma1+r+1)/(beta1/gamma2+beta2/gamma1+(2-c)*r);
    rho12th(i) = (beta1+beta2-mu)-sqrt((beta1+beta2)*(beta1+beta2-mu)*A);
    up_bnd(i) = (beta1+beta2-mu)*(1-1/sqrt(2));

    rho2(i)  = min(max(beta2/(2*s2)*(1-1/R20-c),0),alpha2);
    rho1(i)  = min(max(beta1/(2*s1)*(1-1/R10-c),0),alpha1);
    
    U_12 = U_SISSIS(rho12(i),beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    U_2  = U_SISSIS(rho2(i),beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    
    if U_12<=U_2
        rhohat(i) = rho2(i);
    else 
        rhohat(i) = rho12(i);
    end
    
    i=i+1;
end

figure()
plot(vecC,min(rhohat,alpha2))
hold on
%plot(vecC,min(rho2,alpha2))

tab = [vecC',rho12',rho1',rho2',rhohat',rho12th',up_bnd'];
tab=round(tab,5)
format short


%% s1=s2=1
r=beta1*beta2/(gamma1*gamma2);
A = (beta1/gamma2+beta2/gamma1+r+1)/(beta1/gamma2+beta2/gamma1+(2-c)*r);
rho12th = (beta1+beta2-mu)-sqrt((beta1+beta2)*(beta1+beta2-mu)*A);

%c=0
rho12th = (beta1+beta2-mu)-sqrt((beta1+beta2)*(beta1+beta2-mu)*(1-(r-1)/(2*r+beta1/gamma2+beta2/gamma1)));




%% encadrement de rho12 par rho1'/2 et rho1'
close all
s1=1;s2=1; c=0;

syms beta1 gamma1 beta2 gamma2 mu
% 
% mu=1/30.6;
% gamma1 = 365/14*0.11;
% rho1b = 0.12;
% 
% %gamma2 = 19.2473;% 0.8*365/14; 
% rho2b = 0.12;
% 
% P1b = 5.3/100; R1b = 1/(1-P1b);
% beta1 = R1b*(gamma1+s1*rho1b+mu);
% P2b = 4.7/100; R2b = 1/(1-P2b);
% beta2 = R2b*(gamma2+s2*rho2b+mu);
% 

r=beta1*beta2/(gamma1*gamma2);
A = (beta1/gamma2+beta2/gamma1+r+1)/(beta1/gamma2+beta2/gamma1+(2-c)*r);
rho12th = (beta1+beta2-mu)-sqrt((beta1+beta2)*(beta1+beta2-mu)*A);

alpha1=beta1*(1-(gamma1+mu)/beta1)/s1;


eqn = rho12th==alpha1;

s= solve(eqn,gamma2)
simplify(-(beta1*(- beta1^2 - 2*beta1*beta2 + mu*beta1 + beta2^2 + 3*beta2*gamma1 + mu*beta2 + gamma1^2))/(beta2*gamma1 - 2*beta1*beta2 + beta1*mu + beta2*mu - beta1^2))





