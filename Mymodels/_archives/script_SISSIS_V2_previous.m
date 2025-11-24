%% Script mymodel SISxSIS (V2)
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
beta2 = 0.4;
gamma1 = 0.02;
gamma2 = 0.02;
b = 0.05;
mu = 0.01;
s = 0.9;
delta = 0.2;
e1 = beta1;e2 = beta2;


% Parametres du systeme d'ODE 
tspan = 0:1:5000;
Y0 = [S0; I10; I20; I120];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) SISxSIS_V2(t,Y,b,beta1, beta2,e1,e2,gamma1,gamma2,s,delta, mu, 'frequency'),tspan,Y0, options);

Ys(end,:)
sum(Ys(end,:))

% Plot
figure(1)
plot(ts,Ys)%./(sum(Ys')'));
title('Evolution');
xlabel('Time t');
ylabel('Y');




%% Frequency-dependent transmission case equilibria version 2
R10 = beta1/(gamma1 + s*delta + mu);
R20 = beta2/(gamma2 + s*delta + mu);

Nequ = b/mu;
E0 = [b/mu,0,0,0];
E1 = [Nequ/R10, Nequ - Nequ/R10,0,0];
E2 = [Nequ/R20, 0, Nequ - Nequ/R20,0];

S12 = ((gamma1)/R20 + (gamma2)/R10 + mu + s*delta)/(beta1+beta2-mu-s*delta)*b/mu
E12 = [S12, Nequ/R20-S12, Nequ/R10 - S12, S12 + Nequ*(1-1/R10 - 1/R20)]

Ys(end,:)

%% Utility function
vecDelta = 0:0.005:3.;

%vecGamma1 = gamma1 + s*vecDelta;
%vecGamma2 = gamma2 + s*vecDelta;

[U,vecPrev,vecR10,vecR20] = utilityFunctionSIS2_V2(beta1,beta2,gamma1,gamma2,s, vecDelta, b,mu,beta1,beta2);

%U = vecDelta.*vecPrev;

figure(10)
plot(vecDelta, vecPrev)
title([{'Prevalence of the SISxSIS model'},...
    {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
    num2str(gamma2), ' s=', num2str(s), '\mu=', num2str(mu), ' \pi=',  num2str(b), ' R_0^1=' num2str(round(beta1/(gamma1+mu),2)), ' R_0^2=' num2str(round(beta2/(gamma2+mu),2)) ]}])
xlabel("Voluntary-testing rate \delta","fontweight","bold")
ylabel("I(\delta)","fontweight","bold")


figure(11)
plot(vecDelta, U)
title([{'Utility function of the SISxSIS model'},...
    {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
    num2str(gamma2), ' s=', num2str(s), '\mu=', num2str(mu), ' \pi=',  num2str(b), ' R_0^1=' num2str(round(beta1/(gamma1+mu),2)), ' R_0^2=' num2str(round(beta2/(gamma2+mu),2)) ]}])
xlabel("Voluntary-testing rate \delta","fontweight","bold")
ylabel("U(\delta)","fontweight","bold")
hold on;


[Umax, imax] = max(U)
deltaNum = vecDelta(imax)

%sqrt((gamma1/beta2 + gamma2/beta1+1)*(beta1+beta2-mu)*(gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1+mu))/(gamma1/beta2+gamma2/beta1+1)
%sqrt((beta1+beta2-mu)*(gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1 +mu)/(gamma1/beta2 + gamma2/beta1 +mu) )/s^2


%% resultats de la v1
% 
% A = 3/beta2+3/beta1;
% B = 2*(gamma1+gamma2+mu)*(1/beta1+1/beta2);
% C = -beta1 - beta2 + 2*mu + gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1;
% 
% DISC = B^2 - 4*A*C;
% 
% DiscTh = 4*(gamma1+gamma2+mu)^2*(1/beta1+1/beta2)^2 - 12*(1/beta1 + 1/beta2)*(-beta1-beta2+2*mu + gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1);
% 
% deltaTh = (-B + sqrt(DISC))/(2*A)/s
% 
% truc = 1 - 3*beta1*beta2* (-beta1-beta2+2*mu + gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1 )/(gamma1 + gamma2 + mu)^2/(beta1+beta2)
% deltaTh2 = (gamma1 + gamma2 + mu)/3*(-1 + sqrt(truc))/s

%% resultats provisoires v2
syms gamma1 gamma2 beta1 beta2 s delta mu 
A = (gamma1/beta2+gamma2/beta1 +2);
%B = -3*(beta1+beta2-mu) + 2*(gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1 + mu)...
%    - (gamma1/beta2+gamma2/beta1)*(beta1+beta2-mu) - (gamma1/beta2 + gamma2/beta1+1);
%B = -4*(beta1+beta2-mu) - 2*(beta1+beta2-mu)*(gamma1/beta2+gamma2/beta1);
B = 2*(mu - beta1 - beta2)*(gamma1/beta2 + gamma2/beta1 + 2);
%C = (beta1 + beta2-mu)^2 - (beta1+beta2-mu)*(gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1 + mu);
C = (beta1 + beta2-mu)*(beta1 + beta2-mu - (gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1 + mu))

DELTA = B^2 - 4*A*C;
deltaTh = (-B - sqrt(DELTA))/(2*A)/s %coifection
deltaTh = (-B + sqrt(DELTA))/(2*A)/s

syms x
eqn = (beta1+beta2-mu-x)^2 - (beta1+beta2-mu)*(gamma1*(gamma2 + x+mu)/beta2 + gamma2*(gamma1+x+mu)/beta1 + mu +x)...
    - x*(gamma1/beta2 + gamma2/beta1 +1)*(beta1+beta2-mu-x) + 2*x*(gamma1*(gamma2 + x+mu)/beta2 + gamma2*(gamma1+x+mu)/beta1 + mu +x)==0 
solve(eqn)

eqn = x^2*A + x*B + C ==0
solve(eqn)


deltaTh1 = (beta1 - mu - gamma1)/2/s
deltaTh2 = (beta2 - mu - gamma2)/2/s

syms x

R2 = beta2/(gamma2+x+mu);
R1 = beta1/(gamma1+x+mu);

denom = beta1+beta2-mu - x

eqn = 1 - (gamma1/R2+gamma2/R1+mu+x)/denom - x*( (gamma1/beta2+gamma2/beta1+1)/denom + (gamma1/R2+gamma2/R1+mu+x)/denom^2) ==0
sol = solve(eqn)/s
sol = 0.3294

R2 = beta2/(gamma2+sol+mu);
R1 = beta1/(gamma1+sol+mu);
sol/s* (1 - (gamma1/R2+ gamma2/R1 +mu+sol)/(beta1+beta2-mu-sol)):


