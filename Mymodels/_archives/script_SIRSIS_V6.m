%% Script du modele SIRxSIS version 1 bis
%ce qui change ici, c'est le taux de transmission d'une maladie d'un
%individu coinfecte : epsilon1,2
clear all
close all
%% Parametres
% Populations initiales
S0      = 50;
I10     = 5;
I20     = 3;
I30     = 2;    %I_12 (coinfection)
I40     = 0;
R10     = 0.;
Y0 = [S0; I10; I20; I30; I40; R10];

% Definition des parametres 
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(false,false)
e1 = myrand(1,0,1,2)*beta1;
e2 = myrand(1,0,1,2)*beta2;

gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;

% Parametres du systeme d'ODE 
tspan = 0:1:1000;
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_V5(t,Y,b,beta1,beta2,e1,e2,s1,s2,rho,gamma1,gamma2,mu,'frequency'),tspan,Y0, options);

T=Ys(end,:)


%% Reproduction numbers (doivent tenir compte ded epsilon ?)
syms beta1 beta2 e1 e2 s1 s2 rho gamma1 gamma2 mu
gamma1p=gamma1;gamma2p=gamma2;
I1=0;I2=0;I12=0;IR2=0;R1=0; S=1;

lambda1 = beta1*I1+e1*I12;
lambda2 = beta2*(I2+IR2)+e2*I12;

F = [beta1*S, 0, e1*S,0;...
    0, beta2*S, e2*S, beta2*S;...
    beta1*I2+lambda2, lambda1+beta2*I1, e1*I2+e2*I1,beta2*I2;...
    0,beta2*R1,e2*R1,beta2*R1];
V = [lambda2+mu+gamma1p,beta2*I1,-e2*I1+gamma2,beta2*I1;...
    beta1*I2,lambda1+mu+gamma2p,e1*I2,0;...
    0,0,gamma1p+gamma2p+mu,0;...
    0,0,-gamma1p,gamma2p+mu];

eig(F*V^(-1))
R1p = beta1/(gamma1p+mu);%
R2p = beta2/(gamma2p+mu);%

%% Finding the ES 
%modifier ODE_V5 : c'est pas les bonnes equations
%dS = b - (lambda1+lambda2)*S+gamma2*
