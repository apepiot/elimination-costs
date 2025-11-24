%% Solving the utility
clear all; close all;

mu  = 1/30.6; b=10000;

%SEIIS1 %chlamydia/gono
% P1VT    = 5.3/100;   
% R1VT    = 1/(1-P1VT);          
% sigma1  = 365./11;      
% nu1     = 12/12;      
% gamma1  = 365/14;
% eps1    = 0.11; 
% rhob1   = 1/8.14;
% beta1   = R1VT/((sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rhob1))./((mu + sigma1 + rhob1).*(gamma1 + mu + nu1).*(mu + nu1 + rhob1))); %check
% [R10,~,alpha1] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,0);
% P10     = 1-1/R10;
%SEIIS2 %gono
P1VT     = 4.7/100;
R1VT     = 1/(1-P1VT);
sigma1   = 365/5;
nu1      = 12/6;
gamma1   = 365/14;
eps1     = 0.8;
rhob1    = 1/8.14;
beta1    = R1VT/((sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rhob1))./((mu + sigma1 + rhob1).*(gamma1 + mu + nu1).*(mu + nu1 + rhob1)));
[R10,~,alpha1] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,0); %check
P10      = 1-1/R10;

param.beta=beta1; param.sigma=sigma1; param.gamma=gamma1;
param.nu=nu1;param.eps=eps1;
param.alpha=alpha1;
param.R_base=R1VT;param.R0=R10;
param.rhob = rhob1;
param.modelType = 'SEIIS'; param.mini_d = 'c';
param.P_base = P1VT; param.P0 = P10;
param.mu = mu; param.pi=b;

vecC = linspace(-0.5,0.5,500);
syms rho
i=0; rhohat=[]; Rrho=[]; Ptot=[];
for c=vecC
    i=i+1;
    C_c = @(rho)-U1_SEIISv4(param,mu,b,rho,c);
    rhohat(i) = fmincon(C_c,0.,[-1;1],[0;alpha1]);
    [Rrho(i),~,~,Ptot(i)] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,rhohat(i));
end


figure(1)
hold on
plot(vecC,rhohat)
xlabel('$c$','Interpreter','latex')
ylabel('$\hat\rho$','Interpreter','latex')

tabrecap = [vecC;rhohat;Rrho;Ptot]';

tabrecap(tabrecap(:,2)<0.0001,2)=0;
tabrecap(tabrecap(:,4)<0.0001,4)=0;
