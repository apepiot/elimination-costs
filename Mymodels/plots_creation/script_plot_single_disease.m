
%% SICT
clear all; close all;
%betaI=4;betaC=2;sigma=1/12;gamma=6/7;mu=1/35;nu=1/5;b=2;theta=10;
mu=1/30.6;b=5;
P0=0.143; ratioBeta=9.1; sigma=52/8.2;
theta=1/4.4; R0=1/(1-P0); betaC=R0*(sigma+mu)*(theta+mu)./(ratioBeta*(theta+mu)+sigma);
betaI=betaC*ratioBeta;
param.betaI=betaI; param.betaC=betaC;param.sigma=sigma;param.gamma=0;param.theta=theta;
syms rho

[~,~,alpha] = Rp_SICR_v4(param.betaI,param.betaC,param.theta,param.sigma,param.gamma,mu,b,0);
param.alpha=alpha;
[~,c0] = U1_SICR_v4(param,mu,b,0,0,1);
[~,c1] = U1_SICR_v4(param,mu,b,param.alpha,0,1); 
%c1 = c1*1.02 ;%%%probleme dans dU/drho=0, non convexe (voir notes du 26/07),c1*1.01709
options = optimset('Display','off'); %options for minsearch
i=1; 
vecC = [(c1-(c0-c1)/2):0.0005:(c0+(c0-c1)/2)];
%vecC =[(c1-(c0-c1)/2):0.0001:-0.00965,-0.00965:0.000001:-0.00958,-0.00958:0.000001:(c0+(c0-c1)/2)];
%vecC=linspace(-9.65,-9.5,5)*10^(-3)
vecRhomax1 = [];
for c=vecC    
    fun = @(rho) -U_SICR_v4(param,mu,b,rho,c,1);
    vecRhomax1(i) = min(max(fminsearch(fun,0,options),0),param.alpha);
    %vecRho=0:0.001:2;
    %U = U_SICR_v4(param,mu,b,vecRho,c,1);
    %[~,dU] = U1_SICR_v4(param,mu,b,vecRhomax1(i),c,1)
    %plot(vecRho,max(U,0),'DisplayName',num2str(c)); hold on;
    i=i+1; %c
end

%figure('Renderer', 'painters', 'Position', [10 10 1200 400]) %rectangle 1x3
figure('Renderer', 'painters', 'Position', [10 10 750 600])
%set(gcf, 'PaperSize', [4 2]);
subplot(2,2,1);
plot_paper1_procedure4(vecRhomax1, c1, c0, param.alpha, vecC,'',...
    [-0.03,0,0.03],[0,0.02,0.04,0.08],'A')
title('HIV (SICT)','Interpreter','latex')

%% SEIIIS
clear all;
%beta=4;sigma=1/12;gamma10=6/7;gamma30=1;mu=1/35;nu=1/5;b=2;theta=10; tau=1; 
mu=1/30.6;b=5;
P0=0.077;  R0=1/(1-P0); sigma=365/25; tau=365/46;theta=12/3.6;gamma30=1/20;
beta=R0*(sigma+mu)*(theta+mu)*(gamma30+mu)*(tau+mu)./(sigma*((gamma30+mu)*(tau+theta+mu)+tau*theta));
param.beta=beta;param.sigma=sigma;param.gamma1=0;param.gamma3=gamma30;param.theta=theta;param.tau=tau;
param.nu=0;
syms rho
%let's find rho'
%R0 = sigma*beta*(tau+theta+mu)/((theta+mu)*(gamma10+tau+mu)*(sigma+mu));
[~,~,alphaS] = Rp_SEIIIS_v4(param.beta,param.sigma,param.tau,param.nu,param.gamma1,param.theta,param.gamma3,mu,b,0);
param.alpha=alphaS;
[~,c0] = U1_SEIIIS_v4(param,mu,b,0,0,1);
[~,c1] = U1_SEIIIS_v4(param,mu,b,param.alpha,0,1); 

options = optimset('Display','off'); %options for minsearch
i=1; 
vecC =(c1-(c0-c1)/2):(c0-c1)/200:(c0+(c0-c1)/2);
vecRhomax1 = [];
for c=vecC    
    fun = @(rho) -U_SEIIIS_v4(param,mu,b,rho,c,1);
    vecRhomax1(i) = min(max(fminsearch(fun,0,options),0),param.alpha);
    i=i+1; %c
end
subplot(2,2,2); 
plot_paper1_procedure4(vecRhomax1, c1, c0, param.alpha, vecC,'',...
    [-0.12,0,0.12],[0,0.0025,0.005],'B')
%ytickformat('%.4f')
title('Syphilis (SEIIIS)','Interpreter','latex')

%% SEIIS %Ct
clear all;
mu=1/30.6;b=5;
P0=0.074;R0=1/(1-P0);sigma=365/11;nu=1;gamma=365/14;eps1=0.11;
beta=R0*(sigma+mu)*(gamma+nu+mu)*(nu+mu)./(sigma*(gamma*(1-eps1)+nu+mu));
param.beta=beta;param.eps=eps1;param.sigma=sigma;param.gamma=gamma;param.nu=nu;
%syms rho

[~,~,alpha] = Rp_SEIIS_v4(param.beta,param.nu,param.eps,param.sigma,param.gamma,mu,b,0);
param.alpha=alpha;
[~,c0] = U1_SEIISv4(param,mu,b,0,0,1);
[~,c1] = U1_SEIISv4(param,mu,b,param.alpha,0,1); 

options = optimset('Display','off'); %options for minsearch
i=1; 
vecC =(c1-(c0-c1)/2):(c0-c1)/200:(c0+(c0-c1)/2);
vecRhomax1 = [];
for c=vecC    
    fun = @(rho) -U_SEIIS_v4(param,mu,b,rho,c,1);
    vecRhomax1(i) = min(max(fminsearch(fun,0,options),0),param.alpha);
    i=i+1; %c
end

%figure(1)
subplot(2,2,3);
plot_paper1_procedure4(vecRhomax1, c1, c0, param.alpha, vecC,'',...
    [-0.12,0,0.12],[0,0.03,0.06],'C')
title('Ct (SEIIS)','Interpreter','latex')

%% SEIIS %Ng
clear all;
mu=1/30.6;b=5;
P0=0.07;R0=1/(1-P0);sigma=365/5;nu=12/6;gamma=365/14;eps1=0.8;
beta=R0*(sigma+mu)*(gamma+nu+mu)*(nu+mu)./(sigma*(gamma*(1-eps1)+nu+mu));
param.beta=beta;param.eps=eps1;param.sigma=sigma;param.gamma=gamma;param.nu=nu;
%syms rho

[~,~,alpha] = Rp_SEIIS_v4(param.beta,param.nu,param.eps,param.sigma,param.gamma,mu,b,0);
param.alpha=alpha;
[~,c0] = U1_SEIISv4(param,mu,b,0,0,1);
[~,c1] = U1_SEIISv4(param,mu,b,param.alpha,0,1); 

options = optimset('Display','off'); %options for minsearch
i=1; 
vecC =(c1-(c0-c1)/2):(c0-c1)/200:(c0+(c0-c1)/2);
vecRhomax1 = [];
for c=vecC    
    fun = @(rho) -U_SEIIS_v4(param,mu,b,rho,c,1);
    vecRhomax1(i) = min(max(fminsearch(fun,0,options),0),param.alpha);
    i=i+1; %c
end

%figure(1)
subplot(2,2,4);
plot_paper1_procedure4(vecRhomax1, c1, c0, param.alpha, vecC,'',...
    [-0.1,0,0.1],[0,0.075,0.15],'D');
title('Ng (SEIIS)','Interpreter','latex')
%%
set(gcf, 'PaperSize', [31 30]);
saveas(gcf,'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v3\graphes\subplot_1dis_2.pdf')
saveas(gcf,'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v3\graphes\subplot_1dis_2.png')

