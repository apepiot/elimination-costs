%% Calcul des intervalles de beta
clear all; close all;
% middle
mu = 1/30.6;

%SEIIS1 %chlamydia
P1     = 7.4/100;
R1     = 1/(1-P1);
sigma1 = 365./11;
nu1    = 12/12;
gamma1 = 365/14;
eps1   = 0.11;
rhob1  = 1/8.14; %baseline voluntary testing rate
beta1  = R1/((sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rhob1))./((mu + sigma1 + rhob1).*(gamma1 + mu + nu1).*(mu + nu1 + rhob1)));

%SEIIS2 %gono
P2     = 7/100;
R2     = 1/(1-P2);
sigma2 = 365/5;
nu2    = 12/6;
gamma2 = 365/14;%365/3;
eps2   = 0.8;
rhob2  = 1/8.14; %baseline voluntary testing rate
beta2  = R2/((sigma2*(gamma2*(1-eps2) + mu + nu2 + eps2*rhob2))./((mu + sigma2 + rhob2).*(gamma2 + mu + nu2).*(mu + nu2 + rhob2)));

%HIV
%PHIV = randPERT(10,20,30,1)/100; %prevalence and not undiagnosed prevalence
LambdaSICR = 7/100; %incidence
%RSICR     = 1/(1-PHIV);%2.2; %P=1-1/R so R=1/(1-P);
sigmaSICR = 52/8.2; %(1-0.39)*365/(8.2*7)%0.6*(); %8.2 semaines
thetaSICR = 1/4.4;%1/9.8;
gammaSICR = 0; %0.39*365/(8.2*7); %0.4*(); %proportion de diag au stade precoce : 40
ratioBeta = 9.1; %betaI/betaC
rhobh     = 1/1.9; %baseline voluntary testing rate
betaCSICR = (LambdaSICR+mu)./mu.*(sigmaSICR+rhobh+mu)*(thetaSICR+rhobh+mu)./(ratioBeta*(thetaSICR+rhobh+mu)+sigmaSICR);
betaISICR = ratioBeta*betaCSICR;

%SEIIIS/S(syphilis)
PS      = 7.7/100;
RS      = 1/(1-PS);
sigmaS  = 365/25;
tauS    = 365/46;%(1-0.55*0.31)*365/45; %0.6*
thetaS  = 12/3.6;
gamma1S = 0;%(0.55*0.31)*365/45; %0.2*
gamma3S = 1/20;
nuS     = 0;
rhobs   = 1/2.48; %baseline voluntary testing rate
betaS = RS/((sigmaS.*((gamma3S+mu+nuS+rhobs).*(mu+rhobs+tauS+thetaS)+tauS*thetaS))./((mu+rhobs+sigmaS).*(mu+rhobs+thetaS).*(gamma3S+mu+nuS+rhobs).*(gamma1S+mu+rhobs+tauS)));

%% intervals
clear all;
pHIV=0;
N=2000;
vecBetaC = zeros(N,1);vecBetaI=vecBetaC;
vecBetaCt= vecBetaC; vecBetaNg=vecBetaC;
vecBetaS = vecBetaC; vecR_h=vecBetaC;
vecR_s=vecBetaC; vecR_c=vecBetaC; vecR_g=vecBetaC;
for i=1:N
    [paramTab,mu,~] = sampleParameters_v3_extent(true,true,true,true,2,pHIV);
    vecBetaI(i) = paramTab{3}.betaI;
    vecBetaC(i) = paramTab{3}.betaC;
    vecBetaCt(i) = paramTab{1}.beta;
    vecBetaNg(i) = paramTab{2}.beta;
    vecBetaS(i) = paramTab{4}.beta;
    
    vecR_h(i)   = paramTab{3}.R_prev_base; %R(\rho_h^0) with pHIV=0;
    vecR_s(i)   = paramTab{4}.R_base;
    vecR_c(i)   = paramTab{1}.R_base;
    vecR_g(i)   = paramTab{2}.R_base;
end

[prctile(vecBetaI,2.5),mean(vecBetaI),prctile(vecBetaI,97.5)]
[prctile(vecBetaC,2.5),mean(vecBetaC),prctile(vecBetaC,97.5)]
[prctile(vecBetaS,2.5), mean(vecBetaS),prctile(vecBetaS,97.5)]
[prctile(vecBetaCt,2.5), mean(vecBetaCt),prctile(vecBetaCt,97.5)]
[prctile(vecBetaNg,2.5), mean(vecBetaNg),prctile(vecBetaNg,97.5)]


%% Summary table of parameters
clear all; close all;
%----------------%
roundNos = 1:2;
folderRun = 'Run5';
path  = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\'];
%----------------%
PCt = readtable([path,'_round_',num2str(roundNos(1)),'\allParametersTabCt'],'ReadVariableNames', true);
PNg = readtable([path,'_round_',num2str(roundNos(1)),'\allParametersTabNg'],'ReadVariableNames', true);
PHIV= readtable([path,'_round_',num2str(roundNos(1)),'\allParametersTabHIV'],'ReadVariableNames', true);
PS  = readtable([path,'_round_',num2str(roundNos(1)),'\allParametersTabS'],'ReadVariableNames', true);

for roundNo=roundNos(2:end)
    PCt_new = readtable([path,'_round_',num2str(roundNo),'\allParametersTabCt'],'ReadVariableNames', true);
    PNg_new = readtable([path,'_round_',num2str(roundNo),'\allParametersTabNg'],'ReadVariableNames', true);
    PHIV_new= readtable([path,'_round_',num2str(roundNo),'\allParametersTabHIV'],'ReadVariableNames', true);
    PS_new  = readtable([path,'_round_',num2str(roundNo),'\allParametersTabS'],'ReadVariableNames', true);
    PCt  = [PCt;PCt_new];
    PNg  = [PNg;PNg_new];
    PHIV = [PHIV;PHIV_new];
    PS   = [PS;PS_new];
end
writetable(PCt,[path,'\_results','\allParametersTabCt','_concatenated.txt'])
writetable(PNg,[path,'\_results','\allParametersTabNg','_concatenated.txt'])
writetable(PHIV,[path,'\_results','\allParametersTabHIV','_concatenated.txt'])
writetable(PS,[path,'\_results','\allParametersTabS','_concatenated.txt'])



%%
close all;
%rho' of each disease (HIV:p=0)
fig = figure(1)
BinWidth=0.1%(max(alphaTot(:))-min(alphaTot(:)))/60;
histogram(PHIV.alpha,'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(PS.alpha,'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(PCt.alpha,'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(PNg.alpha,'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
legend('$\rho_{h}^\prime(p=0)$',...
    '$\rho_{s}^\prime$',...
    '$\rho_{c}^\prime$',...
    '$\rho_{g}^\prime$',...
    'Interpreter','latex','FontSize',12,'Box','off')
%%
path = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\graphes\parametersAnalysis2\';
saveas(fig,[path,'hist_alpha_p_0.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[path,'hist_alpha_p_0.pdf'])

%%
fig = figure(2)
BinWidth=0.05%(max(alphaTot(:))-min(alphaTot(:)))/60;
histogram(PHIV.RSICT,'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(PS.R,'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(PCt.R,'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(PNg.R,'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
legend('$\mathtt R_h(\rho_0^h))$',...
    '$\mathtt R_s(\rho_0^s))$',...
    '$\mathtt R_c(\rho_0^c))$',...
    '$\mathtt R_g(\rho_0^g))$',...
    'Interpreter','latex','FontSize',12,'Box','off')
%%
path = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\graphes\parametersAnalysis2\';
saveas(fig,[path,'hist_R_p_0.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[path,'hist_R_p_0.pdf'])
%%
fig = figure(3)
BinWidth=0.01%(max(alphaTot(:))-min(alphaTot(:)))/60;
histogram(PHIV.rhob,'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(PS.rhob,'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(PCt.rhob,'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(PNg.rhob,'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
legend('$\rho_0^h$',...
    '$\rho_0^s$',...
    '$\rho_0^c$',...
    '$\rho_0^g$',...
    'Interpreter','latex','FontSize',12,'Box','off')

path = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\graphes\parametersAnalysis2\';
saveas(fig,[path,'hist_rhos_0_p_0.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[path,'hist_rhos_0_p_0.pdf'])

%% Comparison of R and alpha for different values of p (PrEP coverage)
clear all; close all;
b = 2;
N = 1000;
vecP = [0,0.25,0.5,0.75];
R = zeros(N,length(vecP)); alpha = zeros(N,length(vecP));
k = 0;
zeta = 0.6; eta = 4;
for p=vecP
    k=k+1;
    for i=1:N
        [paramTab,mu,~] = sampleParameters_v3(false,false,true,false,b);
        paramHIV = paramTab{1};
        [R(i,k),~,alpha(i,k)] = Rp_SICTP(paramHIV.betaI,paramHIV.betaC,...
            paramHIV.theta,paramHIV.sigma,zeta,eta,p,mu,b,paramHIV.rhob);
    end    
end

fig = figure(4)
BinWidth=0.1;
histogram(R(:,1),'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth)
hold on
histogram(R(:,2),'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth)
histogram(R(:,3),'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth)
histogram(R(:,4),'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth)
plot([1 1],[0,0.3],'r-')
legend(['p=',num2str(vecP(1))],...
       ['p=',num2str(vecP(2))],...
       ['p=',num2str(vecP(3))],...
       ['p=',num2str(vecP(4))])
xlabel('$\mathtt R(\rho_0^h)$ with PrEP coverage $p$','Interpreter','latex')

path = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\graphes\parametersAnalysis2\';
saveas(fig,[path,'hist_R_rho0_p_diff_sictp.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[path,'hist_R_rho0_p_diff_sictp.pdf'])

fig =figure(5)
BinWidth=0.2;
histogram(alpha(:,1),'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth)
hold on
histogram(alpha(:,2),'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth)
histogram(alpha(:,3),'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth)
histogram(alpha(:,4),'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth)
legend(['p=',num2str(vecP(1))],...
    ['p=',num2str(vecP(2))],...
    ['p=',num2str(vecP(3))],...
    ['p=',num2str(vecP(4))])
xlabel('$\rho^\prime_h$','Interpreter','latex')

path = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\graphes\parametersAnalysis2\';
saveas(fig,[path,'hist_alpha_p_diff_sictp.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[path,'hist_alpha_p_diff_sictp.pdf'])



