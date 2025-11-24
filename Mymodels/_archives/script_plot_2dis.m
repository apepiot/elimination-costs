%RHOHAT_v3_3 for PLOTTING 

% one SICR : HIV
% one SEIIIS : syphilis
% one or two SEIIS : chlamydia, gonorrhea

%% Parameters
clear all;%close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ct=true; Ng=false; HIV=true; Syph=false; plotFigure=true;
biasFactor = 1;
%%
nbSim = 1;

nSEIIS=Ct+Ng; nSICR=HIV; nSEIIIS=Syph;
N = nSICR+nSEIIS+nSEIIIS;
%STIChoice = 'Ct'; STIChoiceMini = 'Ct';%if one STI considered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = 5;

bestStrat={};worstStrat={};tps1=tic;tabRecapStratCurrentSim=[];

if(0)
mu = 0.0335;
%SEIIS1 %chlamydia
R1     = 1.0860;
sigma1 = 40.2979;
nu1    = 0.5246;
gamma1 = 24.8922; 
eps1   = 0.0486;
beta1  = 0.6346;
[~,~,alpha1] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,0);

%SEIIS2 %gono
P2     = randPERT(6,7,8,1)/100;
R2     = 1/(1-P2);
sigma2 = 365/randPERT(1,5,14,1);
nu2    = 12/randPERT(4,6,12,1);
eps2   = randPERT(0.5,0.8,1,1);
gamma2 = 365/randPERT(5,14,20,1);%365/3; 
beta2  = R2/((sigma2*(gamma2*(1-eps2) + mu + nu2))./((mu+sigma2).*(gamma2 + mu + nu2).*(mu + nu2)));
[~,~,alpha2] = Rp_SEIIS_v4(beta2,nu2,eps2,sigma2,gamma2,mu,b,0);

%SICR 
RSICR     = 1.1914;
sigmaSICR = 5.6398;
thetaSICR = 0.1986;
gammaSICR = 0; %0.39*365/(8.2*7); %0.4*(); %proportion de diag au stade precoce : 40
%ratioBeta = randPERT(8.4,9.1,9.6,1);
betaCSICR = 0.2000;
betaISICR = 1.8551;
[~,~,alphaSICR] = Rp_SICR_v4(betaISICR,betaCSICR,thetaSICR,sigmaSICR,gammaSICR,mu,b,0);

%SEIIIS/S(syphilis)
PS      = randPERT(6.7,7.7,8.7,1)/100;
RS      = 1/(1-PS);
sigmaS  = 365/randPERT(10,25,95,1);
tauS    = 365/randPERT(10,46,130,1);%(1-0.55*0.31)*365/45; %0.6*
thetaS  = 12/randPERT(1,3.6,12,1);
gamma1S = 0;%(0.55*0.31)*365/45; %0.2*
gamma3S = 1/randPERT(10,20,30,1);%1/5;
nuS     = 0;
betaS = RS./(sigmaS.*((gamma3S+mu+nuS).*(mu+tauS+thetaS)+tauS*thetaS)./((mu + sigmaS).*(mu + thetaS).*(gamma3S + mu + nuS).*(gamma1S + mu + tauS)));
[~,~,alphaS] = Rp_SEIIIS_v4(betaS,sigmaS,tauS,nuS,gamma1S,thetaS,gamma3S,mu,b,0);
end

% choose the set
run("set4.m")


%paramSEIIS = [beta1,gamma1,nu1,sigma1,eps1,alpha1,R1;...
%              beta2,gamma2,nu2,sigma2,eps2,alpha2,R2];
paramCt.beta = beta1;paramCt.sigma=sigma1;paramCt.gamma=gamma1;paramCt.nu=nu1;paramCt.eps=eps1;paramCt.alpha=alpha1;paramCt.R=R1;
paramNg.beta = beta2;paramNg.sigma=sigma2;paramNg.gamma=gamma2;paramNg.nu=nu2;paramNg.eps=eps2;paramNg.alpha=alpha2;paramNg.R=R2;
paramCt.disease = 'Ct'; paramCt.modelType = 'SEIIS';
paramNg.disease = 'Ng'; paramNg.modelType = 'SEIIS';
%paramSEIIS = {paramCt,paramNg};

%paramSICR = [betaISICR,betaCSICR,gammaSICR,sigmaSICR,thetaSICR,alphaSICR,RSICR,etaH,omegaHIV];
paramHIV.betaI = betaISICR; paramHIV.betaC = betaCSICR; paramHIV.gamma = gammaSICR; paramHIV.theta = thetaSICR;...
    paramHIV.sigma=sigmaSICR; paramHIV.alpha = alphaSICR;
paramHIV.R = RSICR; %paramHIV.eta = etaH; paramHIV.omega=omegaHIV;
paramHIV.disease = 'HIV'; paramHIV.modelType = 'SICR';

%paramS    = [betaS,sigmaS,tauS,gamma1S,thetaS,gamma3S,nuS,alphaS,RS,etaS,omegaS];
paramS.beta = betaS; paramS.sigma=sigmaS; paramS.tau=tauS; paramS.gamma1=gamma1S; paramS.theta=thetaS; paramS.gamma3=gamma3S;
paramS.nu=nuS;paramS.alpha=alphaS;%paramS.eta=etaS;paramS.omega=omegaS;
paramS.R=RS;
paramS.disease = 'Syph.'; paramS.modelType = 'SEIIIS';

vecAlphas = [paramCt.alpha(Ct),paramNg.alpha(Ng),paramHIV.alpha(HIV),paramS.alpha(Syph)];
maxalpha = max(vecAlphas);
paramTabAll{1}=paramCt;
paramTabAll{2}=paramNg;
paramTabAll{3}=paramHIV;
paramTabAll{4}=paramS;

paramTab = paramTabAll([Ct,Ng,HIV,Syph]);
%diseases.Ct=Ct;diseases.Ng=Ng;diseases.HIV=HIV;diseases.S=Syph;

%% PLOT RHOHAT
cleft=-0.2; cright=0.25;
vecC = linspace(cleft,cright,1000); 
%vecC =[linspace(cleft,-0.0155,200),linspace(-0.015,-0.014,20),linspace(-0.0135,cright,100)]; %set7
[tab,tabco,tabcn,tabTimes] = findRhohat_v5(N,paramTab,mu,b,vecC,1);
cs = findThresholds_v4(nSEIIS,nSICR,nSEIIIS, tab, vecAlphas, vecC);
plot_zones_v3ter



%yticks([0,0.05,0.1]) %set7
%yticks([0,0.06,0.12]) %set6
yticks([0,0.025,0.05,0.075]) %set4
% change to paper 2.3
set(gcf, 'PaperSize', [15 15]);
saveas(gcf,'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v3\graphes\SEIISxSICT_set4.pdf')
saveas(gcf,'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v3\graphes\SEIISxSICT_set4.png')


%%
close all
rhohat=[];i=1;
vecRho=0:0.001:0.1
vecC=[-0.01:0.0005:0];%-0.004:0.0001:0;
vecC=-0.0044
c=0
for c=vecC
    [U,dU]=U1_SICR_v4(paramHIV, mu, b, vecRho,c,1);
    hold on
    plot(vecRho,U)
    fun = @(rho) -U1_SICR_v4(paramHIV,mu,b,rho,c,1);
    rhohat(i) = min(max(fminsearch(fun,0),0),paramHIV.alpha);
    i=i+1;
end

figure()
plot(vecC,rhohat)
[~,cnn] = U1_SICR_v4(paramHIV,mu,b,paramHIV.alpha,0,1)



%% plot of the zones
clear all; close all;
col = [132/255, 151/255, 176/255];
set(groot,'defaultAxesTickLabelInterpreter','latex');  

vecC=[-1,1];
c0=0.5; c1=0; c2=-0.5;
limy=2; ex=0.1; ey=0.1;

%Areas
a3=area([vecC(1) vecC(1) c2 c2],[0,limy,limy,0 ],...
    'LineStyle','None','DisplayName',['both infections eliminated']); a3(1).FaceColor = col; a3.FaceAlpha = 0.5;
text((vecC(1)+c2)/2,limy*0.8,'\fontsize{22}{0} \selectfont III','Interpreter','latex','HorizontalAlignment', 'center')
hold on
a2b=area([c2 c2 c1 c1],[0,limy,limy,0 ],...
    'LineStyle','None','DisplayName',['one disease eliminated']); a2b(1).FaceColor = col; a2b.FaceAlpha = 0.2;
text((c2+c1)/2,limy*0.8,'\fontsize{22}{0} \selectfont IIb','Interpreter','latex','HorizontalAlignment', 'center')
a2a=area([c1 c1 c0 c0],[0,limy,limy,0 ],...
    'LineStyle','None','HandleVisibility','off'); a2a(1).FaceColor = col; a2a.FaceAlpha = 0.;
text((c1+c0)/2,limy*0.8,'\fontsize{22}{0} \selectfont IIa','Interpreter','latex','HorizontalAlignment', 'center')
a1=area([c1 c1 vecC(end) vecC(end)],[0,limy,limy,0 ],...
    'LineStyle','None','DisplayName',['no disease eliminated']); a1(1).FaceColor = col; a1.FaceAlpha = 0.;
text((c0+vecC(end))/2,limy*0.8,'\fontsize{22}{0} \selectfont I','Interpreter','latex','HorizontalAlignment', 'center')

%legend('Interpreter','latex')

%Vertical lines
plot([c2,c2],[0,limy],'k:','LineWidth',1.,'HandleVisibility','off')
plot([c1,c1],[0,limy],'k:','LineWidth',1.,'HandleVisibility','off')
plot([c0,c0],[0,limy],'k:','LineWidth',1.,'HandleVisibility','off')

%Labels
xticks([c2,c1,c0])
xticklabels({'$c^2=c^\prime_i$','$c^1$','$c^0$'})
%xlabel('$c$','Interpreter','latex')
text('position',[vecC(end) 0-ey 0],'interpreter','latex',...
    'string','\fontsize{22}{0}\selectfont$c$','HorizontalAlignment', 'center');

yticks([])
text('position',[vecC(1)-1.4*ex limy 0],'interpreter','latex',...
    'string','\fontsize{22}{0}\selectfont$\hat\rho(c)$','HorizontalAlignment', 'center');

ax = gca;
ax.FontSize = 22; 

set(gcf, 'PaperSize', [15 15]);
saveas(gcf,'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v3\graphes\areas_2dis.pdf')
saveas(gcf,'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v3\graphes\areas_2dis.png')




%% don't remember what i wanted to do below
close all;
p = plot(vecRho,U123,'-','Color','black','DisplayName', 'U', 'Linewidth',2);
p123 = plot(vecRho(vecRho<paramTab{3}.alpha),U123(vecRho<paramTab{3}.alpha),'-','Color',[0.4660 0.6740 0.1880],'DisplayName', 'U', 'Linewidth',2); 
hold on
p12 = plot(vecRho(vecRho<=paramTab{3}.alpha),U12(vecRho<=paramTab{3}.alpha),'-','Color',[0 0.4470 0.7410], 'DisplayName', 'U_{c\timesh}','Linewidth',2);
p2 = plot(vecRho(vecRho>=paramTab{3}.alpha),U2(vecRho>=paramTab{3}.alpha),'-','Color',[0.9290 0.6940 0.1250], 'DisplayName','U_c','Linewidth',2);

plot(vecRho(vecRho<=paramTab{3}.alpha),U2(vecRho<=paramTab{3}.alpha),'--','Color',[0.9290 0.6940 0.1250],'Linewidth',2)
plot(vecRho,U1, '--','Color',[0.8500 0.3250 0.0980],'Linewidth',2)
p1 = plot(vecRho(1),U1(1), '-','Color',[0.8500 0.3250 0.0980],'DisplayName','U_h','Linewidth',2);
legend([p,p12,p2,p1])
ylim([0, 1.25*max(U)])

