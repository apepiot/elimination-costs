%RHOHAT_v3_3 for PLOTTING 

% one SICR : HIV
% one SEIIIS : syphilis
% one or two SEIIS : chlamydia, gonorrhea

%% Parameters
clear all;%close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ct=true; Ng=true; HIV=true; Syph=true; plotFigure=false;
biasFactor = 1;
%%
nbSim = 1;

nSEIIS=Ct+Ng; nSICR=HIV; nSEIIIS=Syph;
N = nSICR+nSEIIS+nSEIIIS;
%STIChoice = 'Ct'; STIChoiceMini = 'Ct';%if one STI considered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = 5;

%%
bestStrat={};worstStrat={};tps1=tic;tabRecapStratCurrentSim=[];
mu = 1/randPERT(27.2,30.6,33.7,1);
%SEIIS1 %chlamydia
P1     = randPERT(6.4,7.4,8.4,1)/100;
R1     = 1/(1-P1);
sigma1 = 365./randPERT(7,11,21,1);
nu1    = 12/randPERT(6,12,36,1);
gamma1 = 365/randPERT(5,14,20,1);%365/(31.5-1./sigma1); 
eps1   = randPERT(0,0.11,0.5,1);
beta1  = R1/((sigma1*(gamma1*(1-eps1) + mu + nu1))./((mu+sigma1).*(gamma1 + mu + nu1).*(mu + nu1)));
[~,~,alpha1] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,0);
beta_mean = 1/(1-7.4/100)*(365/11+1/30.6)*(365/14+12/12+1/30.6)*(12/12+1/30.6)/(365/11*(365/14*(1-0.11)+1/30.6+12/12));

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
PHIV = randPERT(12.0,14.3,16.9,1)/100;
RSICR     = 1/(1-PHIV);%2.2; %P=1-1/R so R=1/(1-P);
sigmaSICR = 52/randPERT(6.7,8.2,9.8,1); %(1-0.39)*365/(8.2*7)%0.6*(); %8.2 semaines
thetaSICR = 1/randPERT(4,4.4,10,1);%1/9.8;
gammaSICR = 0; %0.39*365/(8.2*7); %0.4*(); %proportion de diag au stade precoce : 40
ratioBeta = randPERT(8.4,9.1,9.6,1);
betaCSICR = RSICR*(sigmaSICR+gammaSICR+mu)*(thetaSICR+mu)/(ratioBeta*(thetaSICR+mu)+sigmaSICR);
betaISICR = ratioBeta*betaCSICR;
%alphaSICR = betaISICR/2 - gammaSICR/2 - mu - sigmaSICR/2 - thetaSICR/2 + (betaISICR^2 - 2*betaISICR*gammaSICR - 2*betaISICR*sigmaSICR + 2*betaISICR*thetaSICR +...
%    gammaSICR^2 + 2*gammaSICR*sigmaSICR - 2*gammaSICR*thetaSICR + sigmaSICR^2 - 2*sigmaSICR*thetaSICR + 4*betaCSICR*sigmaSICR + thetaSICR^2)^(1/2)/2;
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

%%
% mu = 0.0335;
% %SEIIS1 %chlamydia
% R1     = 2.2280;
% sigma1 = 35.7216;
% nu1    = 1.5610;
% gamma1 = 37.9999; 
% eps1   = 0.2426;
% beta1  = 2.2280;
% [~,~,alpha1] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,0);
% 
% %SEIIS2 %gono
% R2     = 1.0778;
% sigma2 = 82.1555;
% nu2    = 1.2955;
% eps2   = 0.6087;
% gamma2 = 27.8012; 
% beta2  = 3.4154;
% [~,~,alpha2] = Rp_SEIIS_v4(beta2,nu2,eps2,sigma2,gamma2,mu,b,0);
% 
% %SICR 
% RSICR     = 1.1632;
% sigmaSICR = 5.9899;
% thetaSICR = 0.1667;
% gammaSICR = 0;
% betaCSICR = 0.1796;
% betaISICR = 1.5892;
% [~,~,alphaSICR] = Rp_SICR_v4(betaISICR,betaCSICR,thetaSICR,sigmaSICR,gammaSICR,mu,b,0);
% 
% %SEIIIS/S(syphilis)
% RS      = 1.0833;
% sigmaS  = 12.1955;
% tauS    = 5.9349;
% thetaS  = 2.1477;
% gamma1S = 0;
% gamma3S = 0.0437;
% nuS     = 0;
% betaS = 0.0800;
% [~,~,alphaS] = Rp_SEIIIS_v4(betaS,sigmaS,tauS,nuS,gamma1S,thetaS,gamma3S,mu,b,0);

%%
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

%%
vecC = linspace(-0.15,0.4,10);
%biasFactor=10;
%cleft=-0.15; cright=0.4;
%vecC = linspace(cleft,cright,100);
tStart = tic;
tic
[tab,tabco,tabcn,tabTimes] = findRhohat_v5(N,paramTab,mu,b,vecC,biasFactor);
toc
tEnd = toc(tStart)
%load handel
%load chirp 
load splat               
sound(y,Fs)
%gong           

cs = findThresholds_v4(nSEIIS,nSICR,nSEIIIS,tab,vecAlphas,vecC);
%plot_zones_v3bis
%%
computeC0Cnn = 1;
optFindCnnAndC0.errMax = 1e-1;
optFindCnnAndC0.Tmax = 100;
optFindCnnAndC0.aPrioriCnn = -0.09;
optFindCnnAndC0.aPrioriC0 = 0.3;
optFindRhohat.sampleSize = 15;%100;
[bestStrategy,worstStrategy,costStratTable,cnn,c0,costStratTable_sorted,tab2,tabco2,tabcn2] = findBestStrategy(1,optFindCnnAndC0,optFindRhohat,N,paramTab,mu,b,biasFactor,0);
%%
costStratTable_sorted
bestStrategy
%path = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v2\costElimTable_b_10-set11_ech100.txt';
%fileID = fopen(path,'w');
%fprintf(fileID,'%d\n',costStratTable_sorted);
%fclose(fileID);
%writetable(costStratTable_sorted(:,1:5),path,'WriteVariableNames',false)

%%
figure()
vecC=tab2.c;
tab=tab2; 
%1 disease
plot(vecC,tab.one(1).rhohat,':','DisplayName','Ct')
hold on
plot(vecC,tab.one(2).rhohat,':','DisplayName','Ng')
plot(vecC,tab.one(3).rhohat,':','DisplayName','HIV')
plot(vecC,tab.one(4).rhohat,':','DisplayName','syph.')

%2 diseases
plot(vecC,tab.two(1).rhohat,'-.','DisplayName','CtxNg')
plot(vecC,tab.two(2).rhohat,'-.','DisplayName','CtxHIV')
plot(vecC,tab.two(3).rhohat,'-.','DisplayName','Ctxsyph.')
plot(vecC,tab.two(4).rhohat,'-.','DisplayName','NgxHIV')
plot(vecC,tab.two(5).rhohat,'-.','DisplayName','Ngxsyph.')
plot(vecC,tab.two(6).rhohat,'-.','DisplayName','HIVxsyph.')

%3 diseases
plot(vecC,tab.three(1).rhohat,'--','DisplayName','CtxNgxHIV')
plot(vecC,tab.three(2).rhohat,'--','DisplayName','CtxNgxsyph.')
plot(vecC,tab.three(2).rhohat,'--','DisplayName','CtxHIVxsyph.')
plot(vecC,tab.three(3).rhohat,'--','DisplayName','NgxHIVxsyph.')

%4 diseases
plot(vecC,tab.four.rhohat,'-','DisplayName','CtxNgxHIVxsyph.')

legend()

%%
%Ct,Ng,HIV
c=-0.12;
vecRho=0:0.001:0.12;
options.TolFun=1e-4; f=1;
U123 = U123_SEIIS2SICR_v4(paramTab{1},paramTab{2},paramTab{3},mu,b,vecRho,c,biasFactor,'fsolve',options);
U_3d = U_SEIIS2SICR_v4(paramTab{1},paramTab{2},paramTab{3},mu,b,vecRho,c,biasFactor);
U12 = U12_SEIIS2_v4(paramTab{1},paramTab{2},mu,b,vecRho,c,biasFactor,'fsolve',options);
U_2d = U_SEIIS2_v4(paramTab{1},paramTab{2},mu,b,vecRho,c,biasFactor);
U1 = U1_SEIISv4(paramTab{1},mu,b,vecRho,c,biasFactor);
U_1d = U_SEIIS_v4(paramTab{1},mu,b,vecRho,c,biasFactor);
figure(2)
plot(vecRho,U123,'DisplayName','U_{CtNgHIV}')
hold on
plot(vecRho,U123,'DisplayName','U_{123-CtNgHIV}')
plot(vecRho,U_3d,'DisplayName','U_{CtNgHIV}')
plot(vecRho,U12,'DisplayName','U_{12-CtNg}')
plot(vecRho,U_2d,'DisplayName','U_{CtNg}')
plot(vecRho,U1,'DisplayName','U_{1-Ct}')
plot(vecRho,U_1d,'DisplayName','U_{Ct}')
legend()

%% Recherche des couts d'elimination avec fminsearch
syms c
%notation : c_i_ijk : cout de l'elimination de la maladie i dans le modele ijk
alpha1 = paramCt.alpha;
alpha2 = paramNg.alpha;
alpha3 = paramHIV.alpha;
alpha4 = paramS.alpha;
% 1 maladie : directement donné par les equations

% 2 maladies
startPoint=0;
% Ct x Ng 
tStart = tic;
N=2; param2{1}=paramCt; param2{2}=paramNg; 
%minAlpha=min([param2{1}.alpha,param2{2}.alpha]);
fun1 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha1) ;
c_1_12 = fzero(fun1,startPoint)
fun2 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha2) ;
c_2_12 = fzero(fun2,startPoint)

% Ct x HIV
N=2; param2{1}=paramCt; param2{2}=paramHIV; 
fun1 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha1) ;
c_1_13 = fzero(fun1,startPoint)
fun3 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha3) ;
c_3_13 = fzero(fun3,startPoint)

% Ct x Syph
N=2; param2{1}=paramCt; param2{2}=paramS; 
fun1 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha1) ;
c_1_14 = fzero(fun1,startPoint)
fun4 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha4) ;
c_4_14 = fzero(fun4,startPoint)

% Ng x HIV
N=2; param2{1}=paramNg; param2{2}=paramHIV; 
fun2 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha2) ;
c_2_23 = fzero(fun2,startPoint)
fun3 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha3) ;
c_3_23 = fzero(fun3,startPoint)

% Ng x Syph
N=2; param2{1}=paramNg; param2{2}=paramS; 
fun2 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha2) ;
c_2_24 = fzero(fun2,startPoint)
fun4 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha4) ;
c_4_24 = fzero(fun4,startPoint)

% HIV x Syph
N=2; param2{1}=paramNg; param2{2}=paramS; 
fun3 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha3) ;
c_3_34 = fzero(fun3,startPoint)
fun4 = @(c) (findRhohat_v4_min(N,param2,mu,b,c,1) - alpha4) ;
c_4_34 = fzero(fun4,startPoint)

tEnd = toc(tStart)

findRhohat_v4_min(N,param2,mu,b,0,1) - minAlpha