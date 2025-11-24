% Code which solve the utlity maximization pb in a n-diseases model (n<=4
% for now)
%DIFFERENCE WITH RHOHAT_v3 : we implement RHOHAT_v3 in a function to
%perform a sensitivity analysis on the parameters.

% one SICR : HIV
% one SEIIIS : syphilis
% one or two SEIIS : chlamydia, gonorrhea

%2min30 for each sim 
%100s

%modification compared to RHOHAT_v3 : structure based for parameters
%% Parameters
clear all;%close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ct=true; Ng=true; HIV=true; Syph=true;
biasFactor = 1;

nSEIIS=Ct+Ng; nSICR=HIV; nSEIIIS=Syph;
N = nSICR+nSEIIS+nSEIIIS;
%STIChoice = 'Ct'; STIChoiceMini = 'Ct';%if one STI considered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 1/30.6;b = 5;
%%
nbSim = 10;

%%
bestStrat={};worstStrat={};tps1=tic;tabRecapStratCurrentSim=[];
path  = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\numSim.txt';
path2 = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\';
fileID = fopen(path,'r');
formatSpec = '%d';
A = fscanf(fileID,formatSpec);
lastSim = A(end); %a recuperer dans un fichier
fclose(fileID);

for i=(lastSim+1):(nbSim+lastSim+1)
    i
%SEIIS1 %chlamydia
P1     = randPERT(6.4,7.4,8.4,1)/100;
R1     = 1/(1-P1);
sigma1 = 365./randPERT(7,11,21,1);
nu1    = 12/randPERT(6,12,36,1);
gamma1 = 365/randPERT(5,14,20,1);%365/(31.5-1./sigma1); 
eps1   = randPERT(0,0.11,0.5,1);
beta1  = R1*(sigma1+mu)*(gamma1+nu1+mu)*(nu1+mu)/(sigma1*(gamma1*(1-eps1)+mu+nu1));
souslaracine = @(gamma0,beta,nu,p,sigma) (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2);
alpha1 = max((beta1*eps1*sigma1 + sqrt(souslaracine(gamma1,beta1,nu1,eps1,sigma1)))/(2*(gamma1+mu+nu1)) - (2*mu+nu1+sigma1)/2,0);

%SEIIS2 %gono
P2     = randPERT(6,7,8,1)/100;
R2     = 1/(1-P2);
sigma2 = 365/randPERT(1,5,14,1);
nu2    = 12/randPERT(4,6,12,1);
eps2   = randPERT(0.5,0.8,1,1);
gamma2 = 365/randPERT(5,14,20,1);%365/3; 
beta2  = R2*(sigma2+mu)*(gamma2+nu2+mu)*(nu2+mu)/(sigma2*(gamma2*(1-eps2)+mu+nu2));
alpha2 = (beta2*eps2*sigma2 + sqrt(souslaracine(gamma2,beta2,nu2,eps2,sigma2)))/(2*(gamma2+mu+nu2)) - (2*mu+nu2+sigma2)/2;

%SICR 
PHIV = randPERT(12.0,14.3,16.9,1)/100;
RSICR     = 1/(1-PHIV);%2.2; %P=1-1/R so R=1/(1-P);
sigmaSICR = 365/randPERT(6.7,8.2,9.8,1); %(1-0.39)*365/(8.2*7)%0.6*(); %8.2 semaines
thetaSICR = 1/randPERT(4,4.4,10,1);%1/9.8;
gammaSICR = 0;%0.39*365/(8.2*7); %0.4*(); %proportion de diag au stade precoce : 40
ratioBeta = 9.1;
betaCSICR = RSICR*(sigmaSICR+gammaSICR+mu)*(thetaSICR+mu)/(ratioBeta*(thetaSICR+mu)+1);
betaISICR = ratioBeta*betaCSICR;
alphaSICR = betaISICR/2 - gammaSICR/2 - mu - sigmaSICR/2 - thetaSICR/2 + (betaISICR^2 - 2*betaISICR*gammaSICR - 2*betaISICR*sigmaSICR + 2*betaISICR*thetaSICR +...
    gammaSICR^2 + 2*gammaSICR*sigmaSICR - 2*gammaSICR*thetaSICR + sigmaSICR^2 - 2*sigmaSICR*thetaSICR + 4*betaCSICR*sigmaSICR + thetaSICR^2)^(1/2)/2;

%SEIIIS/S(syphilis)
PS      = randPERT(6.7,7.7,8.7,1)/100;
RS      = 1/(1-PS);
sigmaS  = 365/randPERT(10,25,95,1);
tauS    = 365/randPERT(10,46,130,1);%(1-0.55*0.31)*365/45; %0.6*
thetaS  = 12/randPERT(1,3.6,12,1);
gamma1S = (0.55*0.31)*365/45; %0.2*
gamma3S = 1/randPERT(10,20,30,1);%1/5;
nuS     = 0;
betaS = RS*(thetaS+mu)*(gamma1S+tauS+mu)*(sigmaS+mu)/(sigmaS*(tauS+thetaS+mu));
Rpfun  = @(rho) (sigmaS*betaS*(tauS+thetaS+rho+mu)./((thetaS+rho+mu).*(gamma1S+rho+tauS+mu).*(sigmaS+rho+mu))-1);
alphaS = fzero(Rpfun, 0); 


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
computeC0Cnn = 1;
optFindCnnAndC0.errMax = 1e-1;
optFindCnnAndC0.Tmax = 100;
optFindCnnAndC0.aPrioriCnn = -0.1;
optFindCnnAndC0.aPrioriC0 = 0.4;
optFindRhohat.sampleSize = 10;
N = 4;
[bestStrat(i),worstStrat(i),tabRecapStrat]= findBestStrategy(computeC0Cnn,optFindCnnAndC0,optFindRhohat,N,paramTab,mu,b,biasFactor,0);

tabRecapStrat.nbSim = i*ones(15,1);
tabRecapStratCurrentSim = [tabRecapStratCurrentSim;tabRecapStrat];

%on sauvegarde les resultats tous les x simulations
if (mod(i,3)==0)
    fileID = fopen(path,'w');
    fprintf(fileID,'%d\n',i);
    fclose(fileID);
    fileName = ['Sim_',num2str(lastSim+1),'_to_',num2str(i)];
    %fileID2 = fopen(path2,'w');
    %fprintf(fileID2,'%s',tabRecapStratCurrentSim)
    writetable(tabRecapStratCurrentSim,[path2,fileName],'WriteVariableNames',false)
    tabRecapStratCurrentSim=[]; % reinitializing to empty
end

end



tps2 = toc(tps1)
bestStrat
worstStrat
% 
% T=bestStrat;
% [g, T2] = findgroups(T);
% h = splitapply(@numel,T,g);  
% B=table((T2'),(h'))
%1:'Chlam', 2:'Gono', 3:'HIV', 4:'Syphilis'
