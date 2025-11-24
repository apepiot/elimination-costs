function [paramTab,mu,vecAlphas] = sampleParameters(Ct,Ng,HIV,Syph,b)
%create a table of parameters with a class structure
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
gamma2 = 365/randPERT(5,14,20,1);%365/3;
eps2   = randPERT(0.5,0.8,1,1);
beta2  = R2/((sigma2*(gamma2*(1-eps2) + mu + nu2))./((mu+sigma2).*(gamma2 + mu + nu2).*(mu + nu2)));
[~,~,alpha2] = Rp_SEIIS_v4(beta2,nu2,eps2,sigma2,gamma2,mu,b,0);

%SICR
PHIV = randPERT(12.0,14.3,16.9,1)/100; %prevalence and not undiagnosed prevalence
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

%undiagnosed prevalence :
mu*(thetaSICR+sigmaSICR+mu)/(betaISICR*(thetaSICR+mu)+betaCSICR*sigmaSICR)*(RSICR-1);

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

paramTab = paramTabAll([logical(Ct),logical(Ng),logical(HIV),logical(Syph)]);

end

