function [paramTab,mu,vecAlphas] = sampleParameters_v3(Ct,Ng,HIV,Syph,b)
%difference with sampleParameters_v1: parameters (in particular beta) are
%calculated at the baseline voluntary testing rate

%create a table of parameters with a class structure
mu = 1/randPERT(27.2,30.6,33.7,1);      %check

%SEIIS1 %chlamydia
P1VT    = randPERT(4.0,5.3,6.6,1)/100;      %check  %voluntary testing is included, IC has been computed with a binomial distribution
R1VT    = 1/(1-P1VT);                       %check
sigma1  = 365./randPERT(7,11,21,1);         %check
nu1     = 12/randPERT(6,12,36,1);           %check
gamma1  = 365/randPERT(5,14,20,1);          %365/(31.5-1./sigma1); %check
eps1    = randPERT(0,0.11,0.5,1);           %check
rhob1   = 1/randPERT(6.77,8.14,9.51,1);     %check
beta1   = R1VT/((sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rhob1))./((mu + sigma1 + rhob1).*(gamma1 + mu + nu1).*(mu + nu1 + rhob1))); %check
[R10,~,alpha1] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,0); %check
P10     = 1-1/R10;

%SEIIS2 %gono
P2VT     = randPERT(3.3,4.7,6.8,1)/100; %check %voluntary testing is included
R2VT     = 1/(1-P2VT);                  %check
sigma2   = 365/randPERT(1,5,14,1);      %check
nu2      = 12/randPERT(4,6,12,1);       %check
gamma2   = 365/randPERT(5,14,20,1);     %check
eps2     = randPERT(0.5,0.8,1,1);       %check
rhob2    = 1/randPERT(6.77,8.14,9.51,1);%check
beta2    = R2VT/((sigma2*(gamma2*(1-eps2) + mu + nu2 + eps2*rhob2))./((mu + sigma2 + rhob2).*(gamma2 + mu + nu2).*(mu + nu2 + rhob2)));
[R20,~,alpha2] = Rp_SEIIS_v4(beta2,nu2,eps2,sigma2,gamma2,mu,b,0); %check
P20      = 1-1/R20;

%_sict (before the introduction of PrEP)
%prevalence and not undiagnosed prevalence, with voluntary testing
PHIV_prevagay_VT = randPERT(12.5,16.1,20.4,1)./100;                              %check
%Lambda_sict = randPERT(4,7,10,1)/100; %incidence
R_sict_VT   = 1/(1-PHIV_prevagay_VT);           %P=1-1/R so R=1/(1-P);          %check
sigma_sict  = 52/randPERT(6.7,8.2,9.8,1);       %(1-0.39)*365/(8.2*7)%0.6*();   %check 
rhobh       = 1/randPERT(1.6,1.9,2.2,1);                                        %check
theta_sict0 = 1/randPERT(4,4.4,10,1);           %1/9.8;                         %check, theta_sict0: basé sur les symptomes
theta_sict  = theta_sict0+rhobh;
gamma_sict0 = 0;                                %I->T                           %check
gamma_sict  = rhobh+gamma_sict0;
ratioBeta  = randPERT(8.4,9.1,9.6,1);           %betaI/betaC                    %check
betaC_sict  = R_sict_VT*(sigma_sict+gamma_sict+mu)*(theta_sict+mu)/(ratioBeta*(theta_sict+mu)+sigma_sict); %check
%betaC_sict = (Lambda_sict+mu)./mu.*(sigma_sict+rhobh+mu)*(theta_sict+rhobh+mu)./(ratioBeta*(theta_sict+rhobh+mu)+sigma_sict);
betaI_sict = ratioBeta*betaC_sict; %check
%alpha_sict = betaI_sict/2 - gamma_sict/2 - mu - sigma_sict/2 - theta_sict/2 + (betaI_sict^2 - 2*betaI_sict*gamma_sict - 2*betaI_sict*sigma_sict + 2*betaI_sict*theta_sict +...
%    gamma_sict^2 + 2*gamma_sict*sigma_sict - 2*gamma_sict*theta_sict + sigma_sict^2 - 2*sigma_sict*theta_sict + 4*betaC_sict*sigma_sict + theta_sict^2)^(1/2)/2;
[R0_sict,~,alpha_sict] = Rp_SICR_v4(betaI_sict,betaC_sict,theta_sict0,sigma_sict,gamma_sict0,mu,b,0); %check
P0_sict = 1-1/R0_sict;

%undiagnosed prevalence :
Pund_base = mu*(theta_sict+sigma_sict+mu)/(betaI_sict*(theta_sict+mu)+betaC_sict*sigma_sict)*(R_sict_VT-1);
Pund_0    = mu*(theta_sict0+sigma_sict+mu)/(betaI_sict*(theta_sict0+mu)+betaC_sict*sigma_sict)*(R0_sict-1);

%SEIIIS/S(syphilis)
PS_VT     = randPERT(4.9,6.6,8.9,1)/100;                            %check
RS_VT     = 1/(1-PS_VT);                                            %check
rhobs     = 1/randPERT(2.06,2.48,2.89,1);                           %check
sigmaS    = 365/randPERT(10,25,90,1);                               %check
tauS      = 365/randPERT(10,46,130,1);   %(1-0.55*0.31)*365/45;     %check
thetaS    = 12/randPERT(1,3.6,12,1);                                %check
gamma1S0  = 0;                           %(0.55*0.31)*365/45;       %check
gamma1S   = gamma1S0+rhobs;
gamma3S0  = 1/randPERT(10,20,30,1);%1/5;                            %check 
gamma3S   = gamma3S0+rhobs;
nuS       = 0;                                                      %check
betaS = RS_VT/((sigmaS.*((gamma3S+mu+nuS).*(mu+rhobs+tauS+thetaS)+tauS*thetaS))./((mu+rhobs+sigmaS).*(mu+rhobs+thetaS).*(gamma3S+mu+nuS).*(gamma1S+mu+tauS)));
[R0S,~,alphaS] = Rp_SEIIIS_v4(betaS,sigmaS,tauS,nuS,gamma1S0,thetaS,gamma3S0,mu,b,0);
P0S = 1-1/R0S;

%%
%Ct
paramCt.beta=beta1; paramCt.sigma=sigma1; paramCt.gamma=gamma1;
paramCt.nu=nu1;paramCt.eps=eps1;
paramCt.alpha=alpha1;
paramCt.R_base=R1VT;paramCt.R0=R10;
paramCt.rhob = rhob1;
paramCt.disease = 'Ct'; paramCt.modelType = 'SEIIS'; paramCt.mini_d = 'c';
paramCt.P_base = P1VT; paramCt.P0 = P10;
paramCt.mu = mu; paramCt.pi=b;

%Ng
paramNg.beta = beta2;paramNg.sigma=sigma2;paramNg.gamma=gamma2;paramNg.nu=nu2;
paramNg.eps=eps2;paramNg.alpha=alpha2;
paramNg.R_base=R2VT; paramNg.R0=R20;
paramNg.rhob = rhob2;
paramNg.disease = 'Ng'; paramNg.modelType = 'SEIIS'; paramNg.mini_d = 'g';
paramNg.P_base = P2VT; paramNg.P0 = P20;
paramNg.mu = mu; paramNg.pi=b;

%HIV
paramHIV.betaI = betaI_sict; paramHIV.betaC = betaC_sict; 
paramHIV.gamma0 = gamma_sict0; paramHIV.theta0 = theta_sict0;
paramHIV.gamma_base = gamma_sict; paramHIV.theta_base = theta_sict;
paramHIV.sigma=sigma_sict; paramHIV.alpha_prev = alpha_sict;
paramHIV.rhob = rhobh;
paramHIV.R_prev_base = R_sict_VT; paramHIV.R_prev_0 = R0_sict; 
paramHIV.Pund_prev_base = Pund_base; paramHIV.Pund_prev_0 = Pund_0;
paramHIV.Ptot_prev_base = PHIV_prevagay_VT; paramHIV.Ptot_prev_0 = P0_sict;
paramHIV.disease = 'HIV'; paramHIV.modelType = 'SICT'; paramHIV.mini_d = 'h';
paramHIV.mu = mu; paramHIV.pi=b;

%Syphilis
paramS.beta = betaS; paramS.sigma=sigmaS; paramS.tau=tauS; 
paramS.gamma1_base=gamma1S; paramS.gamma10=gamma1S0;
paramS.gamma3_base=gamma3S; paramS.gamma30=gamma3S0;
paramS.theta=thetaS;
paramS.nu=nuS; paramS.alpha=alphaS;
paramS.rhob = rhobs;
paramS.R_base = RS_VT; paramS.R0 = R0S;
paramS.disease = 'syphilis'; paramS.modelType = 'SEIIIS'; paramS.mini_d = 's';
paramS.P_base = PS_VT; paramS.P0 = P0S; 
paramS.mu = mu; paramS.pi=b;

vecAlphas = [paramCt.alpha(Ct),paramNg.alpha(Ng),paramHIV.alpha_prev(HIV),paramS.alpha(Syph)];

paramTabAll{1}=paramCt;
paramTabAll{2}=paramNg;
paramTabAll{3}=paramHIV;
paramTabAll{4}=paramS;

paramTab = paramTabAll([logical(Ct),logical(Ng),logical(HIV),logical(Syph)]);

end
