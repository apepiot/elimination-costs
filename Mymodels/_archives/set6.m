mu = 0.0335;
%SEIIS1 %chlamydia
R1     = 1.0864;
sigma1 = 32.3603;
nu1    = 1.2005;
gamma1 = 21.1798; 
eps1   = 0.3564;
beta1  = 2.0236;
[~,~,alpha1] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,0);
%0.1057

%SICR 
RSICR     = 1.1587;
sigmaSICR = 6.2020;
thetaSICR = 0.1763;
gammaSICR = 0;
betaCSICR = 0.1879;
betaISICR = 1.6689;
[~,~,alphaSICR] = Rp_SICR_v4(betaISICR,betaCSICR,thetaSICR,sigmaSICR,gammaSICR,mu,b,0);
%0.0431

%SEIIS2 %gono
P2     = NaN;
R2     = NaN;
sigma2 = NaN;
nu2    = NaN;
eps2   = NaN;
gamma2 = NaN;
beta2  = NaN;
alpha2 = NaN;

%SEIIIS/S(syphilis)
PS      = NaN;
RS      = NaN;
sigmaS  = NaN;
tauS    = NaN;
thetaS  = NaN;
gamma1S = NaN;
gamma3S = NaN;
nuS     = NaN;
betaS   = NaN;
alphaS  = NaN;

