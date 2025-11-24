mu = 0.0335;
%SEIIS1 %chlamydia
R1     = 1.0763;
sigma1 = 34.7068;
nu1    = 0.4067;
gamma1 = 33.3388; 
eps1   = 0.2728;
beta1  = 0.6479;
[~,~,alpha1] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,0);
%0.0333


%SICR 
RSICR     = 1.1948;
sigmaSICR = 6.0648;
thetaSICR = 0.2458;
gammaSICR = 0;
betaCSICR = 0.2334;
betaISICR = 2.2028;
[~,~,alphaSICR] = Rp_SICR_v4(betaISICR,betaCSICR,thetaSICR,sigmaSICR,gammaSICR,mu,b,0);
%0.0778

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

