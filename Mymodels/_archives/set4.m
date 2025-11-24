mu = 0.0335;
%SEIIS1 %chlamydia
R1     = 1.0860;
sigma1 = 40.2979;
nu1    = 0.5246;
gamma1 = 24.8922; 
eps1   = 0.0486;
beta1  = 0.6346;
[~,~,alpha1] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,0);

%SICR 
RSICR     = 1.1914;
sigmaSICR = 5.6398;
thetaSICR = 0.1986;
gammaSICR = 0; %0.39*365/(8.2*7); %0.4*(); %proportion de diag au stade precoce : 40
%ratioBeta = randPERT(8.4,9.1,9.6,1);
betaCSICR = 0.2000;
betaISICR = 1.8551;
[~,~,alphaSICR] = Rp_SICR_v4(betaISICR,betaCSICR,thetaSICR,sigmaSICR,gammaSICR,mu,b,0);

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

