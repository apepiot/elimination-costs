% Code which solve the utlity maximization pb in a n-diseases model (n<=3
% for now)

%% Parameters
clear all;close all;
nSIS=2; nSIR=1; nSICAT=0;

mu = 1/35;b = 5;

%SIS1
R1SIS     = 1.1;
gamma1SIS = 365/31; 
beta1SIS  = R1SIS*(gamma1SIS+mu);
s1SIS     = 1;
alpha1SIS = beta1SIS/s1SIS*(1-1/R1SIS);

%SIS2
R2SIS     = 1.1;
gamma2SIS = 6; %10 days
beta2SIS  = R2SIS*(gamma2SIS+mu);
s2SIS     = 1;
alpha2SIS = beta2SIS/s2SIS*(1-1/R2SIS);

%SIS3
R3SIS     = 7;
gamma3SIS = 6; %10 days
beta3SIS  = R3SIS*(gamma3SIS+mu);
s3SIS     = 1;
alpha3SIS = beta3SIS/s3SIS*(1-1/R3SIS);

%SIR1
R1SIR     = 3.05;
gamma1SIR = 12; 
beta1SIR  = R1SIR*(gamma1SIR+mu);
s1SIR     = 1;
alpha1SIR = beta1SIR/s1SIR*(1-1/R1SIR);


paramSIS = [beta1SIS,gamma1SIS,s1SIS;beta2SIS,gamma2SIS,s2SIS;beta3SIS,gamma3SIS,s3SIS];
paramSIR = [beta1SIR,gamma1SIR,s1SIR];
%paramSICAT = [beta1SICAT,gamma1SICAT,s1SICAT];
vecAlphaSIS=[alpha1SIS,alpha2SIS,alpha3SIS];
vecAlphaSIR=[alpha1SIR];
vecAlpha = [vecAlphaSIS(1:nSIS),vecAlphaSIR(1:nSIR)];
maxalpha = max(vecAlpha);
N = nSIR+nSIS+nSICAT;


%% This section gives rhohat=argmaxU for a given model
step = 0.001;
vecC = -1:step:1;
%example : SISxSIS 
cond=1; i=0;
while(cond)%interval vecC should be large enough
    vecC = vecC(1)*1.5:step:vecC(end);
    [tab,tabco,tabcn] = findRhohat(nSIS,nSIR,nSICAT,paramSIS(1:nSIS,:),paramSIR(1:nSIR,:),[],mu,b,vecC);
    cond = tab.rhohat(1)~=max(vecAlpha);
    i=i+1
end

%% c thresholds 
%pres = 0.1*abs(vecC(1)-vecC(end));

%% cswitch for the 3 disease model
paramBETA=[paramSIS(1:nSIS,1);paramSIR(1:nSIR,1)]; %add sicat
paramGAMMA=[paramSIS(1:nSIS,2);paramSIR(1:nSIR,2)]; %add sicat
paramS=[paramSIS(1:nSIS,3);paramSIR(1:nSIR,3)]; %add sicat

cs = findThresholds(nSIS,nSIR,nSICAT, tab, vecAlpha, paramBETA, paramGAMMA, paramS, 5, mu,vecC);
if (N>=2)
    ics1 = find(cs.cs1dis==vecC);
    ics2 = find(cs.cs2dis==vecC);
    if (N>=3)
        ics3 = find(cs.cs3dis==vecC);
    end
end


%% plot
plot_zones;



