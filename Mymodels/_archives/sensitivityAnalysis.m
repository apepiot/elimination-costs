%sensitivity analysis (kind of)

nSEIIS=2; nSICR=1; nSEIIIS=1;

mu = 1/35;b = 5;

%SEIIS1 %chlamydia
IC_PC     = [5/100,7/100,10/100]; %borne inf, mean, borne sup
IC_sigmaC = [362/21,365/11,365/7];
IC_nuC    = [365/595,365/497,365/431];
IC_gammaC = [365/23,365/(31.5-11),365/10]; 
IC_epsC   = [0,0.11,0/75];
IC_Ct{1} = IC_PC; IC_Ct{2} = IC_sigmaC; IC_Ct{3} = IC_nuC; IC_Ct{4} = IC_gammaC; IC_Ct{5} = IC_epsC; 

%SEIIS2 %gono
IC_PN     = [5/100,7/100,10/100];
IC_sigmaN = [365/14,365/5,365/1];
IC_nuN    = [12/12,12/6,12/4];
IC_epsN   = [0.85,0.89,0.93];
IC_gammaN = [365/10,365/5,365/2]; 
IC_N{1} = IC_PN; IC_N{2} = IC_sigmaN; IC_N{3} = IC_nuN; IC_N{4} = IC_gammaN; IC_N{5} = IC_epsN; 

%HIV 
IC_PH = [0.125,0.161,0.204];
IC_sigmaH = 365./[9.8,8.2,6.7].*7; %(1-0.39)*365/(8.2*7)%0.6*(); %8.2 semaines
IC_thetaH = 1./[4.28,4.19,4.09];%1/9.8;
IC_gammaH = [0,0,0];%0.39*365/(8.2*7); %0.4*(); %proportion de diag au stade precoce : 40
IC_ratioBeta = [8.4,9.1,9.6];
IC_H{1} = IC_PH; IC_H{2} = IC_sigmaH; IC_H{3} = IC_thetaH; IC_H{4} = IC_gammaH; IC_H{5} = IC_ratioBeta; 

%Syphilis
IC_PS = [0.049,0.066,0.089];
IC_sigmaS  = 365./[90,21,10];
IC_tauS    = 365./[60.2,45,29.8];%(1-0.55*0.31)*365/45; %0.6*
IC_thetaS  = 12./[12,3.6,1];
IC_gamma1S = [0,0,0];%(0.55*0.31)*365/45; %0.2*
IC_gamma3S = 1./[10,20,30];%1/5;
IC_S{1} = IC_PS; IC_S{2} = IC_sigmaS; IC_S{3} = IC_tauS; IC_S{4} = IC_thetaS;  IC_S{5} = IC_gamma1S; IC_S{6} = IC_gamma3S;

N = nSICR+nSEIIS+nSEIIIS;

IC = [IC_Ct,IC_N,IC_H,IC_S];

nbParameters = 21;

%chose some parameters
valuesMatrix = []; n=4;
for i=1:nbParameters
    valuesMatrix = [valuesMatrix;createSample(IC{i}(1),IC{i}(2),IC{i}(3),n,'uniformly')];
end

nbIter = 5;

for i=1:nbIter
    
    %change with a lhs for the future
    PC     = valuesMatrix(1,randsample(1:n,1));
    sigmaC = valuesMatrix(2,randsample(1:n,1));
    nuC     = valuesMatrix(3,randsample(1:n,1));
    gammaC  = valuesMatrix(4,randsample(1:n,1)); 
    epsC    = valuesMatrix(5,randsample(1:n,1));
    
    PN      = valuesMatrix(6,randsample(1:n,1));
    sigmaN  = valuesMatrix(7,randsample(1:n,1));
    nuN     = valuesMatrix(8,randsample(1:n,1));
    gammaN  = valuesMatrix(9,randsample(1:n,1));
    epsN    = valuesMatrix(10,randsample(1:n,1));
    
    PH      = valuesMatrix(11,randsample(1:n,1));
    sigmaH  = valuesMatrix(12,randsample(1:n,1));
    thetaH  = valuesMatrix(13,randsample(1:n,1));
    gammaH  = valuesMatrix(14,randsample(1:n,1));
    ratioBeta = valuesMatrix(15,randsample(1:n,1));
    
    PS      = valuesMatrix(16,randsample(1:n,1));
    sigmaS  = valuesMatrix(17,randsample(1:n,1));
    tauS    = valuesMatrix(18,randsample(1:n,1));
    thetaS  = valuesMatrix(19,randsample(1:n,1));
    gamma1S = valuesMatrix(20,randsample(1:n,1));
    gamma3S = valuesMatrix(21,randsample(1:n,1));
    
    %% mettre tout ça dans un autre code
    % calculated with the chosen parameters
    RC     = 1./(1-PC);
    betaC  = RC*(sigmaC+mu)*(gammaC+nuC+mu)*(nuC+mu)/(sigmaC*(gammaC*(1-epsC)+mu+nuC));
    souslaracine = @(gamma0,beta,nu,p,sigma) (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2);
    alphaC = max((betaC*epsC*sigmaC + sqrt(souslaracine(gammaC,betaC,nuC,epsC,sigmaC)))/(2*(gammaC+mu+nuC)) - (2*mu+nuC+sigmaC)/2,0);
    RN     = 1./(1-PN);
    betaN  = RN*(sigmaN+mu)*(gammaN+nuN+mu)*(nuN+mu)/(sigmaN*(gammaN*(1-epsN)+mu+nuN));
    alphaN = (betaN*epsN*sigmaN + sqrt(souslaracine(gammaN,betaN,nuN,epsN,sigmaN)))/(2*(gammaN+mu+nuN)) - (2*mu+nuN+sigmaN)/2;
    RH     = 1./(1-PH);%2.2; %P=1-1/R so R=1/(1-P);
    betaCH = RH*(sigmaH+gammaH+mu)*(thetaH+mu)/(ratioBeta*(thetaH+mu)+1);
    betaIH = ratioBeta*betaCH;
    alphaH = betaIH/2 - gammaH/2 - mu - sigmaH/2 - thetaH/2 + (betaIH^2 - 2*betaIH*gammaH - 2*betaIH*sigmaH + 2*betaIH*thetaH +...
        gammaH^2 + 2*gammaH*sigmaH - 2*gammaH*thetaH + sigmaH^2 - 2*sigmaH*thetaH + 4*betaCH*sigmaH + thetaH^2)^(1/2)/2;
    RS      = 1/(1-PS);
    nuS     = 0;
    betaS = RS*(thetaS+mu)*(gamma1S+tauS+mu)*(sigmaS+mu)/(sigmaS*(tauS+thetaS+mu));
    Rpfun  = @(rho) (sigmaS*betaS*(tauS+thetaS+rho+mu)./((thetaS+rho+mu).*(gamma1S+rho+tauS+mu).*(sigmaS+rho+mu))-1);
    alphaS = fzero(Rpfun, 0); 
    
    %
    paramSEIIS = [beta1,gamma1,nu1,sigma1,eps1,alpha1,R1;...
                  beta2,gamma2,nu2,sigma2,eps2,alpha2,R2];
    paramSICR = [betaIH,betaCH,gammaH,sigmaH,thetaH,alphaH,RH,0,0];
    paramS    = [betaS,sigmaS,tauS,gamma1S,thetaS,gamma3S,nuS,alphaS,RS,0,0];
    
    
    % This section gives rhohat=argmaxU for a given model
    xleft  = -1.5;
    xright = 1.5;
    vecC = xleft:(xright-xleft)/50:xright;
    tic; biasFactor=1;
    
    [tab,tabco,tabcn,tabTimes] = findRhohat_v2(nSEIIS,nSICR,nSEIIIS,paramSEIIS,paramSICR,paramS,mu,1,vecC,biasFactor);
end
tps = toc;


