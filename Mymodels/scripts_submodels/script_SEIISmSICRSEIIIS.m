%% This code generates the matrix and the ODE system associated to a SICRxSEIIISxSEIIS^m model
clear all;

%there is SICR, one SEIIIS and all the others are SEIIS
m = 1; %nbdiseases SEIIS (max m=4)
syms mu Lambda1 Lambda2 Lambda3 Lambda4 nu1 nu2 nu3 nu4 eps1 eps2 eps3 eps4 sigma1 sigma2 sigma3 sigma4 gamma1 gamma2 gamma3 gamma4 b rho;
syms lambdaHIV thetaHIV gammaHIV sigmaHIV
syms lambdaSyph sigmaSyph tauSyph thetaSyph gamma1Syph gamma3Syph nuSyph

%% below : to create ODE systems for the article 
%(the ODE systems in the matlab functions have been created without the 
% simplifications below)
gamma1Syph = 0;
nuSyph     = 0;
thetaHIV   = gammaHIV; %the variable theta has been changed by gamma
gammaHIV   = 0;        %the old variable gamma(0)=0
%%

%ne pas changer ici
n = 1; %nbdiseases SICR (VIH)
p = 1; %nbdiseases SEIIIS (syphilis)
nDis = m+n+p;
dis = 1:nDis;
nbBoxesSICR = 4;
nbBoxesSEIIS = 4;
nbBoxesSEIIIS = 5;

%parameters of the SEIIS
sigmas = [sigma1,sigma2,sigma3,sigma4];
gammas = [gamma1,gamma2,gamma3,gamma4];
Lambdas = [Lambda1,Lambda2,Lambda3,Lambda4];
nus = [nu1,nu2,nu3,nu4];
epss = [eps1,eps2,eps3,eps4];

nbCompartments = nbBoxesSICR^n*nbBoxesSEIIS^m*nbBoxesSEIIIS^p;
% creating all the compartments
% a compartment is defined bcy a n-tuple (x1,x2,x3,...xn) 
% xi takes value in 1 for S,2 for I,3 for J and concerns infection i
%// Sample data
z = 1:nbBoxesSEIIIS;
y = 1:nbBoxesSICR;   %'SICR';                
x = 1:nbBoxesSEIIS; %1:S,2:E,3:IA,4:IS

%disease 1 is a sicr, disease2 is a seiiis, the m others are seiis
if m==0
    compartments = combPerso(y,z);
elseif m==1
   compartments = combPerso(y,z,x);
elseif m==2
    compartments = combPerso(y,z,x,x);
elseif m==3
    compartments = combPerso(y,z,x,x,x);
elseif m==4
    compartments = combPerso(y,z,x,x,x,x);
end

%matrix containing all the flow rates initialization
M = sym(zeros(nbCompartments) - mu*eye(nbCompartments));

for j=dis %selon la maladie j, on recupere la matrice associee
    %if we look for the progression of disease 1 (SICR)
    if j==1 && m==0
        otherDiseaseStatesConstant = z';
    elseif j==1 && m==1
        otherDiseaseStatesConstant = combPerso(z,x);
    elseif j==1 && m==2
        otherDiseaseStatesConstant = combPerso(z,x,x);
    elseif j==1 && m==3
        otherDiseaseStatesConstant = combPerso(z,x,x,x);
    elseif j==1 && m==4
        otherDiseaseStatesConstant = combPerso(z,x,x,x,x);
    end
    
    %if we look for progression of disease 2 (SEIIIS)
    if j==2 && m==0
        otherDiseaseStatesConstant = y';
    elseif j==2 && m==1
        otherDiseaseStatesConstant = combPerso(y,x);
    elseif j==2 && m==2
        otherDiseaseStatesConstant = combPerso(y,x,x);
    elseif j==2 && m==3
        otherDiseaseStatesConstant = combPerso(y,x,x,x);
    elseif j==2 && m==4
        otherDiseaseStatesConstant = combPerso(y,x,x,x,x);
    end
    
    %if we look for progression of disease 3 and more (SIJS)
    if j>=3
        if  m==1
            otherDiseaseStatesConstant = combPerso(y,z);
        elseif m==2
            otherDiseaseStatesConstant = combPerso(y,z,x);
        elseif m==3
            otherDiseaseStatesConstant = combPerso(y,z,x,x);
        elseif m==4
            otherDiseaseStatesConstant = combPerso(y,z,x,x,x);
        end
    end
    
    if j==1
        mX = ODESICR(0,0,0,lambdaHIV,thetaHIV,gammaHIV,sigmaHIV,0,0,0);
    elseif j==2
        mX = ODESEIIIS(0,0,0,0,0,lambdaSyph,sigmaSyph,tauSyph,thetaSyph,gamma1Syph,gamma3Syph,nuSyph,0,0);
    else
        mX = ODESEIISv2(0,0,0,0,Lambdas(j-2),epss(j-2),nus(j-2),gammas(j-2),sigmas(j-2),0,0,0); %rho=0=mu=b
    end
    %il faut assigner cette matrice dans M, où 1 varie de 1 à 3 et où
    %toutes les autres maladies sont dans un état constant (ex : I3,I13,I3J1)
    for k=otherDiseaseStatesConstant'
        index = find(sum(compartments(:,dis(dis~=j))==k',2)==(nDis-1));
        %attention les etats 1 2 3 de j doivent être sorted
        M(index, index) = M(index, index) + mX;
    end
end

% adding the rho rate to all the compartments except those where ppl are only infected and 
% symptomatics for the STI : chlamydia,gono,etc. (e.g. Jc, Jg, RJc, RJg)
indexS  = find(sum(compartments,2)==nDis); % SSSS
indexR  = find(sum(compartments,2)==nDis+3 & compartments(:,1)==4); %RSSS
%below *version _2* : not infected by HIV (1) & not infected by syphilis (1) & (sympotmatic with
%the infection/s for which the indiviudal is infected (4))
% indexInfSTI = find((compartments(:,1)==1|compartments(:,1)==4)... %HIV : S or R
%     & compartments(:,2)==1 ... %syphilis : S
%     & sum(compartments(:,3:end),2)==4*m); %  STI(s) state : 4

%below *version _3* : once one individual is infected and symptomatic for
%one STI, s/he won't use VT but targeted testing instead based on his/her
%symptoms.

%1. we look for compartments where individuals are symptomatics (4 is for IS)
if m==0
    indexInfSTI=[];
elseif m==1
    indexInfSTI = find(compartments(:,3)==4);
elseif m==2
    indexInfSTI = find(compartments(:,3)==4 | compartments(:,4)==4);
end

%2. if individual are not (a) S (for all infections) or (b) R for HIV and S
% for the other infections or (c) symptomatics for the Ct and/or Ng
% then they can voluntary test
for k=1:nbCompartments
    if (k~=indexS && k~=indexR && ~ismember(k,indexInfSTI)) %k~=indexR a été ajouté après coup, mais ça doit s'annuler
        if (compartments(k,1)==1) %if not infected by HIV then combined testing ppl go to S
            M(k,k) = M(k,k) - rho;
            M(indexS,k) =  M(indexS,k) + rho;
        else % if infected by HIV then combined testing ppl go to RHIV
            M(k,k) = M(k,k) - rho;
            M(indexR,k) =  M(indexR,k) + rho;
        end
    end
end

% matrix to ODE system
[C,dC,eqn] = matToODE(nbCompartments,M);
dC.'
eqn.'
[C;compartments'].'

indexCHIV = find(compartments(:,1)==3);
indexIHIV = find(compartments(:,1)==2);
indexInfSyph = find(compartments(:,2)==3 | compartments(:,2)==4);
indexInfSEIIS = find(compartments(:,3)==3 | compartments(:,3)==4);

%% solving the ODE system
% RHIV = (betaIHIV*(thetaHIV+mu) + betaCHIV*sigmaHIV)/((mu+thetaHIV)*(sigmaHIV+gammaHIV+mu));
% R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rho))./((mu + sigma1+rho).*(gamma1 + mu + nu1).*(mu + nu1 + rho)); %ok
% RSyph = sigmaSyph*betaSyph*(tauSyph+thetaSyph+rho)/((thetaSyph+rho)*(gamma1Syph+rho+tauSyph+mu)*(sigmaSyph+mu+rho));
% b/mu
% 
% tspan = 0:.1:10;
% Y0 = (b/mu)/nbCompartments*ones(nbCompartments,1)';
% Y0 = [(b/mu)*ones(4,1);0.01*ones(80-4,1)]'
% options = odeset('AbsTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
% betaSEIIS = beta1;
% odefun = @(t,Y) ODE_SICRSEIIISSEIIS(t,Y,M,betaIHIV,betaCHIV,betaSyph,betaSEIIS,...
%     indexIHIV,indexCHIV,indexInfSyph,indexInfSEIIS,lambdaHIV,lambdaSyph,Lambda1,b);
% [ts,Ys] = ode45(odefun,tspan,Y0,options);
% T = Ys(end,:)
% sum(T)
% [1:nbCompartments;T] 


%% SICRxSEIIISxSEIJS (1) (versions 2 or 3)
clear all
compartments = combPerso(1:4,1:5,1:4);
indexInfSEIIS = find(compartments(:,3)==3 | compartments(:,3)==4);
nbCompartments = size(compartments,1);

[betaIHIV,betaCHIV,gammaHIV,sigmaHIV,thetaHIV,~,b,mu,rho] = random_parameters(true, true);
[betaSyph,sigmaSyph,gamma1Syph,gamma3Syph,tauSyph,thetaSyph,nuSyph,mu,rho] = random_parameters(true, true);
[beta1,beta2,gamma1,gamma2,nu1,nu2,sigma1,sigma2,rho] = random_parameters(true, true);
eps1=0.2;eps2=0.3;
RHIV = (betaIHIV*(thetaHIV+mu) + betaCHIV*sigmaHIV)/((mu+thetaHIV)*(sigmaHIV+gammaHIV+mu));
R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rho))./((mu + sigma1+rho).*(gamma1 + mu + nu1).*(mu + nu1 + rho)); %ok
RSyph = sigmaSyph*betaSyph*(tauSyph+thetaSyph+rho)/((thetaSyph+rho)*(gamma1Syph+rho+tauSyph+mu)*(sigmaSyph+mu+rho));
b/mu

%rho=0;

tspan = 0:.1:1000;
Y0 = (b/mu)*ones(80,1)';
options = odeset('AbsTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
% odefun = @(t,Y) ODE_SICRSEIIISSEIIS_3(t,Y,betaIHIV,betaCHIV,thetaHIV,gammaHIV,sigmaHIV,...
%                                       betaSyph,sigmaSyph,tauSyph,thetaSyph,gamma1Syph,gamma3Syph,nuSyph,...
%                                       beta1,nu1,eps1,sigma1,gamma1,mu,b,rho,indexInfSyph,indexInfSEIIS,indexIHIV,indexCHIV);
odefun = @(t,Y) ODE_SEIISSICRSEIIIS_3(t,Y,beta1,gamma1,nu1,sigma1,eps1,...
                                     betaIHIV,betaCHIV,gammaHIV,sigmaHIV,thetaHIV,...
                                     betaSyph,sigmaSyph,tauSyph,nuSyph,gamma1Syph,thetaSyph,gamma3Syph,mu,b,rho);
[ts,Ys] = ode45(odefun,tspan,Y0,options);
T = Ys(end,:)
sum(T)
[1:nbCompartments;T]

%% verifications
%HIV
indexSHIV = find(compartments(:,1)==1); indexIHIV = find(compartments(:,1)==2); indexCHIV = find(compartments(:,1)==3); indexRHIV = find(compartments(:,1)==4);
[sum(T(indexSHIV)),sum(T(indexIHIV)),sum(T(indexCHIV)),sum(T(indexRHIV))]

[ts,YsHIV] = ode45(@(t,Y) ODE_SICR(t,Y,b, betaIHIV, betaCHIV, sigmaHIV, thetaHIV, gammaHIV, mu,rho,'frequency'),tspan,[1,1,1,1], options);
THIV = YsHIV(end,:)

%STI 1
indexS1 = find(compartments(:,3)==1); indexE1 = find(compartments(:,3)==2); indexI1 = find(compartments(:,3)==3); indexJ1 = find(compartments(:,3)==4);
[sum(T(indexS1)),sum(T(indexE1)),sum(T(indexI1)),sum(T(indexJ1))]

R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rho))./((mu + sigma1+rho).*(gamma1 + mu + nu1).*(mu + nu1 + rho)); %ok
S1  = b/(mu*max(R1,1));
lambda1 = max(beta1*(R1-1)*(sigma1+rho+mu)/(beta1+(sigma1+rho+mu)*R1),0);
E1  = lambda1*S1/(sigma1+mu+rho);
I1  = (1-eps1)*sigma1*E1/(nu1+rho+mu);
J1  = eps1*sigma1*E1/(gamma1+nu1+mu);
[S1,E1,I1,J1]
[ts,Ys1] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta1,nu1,eps1,sigma1,gamma1,rho,mu),tspan,[1,1,1,1], options);
Ys1(end,:)

%syphilis
indexSs = find(compartments(:,2)==1); indexEs = find(compartments(:,2)==2);
indexI1s = find(compartments(:,2)==3); indexI2s = find(compartments(:,2)==4);
indexI3s = find(compartments(:,2)==5);
[sum(T(indexSs)),sum(T(indexEs)),sum(T(indexI1s)),sum(T(indexI2s)),sum(T(indexI3s))]
RSyph = sigmaSyph*betaSyph*(tauSyph+thetaSyph+rho+mu)/((thetaSyph+rho+mu)*(gamma1Syph+rho+tauSyph+mu)*(sigmaSyph+mu+rho));

numerateur = -(gamma3Syph+mu+nuSyph+rho)*((mu+rho+sigmaSyph)*(mu+rho+thetaSyph)*(gamma1Syph+mu+rho+tauSyph)-betaSyph*sigmaSyph*(rho+tauSyph+thetaSyph+mu));
denominateur = gamma1Syph*mu^2 + gamma3Syph*mu^2 + gamma1Syph*rho^2 + gamma3Syph*rho^2 + mu^2*nuSyph + 3*mu*rho^2 + 3*mu^2*rho + nuSyph*rho^2 + mu^2*sigmaSyph + mu^2*tauSyph + mu^2*thetaSyph + rho^2*sigmaSyph + rho^2*tauSyph +...
    rho^2*thetaSyph + mu^3 + rho^3 + gamma1Syph*gamma3Syph*mu + gamma1Syph*gamma3Syph*rho + gamma1Syph*gamma3Syph*thetaSyph +...
    gamma1Syph*mu*nuSyph + 2*gamma1Syph*mu*rho + 2*gamma3Syph*mu*rho + gamma1Syph*nuSyph*rho + gamma3Syph*mu*sigmaSyph +...
    gamma3Syph*mu*tauSyph + gamma1Syph*mu*thetaSyph + gamma3Syph*mu*thetaSyph + gamma1Syph*nuSyph*thetaSyph + gamma3Syph*rho*sigmaSyph +...
    gamma3Syph*rho*tauSyph + gamma1Syph*rho*thetaSyph + gamma3Syph*rho*thetaSyph + 2*mu*nuSyph*rho + gamma3Syph*sigmaSyph*tauSyph +...
    gamma3Syph*sigmaSyph*thetaSyph + mu*nuSyph*sigmaSyph + gamma3Syph*tauSyph*thetaSyph + mu*nuSyph*tauSyph + mu*nuSyph*thetaSyph +...
    2*mu*rho*sigmaSyph + 2*mu*rho*tauSyph + 2*mu*rho*thetaSyph + nuSyph*rho*sigmaSyph + nuSyph*rho*tauSyph + mu*sigmaSyph*tauSyph + ...
    nuSyph*rho*thetaSyph + mu*sigmaSyph*thetaSyph + mu*tauSyph*thetaSyph + nuSyph*sigmaSyph*tauSyph + nuSyph*sigmaSyph*thetaSyph +...
    nuSyph*tauSyph*thetaSyph + rho*sigmaSyph*tauSyph + rho*sigmaSyph*thetaSyph + rho*tauSyph*thetaSyph + sigmaSyph*tauSyph*thetaSyph;
lambdaS = max(numerateur/denominateur,0);

Ss  = b/(mu*max(RSyph,1)); Es = b*lambdaS/(mu*RSyph*(sigmaSyph+rho+mu));
I1s = sigmaSyph*Es/(tauSyph+gamma1Syph+rho+mu);
I2s = tauSyph*I1s/(thetaSyph+rho+mu);
I3s = thetaSyph*I2s/(nuSyph+gamma3Syph+rho+mu);
[Ss,Es,I1s,I2s,I3s]
[ts,YSyph] = ode45(@(t,Y) SEIIIS(t,Y,b,betaSyph,sigmaSyph, gamma1Syph, gamma3Syph, tauSyph, thetaSyph, nuSyph, mu,rho),tspan,[1,1,1,1,1], options);
YSyph(end,:)

%% SICRxSEIIISxSEIIS^2 (2) (versions_2)
clear all
compartments = combPerso(1:4,1:5,1:4,1:4);
indexCHIV = find(compartments(:,1)==3);
indexIHIV = find(compartments(:,1)==2);
indexInfSyph = find(compartments(:,2)==3 | compartments(:,2)==4);
indexInf1 = find(compartments(:,3)==3 | compartments(:,3)==4);
indexInf2 = find(compartments(:,4)==3 | compartments(:,4)==4);

nbCompartments = size(compartments,1);

[betaIHIV,betaCHIV,gammaHIV,sigmaHIV,thetaHIV,~,b,mu,rho] = random_parameters(true, true);
[betaSyph,sigmaSyph,gamma1Syph,gamma3Syph,tauSyph,thetaSyph,nuSyph,mu,rho] = random_parameters(true, true);
[beta1,beta2,gamma1,gamma2,nu1,nu2,sigma1,sigma2,rho] = random_parameters(true, true);
eps1=0.2;eps2=0.3;
RHIV = (betaIHIV*(thetaHIV+mu) + betaCHIV*sigmaHIV)/((mu+thetaHIV)*(sigmaHIV+gammaHIV+mu));
R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rho))./((mu + sigma1+rho).*(gamma1 + mu + nu1).*(mu + nu1 + rho)); %ok
R2 = (beta2*sigma2*(gamma2*(1-eps2) + mu + nu2 + eps2*rho))./((mu + sigma2+rho).*(gamma2 + mu + nu2).*(mu + nu2 + rho)); %ok
RSyph = sigmaSyph*betaSyph*(tauSyph+thetaSyph+rho)/((thetaSyph+rho)*(gamma1Syph+rho+tauSyph+mu)*(sigmaSyph+mu+rho));
b/mu

%rho=0;

tspan = 0:.1:1000;
Y0 = (b/mu)/80*ones(320,1)';
options = odeset('AbsTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
% odefun = @(t,Y) ODE_SICRSEIIISSEIIS2_3(t,Y,betaIHIV,betaCHIV,thetaHIV,gammaHIV,sigmaHIV,...
%                                       betaSyph,sigmaSyph,tauSyph,thetaSyph,gamma1Syph,gamma3Syph,nuSyph,...
%                                       beta1,nu1,eps1,sigma1,gamma1,...
%                                       beta2,nu2,eps2,sigma2,gamma2,mu,b,rho,...
%                                       indexInfSyph,indexInf1,indexInf2,indexIHIV,indexCHIV);

odefun = @(t,Y) ODE_SEIIS2SICRSEIIIS_3(t,Y,beta1,gamma1,nu1,sigma1,eps1,...
                                         beta2,gamma2,nu2,sigma2,eps2,...
                                         betaIHIV,betaCHIV,gammaHIV,sigmaHIV,thetaHIV,...
                                         betaSyph,sigmaSyph,tauSyph,gamma1Syph,thetaSyph,gamma3Syph,nuSyph,...
                                         mu,b,rho);
[ts,Ys] = ode45(odefun,tspan,Y0,options);
T = Ys(end,:)
sum(T) - b/mu
[1:nbCompartments;T] 

%% verifications

indexSHIV = find(compartments(:,1)==1); indexIHIV = find(compartments(:,1)==2); indexCHIV = find(compartments(:,1)==3); indexRHIV = find(compartments(:,1)==4);
[sum(T(indexSHIV)),sum(T(indexIHIV)),sum(T(indexCHIV)),sum(T(indexRHIV))]
[ts,YsHIV] = ode45(@(t,Y) ODE_SICR(t,Y,b, betaIHIV, betaCHIV, sigmaHIV, thetaHIV, gammaHIV, mu,rho,'frequency'),tspan,[1,1,1,1], options);
THIV = YsHIV(end,:)

%STI 1
indexS1 = find(compartments(:,3)==1); indexE1 = find(compartments(:,3)==2); indexI1 = find(compartments(:,3)==3); indexJ1 = find(compartments(:,3)==4);
[sum(T(indexS1)),sum(T(indexE1)),sum(T(indexI1)),sum(T(indexJ1))]
S1  = b/(mu*max(R1,1));
lambda1 = max(beta1*(R1-1)*(sigma1+rho+mu)/(beta1+(sigma1+rho+mu)*R1),0);
E1  = lambda1*S1/(sigma1+mu+rho);
I1  = (1-eps1)*sigma1*E1/(nu1+rho+mu);
J1  = eps1*sigma1*E1/(gamma1+nu1+mu);
[S1,E1,I1,J1]
[ts,Ys1] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta1,nu1,eps1,sigma1,gamma1,rho,mu),tspan,[1,1,1,1], options);
Ys1(end,:)

%STI 2
indexS2 = find(compartments(:,4)==1); indexE2 = find(compartments(:,4)==2); indexI2 = find(compartments(:,4)==3); indexJ2 = find(compartments(:,4)==4);
[sum(T(indexS2)),sum(T(indexE2)),sum(T(indexI2)),sum(T(indexJ2))]
S2  = b/(mu*max(R2,1));
lambda2 = max(beta2*(R2-1)*(sigma2+rho+mu)/(beta2+(sigma2+rho+mu)*R2),0);
E2  = lambda2*S2/(sigma2+mu+rho);
I2  = (1-eps2)*sigma2*E2/(nu2+rho+mu);
J2  = eps2*sigma2*E2/(gamma2+nu2+mu);
[S2,E2,I2,J2]
[ts,Ys2] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta2,nu2,eps2,sigma2,gamma2,rho,mu),tspan,[1,1,1,1], options);
Ys2(end,:)

%syphilis
indexSs = find(compartments(:,2)==1); indexEs = find(compartments(:,2)==2);
indexI1s = find(compartments(:,2)==3); indexI2s = find(compartments(:,2)==4);
indexI3s = find(compartments(:,2)==5);
[sum(T(indexSs)),sum(T(indexEs)),sum(T(indexI1s)),sum(T(indexI2s)),sum(T(indexI3s))]
RSyph = sigmaSyph*betaSyph*(tauSyph+thetaSyph+rho+mu)/((thetaSyph+rho+mu)*(gamma1Syph+rho+tauSyph+mu)*(sigmaSyph+mu+rho));

numerateur = -(gamma3Syph+mu+nuSyph+rho)*((mu+rho+sigmaSyph)*(mu+rho+thetaSyph)*(gamma1Syph+mu+rho+tauSyph)-betaSyph*sigmaSyph*(rho+tauSyph+thetaSyph+mu));
denominateur = gamma1Syph*mu^2 + gamma3Syph*mu^2 + gamma1Syph*rho^2 + gamma3Syph*rho^2 + mu^2*nuSyph + 3*mu*rho^2 + 3*mu^2*rho + nuSyph*rho^2 + mu^2*sigmaSyph + mu^2*tauSyph + mu^2*thetaSyph + rho^2*sigmaSyph + rho^2*tauSyph +...
    rho^2*thetaSyph + mu^3 + rho^3 + gamma1Syph*gamma3Syph*mu + gamma1Syph*gamma3Syph*rho + gamma1Syph*gamma3Syph*thetaSyph +...
    gamma1Syph*mu*nuSyph + 2*gamma1Syph*mu*rho + 2*gamma3Syph*mu*rho + gamma1Syph*nuSyph*rho + gamma3Syph*mu*sigmaSyph +...
    gamma3Syph*mu*tauSyph + gamma1Syph*mu*thetaSyph + gamma3Syph*mu*thetaSyph + gamma1Syph*nuSyph*thetaSyph + gamma3Syph*rho*sigmaSyph +...
    gamma3Syph*rho*tauSyph + gamma1Syph*rho*thetaSyph + gamma3Syph*rho*thetaSyph + 2*mu*nuSyph*rho + gamma3Syph*sigmaSyph*tauSyph +...
    gamma3Syph*sigmaSyph*thetaSyph + mu*nuSyph*sigmaSyph + gamma3Syph*tauSyph*thetaSyph + mu*nuSyph*tauSyph + mu*nuSyph*thetaSyph +...
    2*mu*rho*sigmaSyph + 2*mu*rho*tauSyph + 2*mu*rho*thetaSyph + nuSyph*rho*sigmaSyph + nuSyph*rho*tauSyph + mu*sigmaSyph*tauSyph + ...
    nuSyph*rho*thetaSyph + mu*sigmaSyph*thetaSyph + mu*tauSyph*thetaSyph + nuSyph*sigmaSyph*tauSyph + nuSyph*sigmaSyph*thetaSyph +...
    nuSyph*tauSyph*thetaSyph + rho*sigmaSyph*tauSyph + rho*sigmaSyph*thetaSyph + rho*tauSyph*thetaSyph + sigmaSyph*tauSyph*thetaSyph;
lambdaS = max(numerateur/denominateur,0);

Ss  = b/(mu*max(RSyph,1)); Es = b*lambdaS/(mu*RSyph*(sigmaSyph+rho+mu));
I1s = sigmaSyph*Es/(tauSyph+gamma1Syph+rho+mu);
I2s = tauSyph*I1s/(thetaSyph+rho+mu);
I3s = thetaSyph*I2s/(nuSyph+gamma3Syph+rho+mu);
[Ss,Es,I1s,I2s,I3s]
[ts,YSyph] = ode45(@(t,Y) SEIIIS(t,Y,b,betaSyph,sigmaSyph, gamma1Syph, gamma3Syph, tauSyph, thetaSyph, nuSyph, mu,rho),tspan,[1,1,1,1,1], options);
YSyph(end,:)
