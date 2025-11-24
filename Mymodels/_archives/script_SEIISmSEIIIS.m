%% This code generates the matrix and the ODE system associated to a SICRxSEIIISxSEIIS^m model
clear all;

%there is only one disease SICR and all the others are SIIS
m = 2; %nbdiseases SEIIS (max m=4)
syms mu Lambda1 Lambda2 Lambda3 Lambda4 nu1 nu2 nu3 nu4 eps1 eps2 eps3 eps4 sigma1 sigma2 sigma3 sigma4 gamma1 gamma2 gamma3 gamma4 b rho;
syms lambdaSyph sigmaSyph tauSyph thetaSyph gamma1Syph gamma3Syph nuSyph
% [betaIHIV,betaCHIV,gammaHIV,sigmaHIV,thetaHIV,~,b,mu,rho] = random_parameters(true, true);
% [betaSyph,sigmaSyph,gamma1Syph,gamma3Syph,tauSyph,thetaSyph,nuSyph,mu,rho] = random_parameters(true, true);
% [beta1,beta2,gamma1,gamma2,nu1,nu2,sigma1,sigma2,rho] = random_parameters(true, true);
% eps1=0.5;eps2=0.3;

%ne pas changer ici
p = 1; %nbdiseases SEIIIS (syphilis)
nDis = m+p;
dis = 1:nDis;
nbBoxesSEIIS = 4;
nbBoxesSEIIIS = 5;

%parameters of the SEIIS
sigmas = [sigma1,sigma2,sigma3,sigma4];
gammas = [gamma1,gamma2,gamma3,gamma4];
Lambdas = [Lambda1,Lambda2,Lambda3,Lambda4];
nus = [nu1,nu2,nu3,nu4];
epss = [eps1,eps2,eps3,eps4];

nbCompartments = nbBoxesSEIIS^m*nbBoxesSEIIIS^p;
% creating all the compartments
% a compartment is defined bcy a n-tuple (x1,x2,x3,...xn) 
% xi takes value in 1 for S,2 for I,3 for J and concerns infection i
z = 1:nbBoxesSEIIIS;
x = 1:nbBoxesSEIIS; %1:S,2:E,3:IA,4:IS

%disease 1 is a sicr, disease2 is a siiis, the m others are seiis

if m==1
    compartments = combPerso(z,x);
elseif m==2
    compartments = combPerso(z,x,x);
elseif m==3
    compartments = combPerso(z,x,x,x);
elseif m==4
    compartments = combPerso(z,x,x,x,x);
end

%matrix containing all the flow rates initialization
M = sym(zeros(nbCompartments) - mu*eye(nbCompartments));

for j=dis %selon la maladie j, on recupere la matrice associee   
    %if we look for progression of disease 3 and more (SIJS)  
    if j==1
        if m==1
            otherDiseaseStatesConstant = x';
        elseif m==2
            otherDiseaseStatesConstant = combPerso(x,x);
        elseif m==3
            otherDiseaseStatesConstant = combPerso(x,x,x);
        elseif m==4
            otherDiseaseStatesConstant = combPerso(x,x,x,x);
        end
    elseif j>=2
        if  m==1
            otherDiseaseStatesConstant = z';
        elseif m==2
            otherDiseaseStatesConstant = combPerso(z,x);
        elseif m==3
            otherDiseaseStatesConstant = combPerso(z,x,x);
        elseif m==4
            otherDiseaseStatesConstant = combPerso(z,x,x,x);
        end
    end
    if j==1
        mX = ODESEIIIS(0,0,0,0,0,lambdaSyph,sigmaSyph, tauSyph, thetaSyph, gamma1Syph,gamma3Syph,nuSyph,0,0);
    else
        mX = ODESEIISv2(0,0,0,0,Lambdas(j-1),epss(j-1),nus(j-1),gammas(j-1),sigmas(j-1),0,0,0);
    end
    %il faut assigner cette matrice dans M, où 1 varie de 1 à 3 et où
    %toutes les autres maladies sont dans un état constant (ex : I3,I13,I3J1)
    for k=otherDiseaseStatesConstant'
        index = find(sum(compartments(:,dis(dis~=j))==k',2)==(nDis-1));
        %attention les etats 1 2 3 de j doivent être sorted
        %a faire
        M(index, index) = M(index, index) + mX;
    end
end

% (version 2)adding the rho rate to all the compartments except those where ppl are only infected and 
% symptomatics for the STI :chlamydia, gono,etc. (e.g. Jc, Jg, RJc, RJg)
indexS  = find(sum(compartments,2)==nDis);
%below (version 2) : not infected by HIV (1) & not infected by syphilis (1) & (sympotmatic with
%the infection/s for which the indiviudal is infected (4))
indexInfSTI = find(sum(compartments(:,2:end),2)==4*m); %
%version 3
if m==0
    indexInfSTI=[];
elseif m==1
    indexInfSTI = find(compartments(:,2)==4);
elseif m==2
    indexInfSTI = find(compartments(:,2)==4 | compartments(:,3)==4);
end

for k=1:nbCompartments
    if (k~=indexS && ~ismember(k,indexInfSTI)) %if notS and if (infected by HIV or by syphilis or asymptomatic for another sti)
        M(k,k) = M(k,k) - rho;
        M(indexS,k) =  M(indexS,k) + rho;
    end
end

% matrix to ODE system
[C,dC,eqn] = matToODE(nbCompartments,M);
dC.'
eqn.'
[C;compartments'].'
