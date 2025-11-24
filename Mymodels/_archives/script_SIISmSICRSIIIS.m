%% This code finds the matrix and the ODE system associated to a SIJS^mxSICR^n model
% attention : in this version there is rho from I to S
clear all;

%example when there is only one disease SICR and all the others are SIIS
m = 1; %nbdiseases SIIS (max m=4)

%ne pas changer ici
n = 1; %nbdiseases SICR (VIH)
p = 1; %nbdiseases SEIIIS (syphilis)
nDis = m+n+p;
dis = 1:nDis;
nbBoxesSICR = 4;
nbBoxesSIIS = 3;
nbBoxesSEIIIS = 5;

syms mu Lambda1 Lambda2 Lambda3 Lambda4 nu1 nu2 nu3 nu4 eps1 eps2 eps3 eps4 sigma1 sigma2 sigma3 sigma4 gamma1 gamma2 gamma3 gamma4 b rho;
syms lambdasicr thetasicr gammasicr sigmasicr
syms lambdasiiis sigmasiiis tausiiis thetasiiis gamma1siiis gamma3siiis nusiiis

%parameters of the SEIIIS
sigmas = [sigma1,sigma2,sigma3,sigma4];
gammas = [gamma1,gamma2,gamma3,gamma4];
Lambdas = [Lambda1,Lambda2,Lambda3,Lambda4];
nus = [nu1,nu2,nu3,nu4];
epss = [eps1,eps2,eps3,eps4];

nbCompartments = nbBoxesSICR^n*nbBoxesSIIS^m*nbBoxesSEIIIS^p;
% creating all the compartments
% a compartment is defined bcy a n-tuple (x1,x2,x3,...xn) 
% xi takes value in 1 for S,2 for I,3 for J and concerns infection i
%// Sample data
z = 1:nbBoxesSEIIIS;
y = 1:nbBoxesSICR;%'SICR';                 %// Set of possible letters                     %// Length of each permutation
x = 1:nbBoxesSIIS;

%disease 1 is a sicr, disease2 is a siiis, the m others are sijs
if m==1
   compartments = combPerso(y,z,x);
elseif m==2
    compartments = combPerso(y,z,x,x);
elseif m==3
    compartments = combPerso(y,z,x,x,x);
elseif m==4
    compartments = combPerso(y,z,x,x,x,x);
end

%matrix containing all the flow rates initialization
M = zeros(nbCompartments) - mu*eye(nbCompartments);

for j=dis %selon la maladie j, on recupere la matrice associee
    %if we look for the progression of disease 1 (SICR)
    if j==1 && m==1
        otherDiseaseStatesConstant = combPerso(z,x);
    elseif j==1 && m==2
        otherDiseaseStatesConstant = combPerso(z,x,x);
    elseif j==1 && m==3
        otherDiseaseStatesConstant = combPerso(z,x,x,x);
    elseif j==1 && m==4
        otherDiseaseStatesConstant = combPerso(z,x,x,x,x);
    end
    
    %if we look for progression of disease 2 (SEIIIS)
    if j==2 && m==1
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
        mX = ODESICR(0,0,0,lambdasicr, thetasicr, gammasicr,sigmasicr,b,0,0);
    elseif j==2
        mX = ODESEIIIS(0,0,0,0,0,lambdasiiis,sigmasiiis, tausiiis, thetasiiis, gamma1siiis,gamma3siiis,nusiiis,0,0);
    else
        mX = ODESIIS(0,0,0,Lambdas(j-2),epss(j-2),gammas(j-2),nus(j-2),sigmas(j-2),b,0,0);
    end
    %il faut assigner cette matrice dans M, où 1 varie de 1 à 3 et oùxs
    %toutes les autres maladies sont dans un état constant (ex : I3,I13,I3J1)
    for k=otherDiseaseStatesConstant'
        index = find(sum(compartments(:,dis(dis~=j))==k',2)==(nDis-1));
        %attention les etats 1 2 3 de j doivent être sorted
        %a faire
        M(index, index) = M(index, index) + mX;
    end
end

% adding rho to all the compartments (even Ii, ? to be discussed)
% faire une fonction pour ça ? 
indexS  = sum(compartments,2)==nDis;
for k=1:nbCompartments
    if (k==indexS) %if notS
        if (compartments(k,1)==1) %if not infected by HIV then combined testing ppl go to S
            M(k,k) = M(k,k) - rho;
            M(indexS,k) =  M(indexS,k) + rho;
        else % if infected by HIV then combined testing ppl go to RHIV
            M(k,k) = M(k,k) - rho;
            M(indexS,k) =  M(indexS,k) + rho;
        end
    end
end

% matrix to ODE system
[C,dC] = matToODE(nbCompartments,M);
dC.'
[C;compartments'].'

