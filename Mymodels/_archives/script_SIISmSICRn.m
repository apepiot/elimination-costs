%% This code finds the matrix and the ODE system associated to a SIJS^mxSICR model
clear all;

%example when there is only one disease SICR and all the others are SIJS
m = 1; %nbdiseases SIJS (max 3)
n = 1; %nbdiseases SICR (do not change that !!)
nDis = n+m;
dis = 1:nDis;
nbBoxesSICR = 4;
nbBoxesSIJS = 3;

syms mu Lambda1 Lambda2 Lambda3 Lambda4 nu1 nu2 nu3 nu4 eps1 eps2 eps3 eps4 sigma1 sigma2 sigma3 sigma4 gamma1 gamma2 gamma3 gamma4 b rho;
syms lambdasicr thetasicr gammasicr sigmasicr

sigmas = [sigma1,sigma2,sigma3,sigma4];
gammas = [gamma1,gamma2,gamma3,gamma4];
Lambdas = [Lambda1,Lambda2,Lambda3,Lambda4];
nus = [nu1,nu2,nu3,nu4];
epss = [eps1,eps2,eps3,eps4];

nbCompartments = nbBoxesSICR^n*nbBoxesSIIS^m;
% creating all the compartments
% a compartment is defined bcy a n-tuple (x1,x2,x3,...xn) 
% xi takes value in 1 for S,2 for I,3 for J and concerns infection i
%// Sample data
y = 1:nbBoxesSICR;%'SICR';                 %// Set of possible letters                     %// Length of each permutation
x = 1:nbBoxesSIIS;

%disease 1 is a sicr
if m==1 && n==1
    compartments = combPerso(y,x); %%manuellement ici
elseif m==2 && n==1
    compartments = combPerso(y,x,x); %%manuellement ici
elseif m==3 && n==1
    compartments = combPerso(y,x,x,x); %%manuellement ici
elseif m==3 && n==0
    compartments = combPerso(x,x,x); %%manuellement ici
end
%matrix containing all the flow rates
M = zeros(nbCompartments) - mu*eye(nbCompartments);


for j=dis %selon la maladie j, on recupere la matrice associee
    %%to simplify, we admit that there is exaclty one SICR
    if m==2 && n==0 %two SIS
        otherDiseaseStatesConstant=x';
    elseif m==1 && n==1 && j==1 %one SIJS, one SICR and the current disease is SICR
        otherDiseaseStatesConstant=x'; 
    elseif m==1 && n==1 && j==2 %one SIJS, one SICR and the current disease is SIJS
        otherDiseaseStatesConstant=y'; 
    
    elseif (n+m)>=3
        if j==1 && n==1 && m==2 %there is a sicr, other diseases are only SIJS
            otherDiseaseStatesConstant = combPerso(x,x);
        elseif j==1 && n==1 && m==3
            otherDiseaseStatesConstant = combPerso(x,x,x);
        elseif j>1 && n==1 && m==2 %sicrxsiis2
            otherDiseaseStatesConstant = combPerso(y,x);
        elseif j>1 && n==1 && m==3 %sicrxsiis3
            otherDiseaseStatesConstant = combPerso(y,x,x);
        elseif n==0 && m==3 %siis3
            otherDiseaseStatesConstant = combPerso(x,x);    
        end
    end

    if (j>1 || j==1 && n==0)    %find tuple SIJ, in a SIJS^2, there are 6 tuples :
        %(S,I1,J1),(I2,I12,I2J1),(J2,...)
        mX = ODESIJS(0,0,0,Lambdas(j),epss(j),gammas(j),nus(j),sigmas(j),b,0,0); %indice décalé de 1 si n=1
    elseif j==1 && n==1
        mX = ODESICR(0,0,0,lambdasicr, thetasicr, gammasicr,sigmasicr,b,0,0);
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
% a state Ii has the form : 1,2,1,1 or 1,1,2,1...
% this case happens iff the sum of all states is nDis+1
indexS  = sum(compartments,2)==nDis;
indexIiSIIS = (sum(compartments,2)==(nDis+1) & compartments(:,1)~=2); %Ii
indexNotIiInSIISm = compartments(:,1)==1 & ~indexIiSIIS & ~indexS; %every other compartment than Ii and S and for which HIV is susceptible (e.g.I23,J2,J3,I3J2,I2J3,J23)
indexRh  = compartments(:,1)==4 & sum(compartments,2)==nDis+(4-1); %RHIV
indexCh = compartments(:,1)==3 & sum(compartments,2)==nDis+(3-1); %Ch

%we add rho to susceptible compartment except S itself and Ii and from HIV
M(indexS,indexNotIiInSIISm) = M(indexS,indexNotIiInSIISm)+rho; %in S
M(indexNotIiInSIISm,indexNotIiInSIISm) = M(indexNotIiInSIISm,indexNotIiInSIISm) - rho*eye(3^m-m-1); %out every compartment (diagonal) 3^m : nb de compartiments siis - m :pour les Ii, -1 pour S
%and we add (1-epsi)rho to every compartment Ii (SIJS)
compIi = find(indexIiSIIS);
for ell=compIi'
    %on recupere i
    i = find(compartments(ell,:)==2);
    M(indexS,ell) = M(indexS,ell) + (1-epss(i))*rho;
    M(ell,ell) = M(ell,ell) - (1-epss(i))*rho;
end
%add rho from every compartment except S, Ii, Sympt SIJS(i.e. Ji,IiJj,Jij,Iij...)
indexHIV = compartments(:,1)~=1; %infectes par le VIH
M(indexRh, indexHIV) = M(indexRh, indexHIV) + rho;
M(indexHIV,indexHIV) = M(indexHIV,indexHIV)- eye(nbCompartments - 3^m)*rho;

% matrix to ODE system
[C,dC] = matToODE(nbCompartments,M);
dC.'

[C;compartments'].'
