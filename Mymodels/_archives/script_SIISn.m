
%% This code finds the matrix A such that dX=Y+AX where Y=(pi,0,0,...)
clear all;
syms mu Lambda1 Lambda2 Lambda3 Lambda4 nu1 nu2 nu3 nu4 eps1 eps2 eps3 eps4 sigma1 sigma2 sigma3 sigma4 gamma1 gamma2 gamma3 gamma4 b rho;

sigmas = [sigma1,sigma2,sigma3,sigma4];
gammas = [gamma1,gamma2,gamma3,gamma4];
Lambdas = [Lambda1,Lambda2,Lambda3,Lambda4];
nus = [nu1,nu2,nu3,nu4];
epss = [eps1,eps2,eps3,eps4];

n = 2; %nbdiseases
dis = 1:n;
nbBoxes = 3;

nbCompartments = nbBoxes^n;
% creating all the compartments
% a compartment is defined bcy a n-tuple (x1,x2,x3,...xn) 
% xi takes value in 1 for S,2 for I,3 for J and concerns infection i
%// Sample data
x = 1:nbBoxes;%'SIJ';                 %// Set of possible letters                     %// Length of each permutation
%// Create all possible permutations (with repetition) of letters stored in x
C = cell(n, 1);             %// Preallocate a cell arracompartments
[C{:}] = ndgrid(x);         %// Create K grids of values
compartments = cellfun(@(x){x(:)}, C); %// Convert grids to column vectors
compartments = [compartments{:}];

%matrix containing all the flow rates
M = zeros(nbCompartments) - mu*eye(nbCompartments);

if n==2
    otherDiseaseStatesConstant=x';
elseif n>=3
    C = cell(n-1, 1);             %// Preallocate a cell arracompartments
    [C{:}] = ndgrid(x);         %// Create K grids of values
    otherDiseaseStatesConstant = cellfun(@(x){x(:)}, C); 
    otherDiseaseStatesConstant = [otherDiseaseStatesConstant{:}];
end
%find tuple SIJ, in a SIIS^2, there are 6 tuples :
%(S,I1,J1),(I2,I12,I2J1),(J2,...)
for j=dis %selon la maladie j, on recupere la matrice associee
    [mX] = ODESIIS(0,0,0,Lambdas(j),epss(j),gammas(j),nus(j),sigmas(j),b,0,0);
    
    %il faut assigner cette matrice dans M, où 1 varie de 1 à 3 et où
    %toutes les autres maladies sont dans un état constant (ex : I3,I13,I3J1)
    for k=otherDiseaseStatesConstant'
        index = find(sum(compartments(:,dis(dis~=j))==k',2)==(n-1));
        %attention les etats 1 2 3 de j doivent être sorted
        %a faire
        M(index, index) = M(index, index) + mX;
        %k
    end
end

% adding rho to all the compartments (even Ii, ? to be discussed)
% a state Ii has the form : 2,1,1 or 1,2,1...
% this case happens iff the sum of all states is nDis+1
indexIi = sum(compartments,2)==(n+1);
indexS = sum(compartments,2)==n;
%we add rho to susceptible compartment (1) except S itself
M(indexS,~indexIi & ~indexS) = M(indexS,~indexIi & ~indexS)+rho; %in S
M(~indexIi & ~indexS,~indexIi & ~indexS) = M(~indexIi & ~indexS,~indexIi & ~indexS) - rho*eye(nbCompartments-n-1); %out every compartment (diagonal)
%and we add (1-epsi)rho to every compartment Ii
compIi = find(indexIi);
for ell=compIi'
    %on recupere i
    i = find(compartments(ell,:)==2);
    M(indexS,ell) = M(indexS,ell) + (1-epss(i))*rho;
    M(ell,ell) = M(ell,ell) - (1-epss(i))*rho;
end

% matrix to ODE system
[C,dC] = matToODE_SIISn(n,M)
dC.'




