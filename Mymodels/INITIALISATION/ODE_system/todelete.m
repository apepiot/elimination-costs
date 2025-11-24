clear all;

% This code generates the ODE system of the SICTPxSEIIISxSEIISxSEIIS model 
% with combined testing.
% It would need to be adapted for another purpose (e.g., generating 
% SICTPxSEIIISxSEIIS).

% calls MToODE, M_SEIIIS, M_SEIIS, M_SICTP, createTableComp

n=1; %number of diseases following the SICTP framework (i.e., HIV)
p=1; %number of diseases following the SEIIIS framework (i.e., syphilis)
m=0; %number of diseases following the SEIIS framework (i.e., Ct, Ng)

ph = 0;

syms mu b;
syms Lambdac; syms nuc; syms epsc; syms sigmac; syms gammac;
syms Lambdag; syms nug; syms epsg; syms sigmag; syms gammag;
syms Lambdah thetah sigmah ph eta_h_prep zetah
syms Lambdas sigmas taus thetas gamma1s gamma3s nus

% targeted testing rates
syms rho_h rho_s rho_c rho_g

% testing rates under PrEP
syms eta_c_prep eta_g_prep eta_s_prep

%tetsing rates under ART
syms eta_c_art eta_g_art eta_s_art

%combined testing rates of 2 diseases and more
syms rho_hs rho_hc rho_hg rho_sc rho_sg rho_cg
syms rho_hsc rho_hsg rho_hcg rho_scg
syms rho_hscg

syms VTunderART %boolean, 1 if VT can be practice under ART, 0 otherwise

nDis = m+n+p;       %number of infections in the model
dis = ["HIV","syph","Ct","Ng"];
dis = dis([n==1,p==1,m>=1,m==2]);

nbBoxesSICT   = 4;  %number of boxes in the baseline SICTP model
nbBoxesSEIIS  = 4;  %number of boxes in the baseline SEIIS model
nbBoxesSEIIIS = 5;  %number of boxes in the baseline SEIIIS model

nbCompartments = nbBoxesSICT^n*nbBoxesSEIIS^m*nbBoxesSEIIIS^p;

% creating all the compartments
% a compartment is defined by a n-tuple (x1,x2,x3,...xnDis)
no_boxesSICTP =  1:nbBoxesSICT;     %SICTP  1:S,2:I,3:C,4:P,5:Ip,6:Cp,7:T
no_boxesSEIIIS = 1:nbBoxesSEIIIS;    %SEIIIS 1:S,2:E,3:I1,4:I2,5:I3
no_boxesSEIIS =  1:nbBoxesSEIIS;     %SEIIS  1:S,2:E,3:IA,4:IS

boxesSICT   = ["S","I","C","T"]; %should be the same order than in the ODE_SICTPrEP.m
boxesSEIIIS = ["S","E","I1","I2","I3"];        %should be the same order than in the ODESEIIIS.m
boxesSEIIS  = ["S","E","IA","IS"];             %should be the same order than in the ODESEIIS.m

%Create the table of compartments
tabComp = createTableComp(m,n,p,boxesSEIIS,boxesSICT,boxesSEIIIS,dis);

tabComp.no = [1:nbCompartments]';


%% Vector pi containing new individuals rates: (pi=b, here)
syms B [nbCompartments 1]; %contains other rate, not variable dependent (e.g., pi)
B(:) = 0;
if n==1
    B(sum(table2array(tabComp(:,1:(n+m+p)))==repmat("S",1,n+m+p),2)==n+m+p) = (1-ph)*b;
    B(sum(table2array(tabComp(:,1:(n+m+p)))==["P",repmat("S",1,m+p)],2)==n+m+p) = ph*b;
else
    B(sum(table2array(tabComp(:,1:(m+p)))==repmat("S",1,n+m+p),2)==n+m+p) = b;
end

%% Fill in the matrix M corresponding to the ODE system dX/dt=MX+B

% Initialization of M
M = sym(zeros(nbCompartments) - mu*eye(nbCompartments));

%Focus on each disease progression
for INF=dis
    otherDis = dis(dis~=INF);
    if contains(INF,'HIV')
        boxesInf = boxesSICT;
        nHIV=n-1; pSYPH=p; mSTI=m;
    elseif contains(INF,'syph')
        boxesInf = boxesSEIIIS;
        nHIV=n; pSYPH = p-1; mSTI=m;
    elseif contains(INF,'Ct') || contains(INF,'Ng')
        boxesInf = boxesSEIIS;
        nHIV=n; pSYPH = p; mSTI = m-1;
    end
    
%     %find among dis\inf how many disease follow each framework
%     idx  = strfind(otherDis(:,:),'Ct','ForceCellOutput',true);   mSTI=max([find([idx{:}]),0]);
%     idx  = strfind(otherDis(:,:),'Ng','ForceCellOutput',true);   mSTI=mSTI+max([find([idx{:}]),0]);
%     idx  = strfind(otherDis(:,:),'HIV','ForceCellOutput',true);  nHIV=max([find([idx{:}]),0]);
%     idx  = strfind(otherDis(:,:),'syph','ForceCellOutput',true); pSYPH=max([find([idx{:}]),0]);
    otherDiseaseStatesConstant = createTableComp(mSTI,nHIV,pSYPH,boxesSEIIS,boxesSICT,boxesSEIIIS,otherDis);
    for i=1:size(otherDiseaseStatesConstant,1)
        %for every combination of disease states (other than INF), applies the matrix generating the ODE system of INF
        otherStates = otherDiseaseStatesConstant(i,:);
        T = innerjoin(otherStates,tabComp); %sub-table with otherStates and all the states of HIV

        %makes sure that the order of the states of INF is good (i.e., same that in the ODE system)
        [~, idx] = ismember(table2array(T(:,INF)), boxesInf');
        [~, sortorder] = sort(idx); newT = T(sortorder,:);

        %matrix of the ODE system for INF without including rho, mu, pi.
        if contains(INF,'HIV')
            Mi = M_SICT(Lambdah,thetah,sigmah,0,zetah,eta_h_prep,0,0,0);
        elseif contains(INF,'syph')
            Mi = M_SEIIIS(Lambdas,sigmas,taus,thetas,0,gamma3s,0,0,0,0);
        elseif contains(INF,'Ct')
            Mi = M_SEIIS(Lambdac,epsc,nuc,gammac,sigmac,0,0,0);
        elseif contains(INF,'Ng')
            Mi = M_SEIIS(Lambdag,epsg,nug,gammag,sigmag,0,0,0);
        end
        M(newT.no,newT.no) = M(newT.no,newT.no) + Mi;
    end
end

[X,dX,eqn,dXright] = MToODE(nbCompartments,M,B);

