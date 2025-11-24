%% This code finds the matrix and the ODE system associated to a SEIIS^mxSICR^nxSEIIIS model
% attention : in this version there is rho from SI_S to SS
clear all;

%there is only one disease SICR and all the others are SEIIS
m = 2; %nbdiseases SEIIS (max m=4)

%ne pas changer ici
n = 1; %nbdiseases SICR (VIH)
nDis = m+n;
dis = 1:nDis;
nbBoxesSICR = 4;
nbBoxesSEIIS = 4;

syms mu Lambda1 Lambda2 Lambda3 Lambda4 nu1 nu2 nu3 nu4 eps1 eps2 eps3 eps4 sigma1 sigma2 sigma3 sigma4 gamma1 gamma2 gamma3 gamma4 b rho;
syms lambdaHIV thetaHIV gammaHIV sigmaHIV

%parameters of the SEIIS
sigmas = [sigma1,sigma2,sigma3,sigma4];
gammas = [gamma1,gamma2,gamma3,gamma4];
Lambdas = [Lambda1,Lambda2,Lambda3,Lambda4];
nus = [nu1,nu2,nu3,nu4];
epss = [eps1,eps2,eps3,eps4];

nbCompartments = nbBoxesSICR^n*nbBoxesSEIIS^m;
% creating all the compartments
% a compartment is defined bcy a n-tuple (x1,x2,x3,...xn) 
% xi takes value in 1 for S,2 for E,3 for I (asymptomatic), 4 for J (symptomatic) and concerns infection i
y = 1:nbBoxesSICR;   %'SICR';                
x = 1:nbBoxesSEIIS; %1:S,2:E,3:IA,4:IS

%disease 1 is a sicr
if m==1
   compartments = combPerso(y,x);
elseif m==2
    compartments = combPerso(y,x,x);
elseif m==3
    compartments = combPerso(y,x,x,x);
elseif m==4
    compartments = combPerso(y,x,x,x,x);
end

%matrix containing all the flow rates initialization
M = zeros(nbCompartments) - mu*eye(nbCompartments);

for j=dis %selon la maladie j, on recupere la matrice associee à la progression de la maladie j
    %if we look for the progression of disease 1 (SICR)
    if j==1 && m==1
        otherDiseaseStatesConstant = x';
    elseif j==1 && m==2
        otherDiseaseStatesConstant = combPerso(x,x);
    elseif j==1 && m==3
        otherDiseaseStatesConstant = combPerso(x,x,x);
    elseif j==1 && m==4
        otherDiseaseStatesConstant = combPerso(x,x,x,x);
    end
    
    %if we look for progression of disease 2 and more (i.e. SIJS)
    if j>=2
        if  m==1
            otherDiseaseStatesConstant = y';
        elseif m==2
            otherDiseaseStatesConstant = combPerso(y,x);
        elseif m==3
            otherDiseaseStatesConstant = combPerso(y,x,x);
        elseif m==4
            otherDiseaseStatesConstant = combPerso(y,x,x,x);
        end
    end
    
    if j==1
        mX = ODESICR(0,0,0,lambdaHIV, thetaHIV, gammaHIV,sigmaHIV,b,0,0);
    else
        mX = ODESEIISv2(0,0,0,0,Lambdas(j-1),epss(j-1),nus(j-1),gammas(j-1),sigmas(j-1),b,0,0);
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

% adding the rho rate to all the compartments except those where ppl are only infected and 
% symptomatics for the STI :chlamydia, gono,etc. (e.g. Jc, Jg, RJc, RJg)
indexS  = find(sum(compartments,2)==nDis); %each disease is in state 1 (nDisx1=nDis)
indexR  = find(sum(compartments,2)==(nDis-1)+4 & compartments(:,1)==4);
%below : not infected by HIV (1) & (sympotmatic with
%the infection/s for which the indiviudal is infected (4))
%indexInfSTI = find((compartments(:,1)==1|compartments(:,1)==4) & sum(compartments(:,2:end),2)==4*m); %

%below *version _3* : once one individual is infected and symptomatic for
%one STI, s/he won't use VT but targeted testing instead based on his/her
%symptoms.
if m==0
    indexInfSTI=[];
elseif m==1
    indexInfSTI = find(compartments(:,2)==4);
elseif m==2
    indexInfSTI = find(compartments(:,2)==4 | compartments(:,3)==4);
end

for k=1:nbCompartments
    if (k~=indexS && k~=indexR && ~ismember(k,indexInfSTI)) % version 2 : if notS and if (infected by HIV or asymptomatic for another sti)
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



%% Rhohat in function of the factor f
clear all;close all; %voir param dans RHOHAT_v2
%%
vecC =  -0.3:0.01:0.3;%[-0.1:0.005:0, 0.001:0.001:0.01, 0.015:0.005:0.1];
%vecC=0.2;
figure()
%plot(vecC,alpha1*ones(1,length(vecC)),'-','LineWidth',0.5,'Color',[200/255, 200/255, 200/255])
hold on
%plot(vecC,alphaSICR*ones(1,length(vecC)),'-','LineWidth',0.1,'Color',[200/255, 200/255, 200/255])
TAB=[];CostElim=[];
hold on;
i=1;
%style = {':b','--b','-k','--r',':r','--y',':k'};
style = {':b','-k','--r',':r',':r',':k'};

for f=[0.5 1 5 10] %continuous line ?
    disp(['b=',num2str(f)])
    paramSEIIS_chosen = paramSEIIS((1:nSEIIS)+double(strcmp(STIChoice,'Ng')*nSEIIS==1),:); 
    [tab,tabco,tabcn] = findRhohat_v2(1,1,0,paramSEIIS_chosen,paramSICR,[], mu, 1,vecC,f);
    plot(vecC,tab.rhohat,style{i},'LineWidth',2,'DisplayName',['b=',num2str(f)])
    cs = findThresholds_v2(1,1,0,tab,vecAlpha,vecC);
    TAB = [TAB,tab.rhohat];
    CostElim = [CostElim;f,cs.cs2dis];
     %plot(vecC,TAB(:,i),style{i},'LineWidth',2,'DisplayName',['b=',num2str(f)])
    hold on;
    i=i+1;
end

[tabCt,~,~] = findRhohat_v2(1,0,0,paramSEIIS_chosen,[],[], mu, 1, vecC, 1);  
plot(vecC,tabCt.rhohat,'-b','LineWidth',1,'DisplayName','Ct')
[tabHIV,~,~] = findRhohat_v2(0,1,0,[],paramSICR,[], mu, 1, vecC, 1);  
hold on
plot(vecC,tabHIV.rhohat,'-r','LineWidth',1,'DisplayName','HIV')
legend
xlabel('voluntary testing cost','Interpreter','latex')
ylabel('voluntary testing rate $\hat\rho$','Interpreter','latex')
box on;
yticks([ 0 vecAlpha(cs.order(1)) 0.5 vecAlpha(cs.order(2))])
yticklabels({' ',['$\rho_{Ct}^\prime$'],'0.5',...
             ['$\rho_{HIV}^\prime$']})
legend boxoff
%save('ws_set1_HIVCt_f.mat')
% cs = findThresholds_v2(nSEIIS,nSICR,nSEIIIS, tab, vecAlpha, vecC);
% if (N>=2)
%     ics1 = find(cs.cs1dis==vecC);
%     ics2 = find(cs.cs2dis==vecC);
%     if (N>=3)
%         ics3 = find(cs.cs3dis==vecC);
%         if (N>=4)
%             ics4 = find(cs.cs4dis==vecC);
%         end
%     end
% end


%% verifications of the ODE system SICRxSEIIS
clear all
compartments = combPerso(1:4,1:4);
nbCompartments = size(compartments,1);

[betaIHIV,betaCHIV,gammaHIV,sigmaHIV,thetaHIV,~,b,mu,rho] = random_parameters(true, true);
[beta1,beta2,gamma1,gamma2,nu1,nu2,sigma1,sigma2,rho] = random_parameters(true, true);
eps1=0.2;
RHIV = (betaIHIV*(thetaHIV+mu) + betaCHIV*sigmaHIV)/((mu+thetaHIV)*(sigmaHIV+gammaHIV+mu));
R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rho))./((mu + sigma1+rho).*(gamma1 + mu + nu1).*(mu + nu1 + rho)); %ok
b/mu

tspan = 0:.1:10000;
Y0 = (b/mu)*ones(16,1)';
options = odeset('AbsTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
odefun = @(t,Y) ODE_SEIISSICR_3(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                 betaIHIV,betaCHIV,gammaHIV,sigmaHIV,thetaHIV,mu,b,rho);
[ts,Ys] = ode45(odefun,tspan,Y0,options);
T = Ys(end,:)
sum(T) - b/mu
[1:nbCompartments;T];

%verif for HIV
indexSHIV = find(compartments(:,1)==1); indexIHIV = find(compartments(:,1)==2); indexCHIV = find(compartments(:,1)==3); indexRHIV = find(compartments(:,1)==4);
[sum(T(indexSHIV)),sum(T(indexIHIV)),sum(T(indexCHIV)),sum(T(indexRHIV))]

[ts,YsHIV] = ode45(@(t,Y) ODE_SICR(t,Y,b, betaIHIV, betaCHIV, sigmaHIV, thetaHIV, gammaHIV, mu,rho,'frequency'),tspan,[1,1,1,1], options);
THIV = YsHIV(end,:) 
%if rho is different from 0, then it is not gonna be equal since VT is not 
%always performed when individual are symptomatics for STI 1 (e.g. in 
%stage IIS) whereas VT is always performed in I and C stages of HIV.

%verif for STI1
indexSSTI = find(compartments(:,2)==1); indexESTI = find(compartments(:,2)==2);
indexIASTI = find(compartments(:,2)==3); indexISSTI = find(compartments(:,2)==4);
[sum(T(indexSSTI)),sum(T(indexESTI)),sum(T(indexIASTI)),sum(T(indexISSTI))]

[ts,YsSTI] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta1,nu1,eps1,sigma1,gamma1,rho,mu),tspan,[1,1,1,1], options);
TSTI = YsSTI(end,:) 

%% verifications of the ODE system SICRxSEIIS^2
clear all
compartments = combPerso(1:4,1:4,1:4);
nbCompartments = size(compartments,1);

[betaIHIV,betaCHIV,gammaHIV,sigmaHIV,thetaHIV,~,b,mu,rho] = random_parameters(true, true);
[beta1,beta2,gamma1,gamma2,nu1,nu2,sigma1,sigma2,rho] = random_parameters(true, true);
eps1=0.2;eps2=0.5;
RHIV = (betaIHIV*(thetaHIV+mu) + betaCHIV*sigmaHIV)/((mu+thetaHIV)*(sigmaHIV+gammaHIV+mu));
R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rho))./((mu + sigma1 + rho).*(gamma1 + mu + nu1).*(mu + nu1 + rho)); %ok
R2 = (beta2*sigma2*(gamma2*(1-eps2) + mu + nu2 + eps2*rho))./((mu + sigma2 + rho).*(gamma2 + mu + nu2).*(mu + nu2 + rho)); %ok

tspan = 0:.1:1000;
Y0 = (b/mu)*ones(4*4*4,1)';
options = odeset('AbsTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
odefun = @(t,Y) ODE_SEIIS2SICR_3(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                     beta2,nu2,eps2,sigma2,gamma2,...
                                     betaIHIV,betaCHIV,gammaHIV,sigmaHIV,thetaHIV,mu,b,rho);
[ts,Ys] = ode45(odefun,tspan,Y0,options);
T = Ys(end,:)
sum(T) - b/mu
[1:nbCompartments;T];
%%
%verif for HIV
indexSHIV = find(compartments(:,1)==1); indexIHIV = find(compartments(:,1)==2); indexCHIV = find(compartments(:,1)==3); indexRHIV = find(compartments(:,1)==4);
[sum(T(indexSHIV)),sum(T(indexIHIV)),sum(T(indexCHIV)),sum(T(indexRHIV))]

[ts,YsHIV] = ode45(@(t,Y) ODE_SICR(t,Y,b, betaIHIV, betaCHIV, sigmaHIV, thetaHIV, gammaHIV, mu,rho,'frequency'),tspan,[1,1,1,1], options);
THIV = YsHIV(end,:) 
%if rho is different from 0, then it is not gonna be equal since VT is not 
%always performed when individual are symptomatics for STI 1 (e.g. in 
%stage IIS) whereas VT is always performed in I and C stages of HIV.

%verif for STI1
indexSSTI = find(compartments(:,2)==1); indexESTI = find(compartments(:,2)==2);
indexIASTI = find(compartments(:,2)==3); indexISSTI = find(compartments(:,2)==4);
[sum(T(indexSSTI)),sum(T(indexESTI)),sum(T(indexIASTI)),sum(T(indexISSTI))]

[ts,YsSTI] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta1,nu1,eps1,sigma1,gamma1,rho,mu),tspan,[1,1,1,1], options);
TSTI = YsSTI(end,:) 

%verif for STI 2
indexSSTI = find(compartments(:,3)==1); indexESTI = find(compartments(:,3)==2);
indexIASTI = find(compartments(:,3)==3); indexISSTI = find(compartments(:,3)==4);
[sum(T(indexSSTI)),sum(T(indexESTI)),sum(T(indexIASTI)),sum(T(indexISSTI))]

[ts,YsSTI] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta2,nu2,eps2,sigma2,gamma2,rho,mu),tspan,[1,1,1,1], options);
TSTI = YsSTI(end,:) 

%% plot U in function of c
clear all; close all;
ecarttropgrand=1;k=0;
while (ecarttropgrand & k<10000)

b=1/5;mu = 1/randPERT(27.2,30.6,33.7,1);

%SEIIS2 %STI random
P2     = randPERT(1,10,50,1)/100;
R2     = 1/(1-P2);
sigma2 = 365/randPERT(1,5,30,1);
nu2    = 12/randPERT(1,3,12,1);
eps2   = randPERT(0.5,0.8,1,1);
gamma2 = 365/randPERT(5,14,20,1);%365/3; 
beta2  = R2/((sigma2*(gamma2*(1-eps2) + mu + nu2))./((mu+sigma2).*(gamma2 + mu + nu2).*(mu + nu2)));
[~,~,alpha2] = Rp_SEIIS_v4(beta2,nu2,eps2,sigma2,gamma2,mu,b,0);
paramNg.beta = beta2;paramNg.sigma=sigma2;paramNg.gamma=gamma2;paramNg.nu=nu2;paramNg.eps=eps2;paramNg.alpha=alpha2;paramNg.R=R2;
paramNg.disease = 'Ct'; paramNg.modelType = 'SEIIS';

%SICR 
PHIV = randPERT(5,12,20,1)/100;
RSICR     = 1/(1-PHIV);%2.2; %P=1-1/R so R=1/(1-P);
sigmaSICR = 365/randPERT(3,7,10,1); %(1-0.39)*365/(8.2*7)%0.6*(); %8.2 semaines
thetaSICR = 1/randPERT(1,4,10,1);%1/9.8;
gammaSICR = 0; %0.39*365/(8.2*7); %0.4*(); %proportion de diag au stade precoce : 40
ratioBeta = randPERT(2,50,100,1);
betaCSICR = RSICR*(sigmaSICR+gammaSICR+mu)*(thetaSICR+mu)/(ratioBeta*(thetaSICR+mu)+1);
betaISICR = ratioBeta*betaCSICR;
[~,~,alphaSICR] = Rp_SICR_v4(betaISICR,betaCSICR,thetaSICR,sigmaSICR,gammaSICR,b,mu,0);
paramHIV.betaI = betaISICR; paramHIV.betaC = betaCSICR; paramHIV.gamma = gammaSICR; paramHIV.theta = thetaSICR;...
    paramHIV.sigma=sigmaSICR; paramHIV.alpha = alphaSICR;
paramHIV.R = RSICR; %paramHIV.eta = etaH; paramHIV.omega=omegaHIV;
paramHIV.disease = 'HIV'; paramHIV.modelType = 'SICR';
%
paramTab{1}=paramNg;
paramTab{2}=paramHIV;

ecarttropgrand = (alphaSICR-6*alpha2)<0;
k=k+1
end
paramNg.alpha
paramHIV.alpha
%%
%figure()
plotRhohat=1;
computeC0Cnn = 0;
optFindCnnAndC0.errMax = 1e-1;
optFindCnnAndC0.Tmax = 500;
optFindCnnAndC0.aPrioriCnn = -0.1;
optFindCnnAndC0.aPrioriC0 = 0.2;
optFindRhohat.sampleSize = 1000;
[bestStrategy,worstStrategy,costStratTable] = ...
    findBestStrategy(computeC0Cnn,optFindCnnAndC0,optFindRhohat,2,paramTab,mu,b,1,plotRhohat)

%%
figure
vecC = -[0.15,0.12]%-1:0.1:1
param = paramTab{1}
vecRho = 0:0.1:5
for c=vecC
U = U1_SEIISv4(param,mu,b,vecRho,c,1)
plot(vecRho,U)
hold on
end

fun = @(rho) -U1_SEIISv4(param,mu,b,rho,c,1)
min(max(fminsearch(fun,0),0),param.alpha)

