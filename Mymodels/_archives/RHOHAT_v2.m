% Code which solve the utlity maximization pb in a n-diseases model (n<=4
% for now)

% one SICR : HIV
% one SEIIIS : syphilis
% one or two SEIIS : chlamydia, gonorrhea
%% Parameters
clear all;%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSEIIS=2; nSICR=1; nSEIIIS=1;
STIChoice = 'Ct'; STIChoiceMini = 'Ct';%if one STI considered
biasFactor = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 1/35;b = 5;

%SEIIS1 %chlamydia
R1     = 1.08;
sigma1 = 365/11;
nu1    = 1/1.36;
gamma1 = 365/(31.5-11); 
eps1   = 0.11;
beta1  = R1*(sigma1+mu)*(gamma1+nu1+mu)*(nu1+mu)/(sigma1*(gamma1*(1-eps1)+mu+nu1));
souslaracine = @(gamma0,beta,nu,p,sigma) (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2);
alpha1 = max((beta1*eps1*sigma1 + sqrt(souslaracine(gamma1,beta1,nu1,eps1,sigma1)))/(2*(gamma1+mu+nu1)) - (2*mu+nu1+sigma1)/2,0);

%SEIIS2 %gono
R2     = 1.08;
sigma2 = 365/5;
nu2    = 12/6;
eps2   = 0.5;
gamma2 = 365/5; 
beta2  = R2*(sigma2+mu)*(gamma2+nu2+mu)*(nu2+mu)/(sigma2*(gamma2*(1-eps2)+mu+nu2));
alpha2 = (beta2*eps2*sigma2 + sqrt(souslaracine(gamma2,beta2,nu2,eps2,sigma2)))/(2*(gamma2+mu+nu2)) - (2*mu+nu2+sigma2)/2;

%SICR 
PHIV = 0.161;
RSICR     = 1/(1-PHIV);%2.2; %P=1-1/R so R=1/(1-P);
sigmaSICR = 365/(8.2*7); %(1-0.39)*365/(8.2*7)%0.6*(); %8.2 semaines
thetaSICR = 1/7;%1/9.8;
gammaSICR = 0;%0.39*365/(8.2*7); %0.4*(); %proportion de diag au stade precoce : 40
ratioBeta = 9;
betaCSICR = RSICR*(sigmaSICR+gammaSICR+mu)*(thetaSICR+mu)/(ratioBeta*(thetaSICR+mu)+1);
betaISICR = ratioBeta*betaCSICR;
alphaSICR = betaISICR/2 - gammaSICR/2 - mu - sigmaSICR/2 - thetaSICR/2 + (betaISICR^2 - 2*betaISICR*gammaSICR - 2*betaISICR*sigmaSICR + 2*betaISICR*thetaSICR +...
    gammaSICR^2 + 2*gammaSICR*sigmaSICR - 2*gammaSICR*thetaSICR + sigmaSICR^2 - 2*sigmaSICR*thetaSICR + 4*betaCSICR*sigmaSICR + thetaSICR^2)^(1/2)/2;

%SEIIIS/S(syphilis)
PS = 0.077;
RS      = 1/(1-PS);
sigmaS  = 365/25;
tauS    = 365/45;%(1-0.55*0.31)*365/45; %0.6*
thetaS  = 12/(3.6);
gamma1S = 0;%(0.55*0.31)*365/45; %0.2*
gamma3S = 1/20;%1/5;
nuS     = 0;
betaS = RS*(thetaS+mu)*(gamma1S+tauS+mu)*(sigmaS+mu)/(sigmaS*(tauS+thetaS+mu));
Rpfun  = @(rho) (sigmaS*betaS*(tauS+thetaS+rho+mu)./((thetaS+rho+mu).*(gamma1S+rho+tauS+mu).*(sigmaS+rho+mu))-1);
alphaS = fzero(Rpfun, 0); 

%
etaS = 2;%increase of transmissbility because of HIV (>1) betaHIV^S = betaHIV*eta
omegaS = 0.5; %because of HIV in pctage (e.g. 1/sigmaHIV,S = 1/sigmaHIV*omega)

etaH = 2; %because of syphilis
omegaHIV = 0.5; %because of syphilis
%
paramSEIIS = [beta1,gamma1,nu1,sigma1,eps1,alpha1,R1;...
              beta2,gamma2,nu2,sigma2,eps2,alpha2,R2];
paramSICR = [betaISICR,betaCSICR,gammaSICR,sigmaSICR,thetaSICR,alphaSICR,RSICR,etaH,omegaHIV];
paramS    = [betaS,sigmaS,tauS,gamma1S,thetaS,gamma3S,nuS,alphaS,RS,etaS,omegaS];
vecAlphaSEIIS = [alpha1,alpha2];
vecAlpha = vecAlphaSEIIS((1:nSEIIS)+double(strcmp(STIChoice,'Ng')*nSEIIS==1)); 
if nSICR==1
    vecAlpha = [vecAlpha,alphaSICR];
else 
    paramSICR = [];
end
if nSEIIIS==1
    vecAlpha = [vecAlpha,alphaS];
else 
    paramS = [];
end
maxalpha = max(vecAlpha);

%nSICR=0;
%nSEIIIS=0;
N = nSICR+nSEIIS+nSEIIIS;


%% This section gives rhohat=argmaxU for a given model
xleft  = -1.5;
xright = 1.5;
vecC = linspace(xleft,xright,20);
tic;
%vecC = -0.25:0.0001:0.5
%vecC = [-0.2:0.01:-0.03,-0.02:0.0001:0,...
%        0.0001:0.0001:0.02,0.025:0.01:0.15]; %SICR xSEIIS Ng
%vecC = [-0.02:0.0001:0,0.01:0.01:0.15]; %SICR xSEIIS Ct (a voir)
%vecC = -1:0.0001:1; %4 dis model
%vecC = [-1:0.005:1]; %SICR xSEIIS Ct (Rc=2)

disp('attention au parametre de biais selon le modele dutilite utilisee')
%example : SISxSIS 
cond=1; i=0;
%while(cond)%interval vecC should be large enough
    %vecC = vecC(1)*1.5:step:vecC(end);
    paramSEIIS_chosen = paramSEIIS((1:nSEIIS)+double(strcmp(STIChoice,'Ng')*nSEIIS==1),:); 
    [tab,tabco,tabcn,tabTimes] = findRhohat_v2(nSEIIS,nSICR,nSEIIIS,paramSEIIS_chosen,paramSICR,paramS,mu,1,vecC,biasFactor);
    %[tab,tabco,tabcn] = findRhohat_v2(1,0,0,paramSEIIS(1,:),[],[], mu,1,vecC); 
    cond = tab.rhohat(1)~=max(vecAlpha);
    i=i+1
%end

tps = toc;


 %% cswitch for the n disease model (swith from argmax Uijk to Uij)
if N==1
    if nSICR==1
        [~,cs.c0] = U_SICR(betaISICR, betaCSICR, sigmaSICR, gammaSICR, thetaSICR, mu, 0, 0, vecAlpha,1); %c2 = c0
        [~,cs.cs1dis] = U_SICR(betaISICR, betaCSICR, sigmaSICR, gammaSICR, thetaSICR, mu, vecAlpha, 0,vecAlpha,1); %c1=c1
    end
    if nSEIIS==1 && strcmp(STIChoice,'Ct') %Ct ici
        [~,cs.c0] = U_SEIISv2(eps1, beta1, sigma1, gamma1, mu, nu1, 0, 0, vecAlpha); %c2 = c0
        [~,cs.cs1dis] = U_SEIISv2(eps1, beta1, sigma1, gamma1, mu, nu1, vecAlpha, 0,vecAlpha); %c1=c1
    elseif nSEIIS==1 && strcmp(STIChoice,'Ng') %Ng ici
        [~,cs.c0] = U_SEIISv2(eps2, beta2, sigma2, gamma2, mu, nu2, 0, 0, vecAlpha); %c2 = c0
        [~,cs.cs1dis] = U_SEIISv2(eps2, beta2, sigma2, gamma2, mu, nu2, vecAlpha, 0,vecAlpha); %c1=c1  
    end
    if nSEIIIS==1
        [~,cs.c0] = U_SEIIIS(b, betaS, sigmaS, gamma1S, gamma3S, tauS, thetaS, nuS, mu, 0, 0,vecAlpha)  ;
        [~,cs.cs1dis] = U_SEIIIS(b, betaS, sigmaS, gamma1S, gamma3S, tauS, thetaS, nuS, mu, vecAlpha, 0, vecAlpha) ; %c1=c1
    end
elseif N>=2
    cs = findThresholds_v2(nSEIIS,nSICR,nSEIIIS, tab, vecAlpha, vecC);
end
%cs = findThresholds_v2(nSEIIS,0,0, tab, vecAlpha, vecC);

if (N>=2)
    ics1 = find(cs.cs1dis==vecC);
    ics2 = find(cs.cs2dis==vecC);
    if (N>=3)
        ics3 = find(cs.cs3dis==vecC);
        if (N>=4)
            ics4 = find(cs.cs4dis==vecC);
        end
    end
end


%% plot
%save('ws_set1_4dis_v22v32.mat') %saveworkspace

%dis1 = 'Ct'; dis2 = 'HIV';dis3='Syphilis';
plot_zones_v3;
%xlim([-1,1])
%ylim([0,1.566])
%load('ws6.mat')
%%
plot(vecC,tab.two(:,2))

tab2 = tab;
clear tab;

%plot some particular combination
tab.one = tab2.one(:,[1,2]);
tab.single = tab2.single(:,[1,2]);
tab.two = tab2.two(:,1);
tab.nm = tab2.nm(:,1);
tab.rhohat = tab2.two(:,1);
N=2;
plot_zones_v2
%plot_paper1_procedure2(tab, -0.08, cs.c0,  max(tab.two(:,1)), vecC)


%% thresholds of cost for disease elimination
%rowNames = split(num2str(tab.single),' ');
%costElimTable = zeros(4^2-1,5);
costElimTable = table('Size',[4^2-1,5], 'VariableNames',{'model', 'Chlam', 'Gono', 'HIV', 'Syphilis'},...
    'VariableTypes', {'string','double','double','double','double'});

costElimTable = standardizeMissing(costElimTable,0);
%1 disease model
for i=1:4
    costElimTable{i,1} = {num2str(i)};
    costElimTable{i,i+1} = round(tabcn.one(i),3);
end

%2 disease momdel
for k=1:6
    nameModel = num2str(tab.duos(k,:));
    costElimTable{k+4,1} = {nameModel(find(~isspace(nameModel)))};
    
    for dis=tab.duos(k,:)
        j=max(find(tab.two(:,k)>=vecAlpha(dis)));
        if ~isempty(j)
            costElimTable{k+4,dis+1} = round(vecC(j),3);
        end
    end
end

%3 disease momdel
for k=1:4
    nameModel = num2str(tab.trios(k,:));
    costElimTable{k+4+6,1} = {nameModel(find(~isspace(nameModel)))};
    
    for dis=tab.trios(k,:)
        j=max(find(tab.three(:,k)>=vecAlpha(dis)));
        if ~isempty(j)
            costElimTable{k+4+6,dis+1} = vecC(j);
        end      
    end
end

%4 disease model
nameModel = num2str('1234');
costElimTable{4+6+4+1,1} = {nameModel};
for dis=1:4
    j=max(find(tab.four>=vecAlpha(dis)));
    if ~isempty(j)
        costElimTable{1+4+6+4,dis+1} = vecC(j);
    end      
end


%% Odering strategies by cost
%e.g. 1x2 + 3 + 4 ou 1x2 + 3x4

costStratTable = table('Size',[15,5], 'VariableNames',{'strategies', 'Chlam', 'Gono', 'HIV', 'Syphilis'},...
'VariableTypes', {'string','double','double','double','double'});
costStratTable.strategies = {{'1','2','3','4'}, {'1', '2', '34'}, {'1', '234'},...
             {'1', '23', '4'}, {'1', '24', '3'},...
             {'12', '3', '4'}, {'12', '34'},...
             {'13', '2', '4'}, {'13', '24'},...
             {'14', '2', '3'}, {'14', '23'},...
             {'123','4'},{'124','3'},{'134','2'},...
             {'1234'}}';
strategies = costStratTable.strategies;
 
 %costs associated to each strategy
 nbStrat = size(strategies,1);
 for i=1:nbStrat
     nbKits = size(strategies{i},2);
     for j=1:nbKits
        kitCost = costElimTable{costElimTable.model == strategies{i}{j},2:5};
        costStratTable{i,2:5} = sum([kitCost; costStratTable{i,2:5}],1,'omitnan') ;
     end
 end
 
 %detailing strategies
  for i=1:nbStrat
    costStratTable.strategiesdetailed{i} = strjoin(costStratTable.strategies{i},'+');
    costStratTable.strategiesdetailed{i}
  end
  
%% sorting by costs (HIV, syphilis, chlamydia, gonorrhea)
costStratTable_sorted = sortrows(costStratTable,[{'HIV'},{'Syphilis'},{'Chlam'},{'Gono'}],'descend');
costStratTable_sorted = costStratTable_sorted(:,[{'strategiesdetailed','HIV'},{'Syphilis'},{'Chlam'},{'Gono'}])

writetable(costStratTable_sorted,'myData.txt','Delimiter',' ')  
