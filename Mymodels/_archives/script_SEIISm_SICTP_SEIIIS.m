%% This code generates the matrix and the ODE system associated to a SICR^n x SEIIIS^p x SEIIS^m model
% i forget to include mandatory routine testing 
clear all;

n = 1; %number of SICTPrEP (min=?, max=1)
p = 1; %number of SEIIIS (min=?,max=1)
m = 2; %number of SEIIS (min=?, max m=4)

%% Parameter initialization
syms mu b rho;
syms LambdaSTI [1 m]; syms nuSTI [1 m]; syms epsSTI [1 m]; syms sigmaSTI [1 m]; syms gammaSTI [1 m];
syms Lambdah thetah sigmah ph etah zetah
syms Lambdas sigmas taus thetas gamma1s gamma3s nus 

%(the ODE systems in the matlab functions have been created without considering the 
% simplifications below)
gamma1s = 0;
nus     = 0;
syms gammah; thetah = gammah; %the variable theta has been changed by gamma (to fit the manuscript convention)

nDis = m+n+p;       %number of infections in the model
no_dis = 1:nDis;    %numbering the infections
dis = ["HIV"+(1:n),"syph"+(1:p),"STI"+(1:m)];
nbBoxesSICTP = 7;   %number of boxes in the baseline SICTP model
nbBoxesSEIIS = 4;   %number of boxes in the baseline SEIIS model
nbBoxesSEIIIS = 5;  %number of boxes in the baseline SEIIIS model

%parameters of the SEIIS
%sigmaSTI = [sigmaSTI1,sigmaSTI2,sigmaSTI3,sigmaSTI4];
%gammaSTI = [gammaSTI1,gammaSTI2,gammaSTI3,gammaSTI4];
%LambdaSTI = [LambdaSTI1,LambdaSTI2,LambdaSTI3,LambdaSTI4];
%nuSTI = [nuSTI1,nuSTI2,nuSTI3,nuSTI4];
%epsSTI = [epsSTI1,epsSTI2,epsSTI3,epsSTI4];

nbCompartments = nbBoxesSICTP^n*nbBoxesSEIIS^m*nbBoxesSEIIIS^p;
% creating all the compartments
% a compartment is defined by a n-tuple (x1,x2,x3,...xnDis) 
no_boxesSICTP =  1:nbBoxesSICTP;     %SICTP  1:S,2:I,3:C,4:P,5:Ip,6:Cp,7:T
no_boxesSEIIIS = 1:nbBoxesSEIIIS;    %SEIIIS 1:S,2:E,3:I1,4:I2,5:I3
no_boxesSEIIS =  1:nbBoxesSEIIS;     %SEIIS  1:S,2:E,3:IA,4:IS

boxesSICTP  = ["S","I","C","P","Ip","Cp","T"]; %should be the same order than in the ODESICTP.m
boxesSEIIIS = ["S","E","I1","I2","I3"];        %should be the same order than in the ODESEIIIS.m
boxesSEIIS  = ["S","E","IA","IS"];             %should be the same order than in the ODESEIIS.m

%Create the table of compartments
tabComp = createTableComp(m,n,p,boxesSEIIS,boxesSICTP,boxesSEIIIS,dis);
tabComp.no = [1:nbCompartments]';

%% Fill in the matrix M corresponding to the ODE system dX/dt=MX+B

%Initialization
M = sym(zeros(nbCompartments) - mu*eye(nbCompartments));

if m+n+p==1
    if (n==1 && (m+p)==0)
        %M given by the ODE system of SICTP
    elseif (p==1 && (n+p)==0)
        %M given by the ODE system of SEIIIS
    elseif (m==1 && (p+n)==0)
        %M given by the ODE system of SEIIS
    end
else
    %Focus on each disease progression
    for inf=dis
        if contains(inf,'HIV')
            boxesInf = boxesSICTP;
        elseif contains(inf,'syph')
            boxesInf = boxesSEIIIS;
        elseif contains(inf,'STI')
            boxesInf = boxesSEIIS;
        end

        otherDis = dis(dis~=inf);
        idx  = strfind(otherDis(:,:),'STI','ForceCellOutput',true); mSTI=max([find([idx{:}]),0]);
        idx  = strfind(otherDis(:,:),'HIV','ForceCellOutput',true); nHIV=max([find([idx{:}]),0]);
        idx  = strfind(otherDis(:,:),'syph','ForceCellOutput',true); pSYPH=max([find([idx{:}]),0]);
        otherDiseaseStatesConstant = createTableComp(mSTI,nHIV,pSYPH,boxesSEIIS,boxesSICTP,boxesSEIIIS,otherDis);
        for i=1:size(otherDiseaseStatesConstant,1)
            %for every combination of disease states (other than HIV), applies the matrix
            %of the ODE system of HIV
            otherStates = otherDiseaseStatesConstant(i,:);
            T = innerjoin(otherStates,tabComp); %subtable with otherStates and all the states of HIV

            %makes sure that the order of the states of inf is good (i.e. same that in the ODE system)
            [found, idx] = ismember(table2array(T(:,inf)), boxesInf');
            [~, sortorder] = sort(idx); newT = T(sortorder,:);

            %matrix of the ODE system for HIV
            %we don't include: rho, mu, b.
            if contains(inf,'HIV')
                Mh = M_SICTP(Lambdah,thetah,sigmah,ph,zetah,etah,0,0,0);
            elseif contains(inf,'syph')
                Mh = M_SEIIIS(Lambdas,sigmas,taus,thetas,gamma1s,gamma3s,nus,0,0,0);
            elseif contains(inf,'STI1')
                Mh = M_SEIIS(LambdaSTI1,epsSTI1,nuSTI1,gammaSTI1,sigmaSTI1,0,0,0);
            elseif contains(inf,'STI2')
                Mh = M_SEIIS(LambdaSTI2,epsSTI2,nuSTI2,gammaSTI2,sigmaSTI2,0,0,0);
            elseif contains(inf,'STI3')
                Mh = M_SEIIS(LambdaSTI3,epsSTI3,nuSTI3,gammaSTI3,sigmaSTI3,0,0,0);
            elseif contains(inf,'STI4')
                Mh = M_SEIIS(LambdaSTI4,epsSTI4,nuSTI4,gammaSTI4,sigmaSTI4,0,0,0);
            end
            M(newT.no,newT.no) = M(newT.no,newT.no) + Mh;
        end
    end
end

%% Adds additionnal rate: rho (combined voluntary testing) rates and pi (input)
syms B [nbCompartments 1]; %contains other rate, not variable dependent (e.g. pi)
B(:) = 0;

%Adds input parameters
if n==0 %HIV is not included in the model
    B(sum(table2array(tabComp(:,1:(n+m+p)))==repmat("S",1,n+m+p),2)==n+m+p) = b;
elseif n==1
    B(sum(table2array(tabComp(:,1:(n+m+p)))==repmat("S",1,n+m+p),2)==n+m+p) = (1-ph)*b;
    B(sum(table2array(tabComp(:,1:(n+m+p)))==["P",repmat("S",1,m+p)],2)==n+m+p) = ph*b;
end

%% Adds voluntary testing rate rho
%Symptomatic individuals for STIs test because of symptoms and not with
%voluntary testing, s/he will use a combined testing test. 
if n==0 %HIV is not included in the model, then all testing through combined testing go to SSSS... 
    %1. Compartments where individuals are symptomatic (IS) for at least one STI
    symptomatics  = tabComp(sum(table2array(tabComp(:,1:n+m+p))=="IS",2)>=1,:);
    
    %2. Compartments of asymptomatics and infected by one or several infections
    asymptomatics = tabComp(sum(table2array(tabComp(:,1:n+m+p))=="IS",2)==0,:); %S, E, I1, I2, I3 or IA
    asympt_inf    = asymptomatics(sum(table2array(asymptomatics(:,1:n+m+p))==repmat("S",1,n+m+p),2)<n+m+p,:);
    susceptible   = asymptomatics(sum(table2array(asymptomatics(:,1:n+m+p))==repmat("S",1,n+m+p),2)==n+m+p,:);
    
    if size(symptomatics,1)+size(asympt_inf,1)+size(susceptible,1)~=nbCompartments
        error('le nombre de compartiments ne correspond pas')
    end
    
    %3. Adds voluntary testing rate rho to each asympt_inf compartments to the susceptible compartment S
    M(susceptible.no, asympt_inf.no) = M(susceptible.no, asympt_inf.no)+rho;
    for i=asympt_inf.no
        M(i,i) = M(i,i) - rho;
    end    
elseif n==1 %HIV is in the model, then all testing through combined testing 
            %go to SSSS, PSSS or TSSS
    compSSS = tabComp(sum(table2array(tabComp(:,1:n+m+p))==repmat("S",1,n+m+p),2)==n+m+p,:);
    compPSS = tabComp(sum(table2array(tabComp(:,1:n+m+p))==["P",repmat("S",1,m+p)],2)==n+m+p,:);
    compTSS = tabComp(sum(table2array(tabComp(:,1:n+m+p))==["T",repmat("S",1,m+p)],2)==n+m+p,:);
    
    %si : not SSS et TSS et not PSS
    %%si : not symptomatic (IS) pour aucune des maladies
    %%%elseif state(HIV)=S, alors va dans SSS
    %%%elseif state(HIV)=P, alors va dans PSS
    %%else: alors va dans TSS
    for k=1:nbCompartments
        currentComp = tabComp(k,:);
        if sum(currentComp.no==[compSSS.no,compPSS.no,compTSS.no])==0
            if sum(table2array(currentComp(:,1:n+m+p))=="IS")==0 %not symptomatic
                if currentComp.HIV1=="S" %goes to S
                    M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-rho;
                    M(compSSS.no,currentComp.no)     = M(compSSS.no,currentComp.no)+rho;
                    disp(join(['goes to S:',table2array(currentComp(:,1:n+m+p))]))
                elseif currentComp.HIV1=="P" %goes to P
                    M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-rho;
                    M(compPSS.no,currentComp.no)     = M(compPSS.no,currentComp.no)+rho;
                    disp(join(['goes to P:',table2array(currentComp(:,1:n+m+p))]))
                else %goes to T
                    M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-rho;
                    M(compTSS.no,currentComp.no)     = M(compTSS.no,currentComp.no)+rho;
                    disp(join(['goes to T:',table2array(currentComp(:,1:n+m+p))]))
                end
            else
                %not voluntary testing
                disp(join(['no VT:',table2array(currentComp(:,1:n+m+p))]))
            end
        end
    end
end

%% Converting the matrix product to ODE system
[C,dC,eqn] = matToODE_v2(nbCompartments,M,B);
eqn.'
tabComp.X = C;


%% Verifications of the calculations
clearvars -except tabComp; close all;
% verifier que l'ES vaut b/mu au total
% faire des verifs: avec rho=0, on doit retomber sur les sous modeles
% faire des tests avec Lambda(rho) et ce qu'on obtient en sommant 
b=2;
[paramTab,mu,~] = sampleParameters_v3(true,true,true,true,b);
paramTab{3}.p = 0.4; paramTab{3}.zeta=0.3;paramTab{3}.eta=4; rho=0.;
%check if R's are bigger than one
[RS,LambdaS,~] = Rp_SEIIIS_v4(paramTab{4}.beta,paramTab{4}.sigma,...
    paramTab{4}.tau,paramTab{4}.nu,paramTab{4}.gamma1,paramTab{4}.theta,...
    paramTab{4}.gamma3,mu,b,rho);
[RC,LambdaC,~] = Rp_SEIIS_v4(paramTab{1}.beta,paramTab{1}.nu,paramTab{1}.eps,...
    paramTab{1}.sigma,paramTab{1}.gamma,mu,b,rho);
[RN,LambdaN,~] = Rp_SEIIS_v4(paramTab{2}.beta,paramTab{2}.nu,paramTab{2}.eps,...
    paramTab{2}.sigma,paramTab{2}.gamma,mu,b,rho);
[RH,LambdaH,~] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,paramTab{3}.theta,...
    paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,paramTab{3}.p,mu,b,rho);

Y0 = ones(size(tabComp,1),1); tspan = 0:1:500;
options = odeset('RelTol',1e-5,'Stats','on')%,'OutputFcn',@odeplot);
if (RH>1 || RC>1 || RN>1 || RH>1)   
    tic
%     [ts,Ys] = ode45(@(t,Y) ODE_SICTPSEIIISSEIIS2(t,Y,paramTab{3}.betaI,paramTab{3}.betaC,...
%         paramTab{3}.sigma,paramTab{3}.theta,paramTab{3}.zeta,paramTab{3}.eta,paramTab{3}.p,...
%                             paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.gamma3,paramTab{4}.tau,paramTab{4}.theta,...
%                             paramTab{1}.beta,paramTab{1}.gamma,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,...
%                             paramTab{2}.beta,paramTab{2}.gamma,paramTab{2}.nu,paramTab{2}.eps,paramTab{2}.sigma,...
%                             tabComp,mu,b,rho),...
%         tspan,Y0,options);
    [ts,Ys] = ode45(@(t,Y) ODE_SICTPSEIIISSEIIS(t,Y,paramTab{3}.betaI,paramTab{3}.betaC,...
            paramTab{3}.sigma,paramTab{3}.theta,paramTab{3}.zeta,paramTab{3}.eta,paramTab{3}.p,...
            paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.gamma3,paramTab{4}.tau,paramTab{4}.theta,...
            paramTab{1}.beta,paramTab{1}.gamma,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,...
        tabComp,mu,b,rho),...
        tspan,Y0,options);
    toc
    popEnd = Ys(end,:);
    Nnum   = sum(popEnd)
end

%% Comparing to the single-disease models
N=b/mu

Snum = sum(popEnd(table2array(tabComp(tabComp.HIV1=="S","no"))))
SES_H = b*(1-paramTab{3}.p)./(LambdaH+mu)
Pnunm = sum(popEnd(table2array(tabComp(tabComp.HIV1=="P","no"))))
PES_H = paramTab{3}.p*b./(LambdaH*(1-paramTab{3}.zeta)+mu);
PrevHIV = sum(popEnd(table2array(tabComp(tabComp.HIV1=="I" | tabComp.HIV1=="Ip" | tabComp.HIV1=="C" | tabComp.HIV1=="Cp" |...
    tabComp.HIV1=="T" ,"no"))))./N
1-(SES_H+PES_H)*mu/b

%Inum = sum(popEnd(table2array(tabComp(tabComp.HIV1=="I","no"))))
%Lambdah*Snum./(paramTab{3}.sigma+rho+mu)
%b*(1-paramTab{3}.p)*Lambdah./((Lambdah+mu)*(paramTab{3}.sigma+rho+mu))

Y=popEnd;
Lambdah = paramTab{3}.betaI*sum(Y(table2array(tabComp(tabComp.HIV1=="I" | tabComp.HIV1=="Ip","no"))))./N +...
        paramTab{3}.betaC*sum(Y(table2array(tabComp(tabComp.HIV1=="C" | tabComp.HIV1=="Cp","no"))))./N
LambdaH

PrevS = sum(popEnd(table2array(tabComp(tabComp.syph1=="E" | tabComp.syph1=="I1" | tabComp.syph1=="I2" | tabComp.syph1=="I3" ,"no"))))./N
1-1./RS

PrevC = sum(popEnd(table2array(tabComp(tabComp.STI1=="IA" | tabComp.STI1=="IS" | tabComp.STI1=="E" ,"no"))))./N
1-1./RC

PrevN = sum(popEnd(table2array(tabComp(tabComp.STI2=="IA" | tabComp.STI2=="IS" | tabComp.STI2=="E" ,"no"))))./N
1-1./RN



%%
% 
% tspan = 0:.1:1000;
% Y0 = (b/mu)*ones(80,1)';
% options = odeset('AbsTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
% % odefun = @(t,Y) ODE_SICRSEIIISSEIIS_3(t,Y,betaIHIV,betaCHIV,thetah,gammah,sigmah,...
% %                                       betaSyph,sigmas,taus,thetas,gamma1s,gamma3s,nus,...
% %                                       beta1,nuSTI1,epsSTI1,sigmaSTI1,gammaSTI1,mu,b,rho,indexInfSyph,indexInfSEIIS,indexIHIV,indexCHIV);
% odefun = @(t,Y) ODE_SEIISSICRSEIIIS_3(t,Y,beta1,gammaSTI1,nuSTI1,sigmaSTI1,epsSTI1,...
%                                      betaIHIV,betaCHIV,gammah,sigmah,thetah,...
%                                      betaSyph,sigmas,taus,nus,gamma1s,thetas,gamma3s,mu,b,rho);
% [ts,Ys] = ode45(odefun,tspan,Y0,options);
% T = Ys(end,:)
% sum(T)
% [1:nbCompartments;T]
% 
% %% verifications
% %HIV
% indexSHIV = find(compartments(:,1)==1); indexIHIV = find(compartments(:,1)==2); indexCHIV = find(compartments(:,1)==3); indexRHIV = find(compartments(:,1)==4);
% [sum(T(indexSHIV)),sum(T(indexIHIV)),sum(T(indexCHIV)),sum(T(indexRHIV))]
% 
% [ts,YsHIV] = ode45(@(t,Y) ODE_SICR(t,Y,b, betaIHIV, betaCHIV, sigmah, thetah, gammah, mu,rho,'frequency'),tspan,[1,1,1,1], options);
% THIV = YsHIV(end,:)
% 
% %STI 1
% indexS1 = find(compartments(:,3)==1); indexE1 = find(compartments(:,3)==2); indexI1 = find(compartments(:,3)==3); indexJ1 = find(compartments(:,3)==4);
% [sum(T(indexS1)),sum(T(indexE1)),sum(T(indexI1)),sum(T(indexJ1))]
% 
% R1 = (beta1*sigmaSTI1*(gammaSTI1*(1-epsSTI1) + mu + nuSTI1 + epsSTI1*rho))./((mu + sigmaSTI1+rho).*(gammaSTI1 + mu + nuSTI1).*(mu + nuSTI1 + rho)); %ok
% S1  = b/(mu*max(R1,1));
% lambda1 = max(beta1*(R1-1)*(sigmaSTI1+rho+mu)/(beta1+(sigmaSTI1+rho+mu)*R1),0);
% E1  = lambda1*S1/(sigmaSTI1+mu+rho);
% I1  = (1-epsSTI1)*sigmaSTI1*E1/(nuSTI1+rho+mu);
% J1  = epsSTI1*sigmaSTI1*E1/(gammaSTI1+nuSTI1+mu);
% [S1,E1,I1,J1]
% [ts,Ys1] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta1,nuSTI1,epsSTI1,sigmaSTI1,gammaSTI1,rho,mu),tspan,[1,1,1,1], options);
% Ys1(end,:)
% 
% %syphilis
% indexSs = find(compartments(:,2)==1); indexEs = find(compartments(:,2)==2);
% indexI1s = find(compartments(:,2)==3); indexI2s = find(compartments(:,2)==4);
% indexI3s = find(compartments(:,2)==5);
% [sum(T(indexSs)),sum(T(indexEs)),sum(T(indexI1s)),sum(T(indexI2s)),sum(T(indexI3s))]
% RSyph = sigmas*betaSyph*(taus+thetas+rho+mu)/((thetas+rho+mu)*(gamma1s+rho+taus+mu)*(sigmas+mu+rho));
% 
% numerateur = -(gamma3s+mu+nus+rho)*((mu+rho+sigmas)*(mu+rho+thetas)*(gamma1s+mu+rho+taus)-betaSyph*sigmas*(rho+taus+thetas+mu));
% denominateur = gamma1s*mu^2 + gamma3s*mu^2 + gamma1s*rho^2 + gamma3s*rho^2 + mu^2*nus + 3*mu*rho^2 + 3*mu^2*rho + nus*rho^2 + mu^2*sigmas + mu^2*taus + mu^2*thetas + rho^2*sigmas + rho^2*taus +...
%     rho^2*thetas + mu^3 + rho^3 + gamma1s*gamma3s*mu + gamma1s*gamma3s*rho + gamma1s*gamma3s*thetas +...
%     gamma1s*mu*nus + 2*gamma1s*mu*rho + 2*gamma3s*mu*rho + gamma1s*nus*rho + gamma3s*mu*sigmas +...
%     gamma3s*mu*taus + gamma1s*mu*thetas + gamma3s*mu*thetas + gamma1s*nus*thetas + gamma3s*rho*sigmas +...
%     gamma3s*rho*taus + gamma1s*rho*thetas + gamma3s*rho*thetas + 2*mu*nus*rho + gamma3s*sigmas*taus +...
%     gamma3s*sigmas*thetas + mu*nus*sigmas + gamma3s*taus*thetas + mu*nus*taus + mu*nus*thetas +...
%     2*mu*rho*sigmas + 2*mu*rho*taus + 2*mu*rho*thetas + nus*rho*sigmas + nus*rho*taus + mu*sigmas*taus + ...
%     nus*rho*thetas + mu*sigmas*thetas + mu*taus*thetas + nus*sigmas*taus + nus*sigmas*thetas +...
%     nus*taus*thetas + rho*sigmas*taus + rho*sigmas*thetas + rho*taus*thetas + sigmas*taus*thetas;
% lambdaS = max(numerateur/denominateur,0);
% 
% Ss  = b/(mu*max(RSyph,1)); Es = b*lambdaS/(mu*RSyph*(sigmas+rho+mu));
% I1s = sigmas*Es/(taus+gamma1s+rho+mu);
% I2s = taus*I1s/(thetas+rho+mu);
% I3s = thetas*I2s/(nus+gamma3s+rho+mu);
% [Ss,Es,I1s,I2s,I3s]
% [ts,YSyph] = ode45(@(t,Y) SEIIIS(t,Y,b,betaSyph,sigmas, gamma1s, gamma3s, taus, thetas, nus, mu,rho),tspan,[1,1,1,1,1], options);
% YSyph(end,:)
% 
% %% SICRxSEIIISxSEIIS^2 (2) (versions_2)
% clear all
% compartments = combPerso(1:4,1:5,1:4,1:4);
% indexCHIV = find(compartments(:,1)==3);
% indexIHIV = find(compartments(:,1)==2);
% indexInfSyph = find(compartments(:,2)==3 | compartments(:,2)==4);
% indexInf1 = find(compartments(:,3)==3 | compartments(:,3)==4);
% indexInf2 = find(compartments(:,4)==3 | compartments(:,4)==4);
% 
% nbCompartments = size(compartments,1);
% 
% [betaIHIV,betaCHIV,gammah,sigmah,thetah,~,b,mu,rho] = random_parameters(true, true);
% [betaSyph,sigmas,gamma1s,gamma3s,taus,thetas,nus,mu,rho] = random_parameters(true, true);
% [beta1,beta2,gammaSTI1,gammaSTI2,nuSTI1,nuSTI2,sigmaSTI1,sigmaSTI2,rho] = random_parameters(true, true);
% epsSTI1=0.2;epsSTI2=0.3;
% RHIV = (betaIHIV*(thetah+mu) + betaCHIV*sigmah)/((mu+thetah)*(sigmah+gammah+mu));
% R1 = (beta1*sigmaSTI1*(gammaSTI1*(1-epsSTI1) + mu + nuSTI1 + epsSTI1*rho))./((mu + sigmaSTI1+rho).*(gammaSTI1 + mu + nuSTI1).*(mu + nuSTI1 + rho)); %ok
% R2 = (beta2*sigmaSTI2*(gammaSTI2*(1-epsSTI2) + mu + nuSTI2 + epsSTI2*rho))./((mu + sigmaSTI2+rho).*(gammaSTI2 + mu + nuSTI2).*(mu + nuSTI2 + rho)); %ok
% RSyph = sigmas*betaSyph*(taus+thetas+rho)/((thetas+rho)*(gamma1s+rho+taus+mu)*(sigmas+mu+rho));
% b/mu
% 
% %rho=0;
% 
% tspan = 0:.1:1000;
% Y0 = (b/mu)/80*ones(320,1)';
% options = odeset('AbsTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
% % odefun = @(t,Y) ODE_SICRSEIIISSEIIS2_3(t,Y,betaIHIV,betaCHIV,thetah,gammah,sigmah,...
% %                                       betaSyph,sigmas,taus,thetas,gamma1s,gamma3s,nus,...
% %                                       beta1,nuSTI1,epsSTI1,sigmaSTI1,gammaSTI1,...
% %                                       beta2,nuSTI2,epsSTI2,sigmaSTI2,gammaSTI2,mu,b,rho,...
% %                                       indexInfSyph,indexInf1,indexInf2,indexIHIV,indexCHIV);
% 
% odefun = @(t,Y) ODE_SEIIS2SICRSEIIIS_3(t,Y,beta1,gammaSTI1,nuSTI1,sigmaSTI1,epsSTI1,...
%                                          beta2,gammaSTI2,nuSTI2,sigmaSTI2,epsSTI2,...
%                                          betaIHIV,betaCHIV,gammah,sigmah,thetah,...
%                                          betaSyph,sigmas,taus,gamma1s,thetas,gamma3s,nus,...
%                                          mu,b,rho);
% [ts,Ys] = ode45(odefun,tspan,Y0,options);
% T = Ys(end,:)
% sum(T) - b/mu
% [1:nbCompartments;T] 
% 
% %% verifications
% 
% indexSHIV = find(compartments(:,1)==1); indexIHIV = find(compartments(:,1)==2); indexCHIV = find(compartments(:,1)==3); indexRHIV = find(compartments(:,1)==4);
% [sum(T(indexSHIV)),sum(T(indexIHIV)),sum(T(indexCHIV)),sum(T(indexRHIV))]
% [ts,YsHIV] = ode45(@(t,Y) ODE_SICR(t,Y,b, betaIHIV, betaCHIV, sigmah, thetah, gammah, mu,rho,'frequency'),tspan,[1,1,1,1], options);
% THIV = YsHIV(end,:)
% 
% %STI 1
% indexS1 = find(compartments(:,3)==1); indexE1 = find(compartments(:,3)==2); indexI1 = find(compartments(:,3)==3); indexJ1 = find(compartments(:,3)==4);
% [sum(T(indexS1)),sum(T(indexE1)),sum(T(indexI1)),sum(T(indexJ1))]
% S1  = b/(mu*max(R1,1));
% lambda1 = max(beta1*(R1-1)*(sigmaSTI1+rho+mu)/(beta1+(sigmaSTI1+rho+mu)*R1),0);
% E1  = lambda1*S1/(sigmaSTI1+mu+rho);
% I1  = (1-epsSTI1)*sigmaSTI1*E1/(nuSTI1+rho+mu);
% J1  = epsSTI1*sigmaSTI1*E1/(gammaSTI1+nuSTI1+mu);
% [S1,E1,I1,J1]
% [ts,Ys1] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta1,nuSTI1,epsSTI1,sigmaSTI1,gammaSTI1,rho,mu),tspan,[1,1,1,1], options);
% Ys1(end,:)
% 
% %STI 2
% indexS2 = find(compartments(:,4)==1); indexE2 = find(compartments(:,4)==2); indexI2 = find(compartments(:,4)==3); indexJ2 = find(compartments(:,4)==4);
% [sum(T(indexS2)),sum(T(indexE2)),sum(T(indexI2)),sum(T(indexJ2))]
% S2  = b/(mu*max(R2,1));
% lambda2 = max(beta2*(R2-1)*(sigmaSTI2+rho+mu)/(beta2+(sigmaSTI2+rho+mu)*R2),0);
% E2  = lambda2*S2/(sigmaSTI2+mu+rho);
% I2  = (1-epsSTI2)*sigmaSTI2*E2/(nuSTI2+rho+mu);
% J2  = epsSTI2*sigmaSTI2*E2/(gammaSTI2+nuSTI2+mu);
% [S2,E2,I2,J2]
% [ts,Ys2] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta2,nuSTI2,epsSTI2,sigmaSTI2,gammaSTI2,rho,mu),tspan,[1,1,1,1], options);
% Ys2(end,:)
% 
% %syphilis
% indexSs = find(compartments(:,2)==1); indexEs = find(compartments(:,2)==2);
% indexI1s = find(compartments(:,2)==3); indexI2s = find(compartments(:,2)==4);
% indexI3s = find(compartments(:,2)==5);
% [sum(T(indexSs)),sum(T(indexEs)),sum(T(indexI1s)),sum(T(indexI2s)),sum(T(indexI3s))]
% RSyph = sigmas*betaSyph*(taus+thetas+rho+mu)/((thetas+rho+mu)*(gamma1s+rho+taus+mu)*(sigmas+mu+rho));
% 
% numerateur = -(gamma3s+mu+nus+rho)*((mu+rho+sigmas)*(mu+rho+thetas)*(gamma1s+mu+rho+taus)-betaSyph*sigmas*(rho+taus+thetas+mu));
% denominateur = gamma1s*mu^2 + gamma3s*mu^2 + gamma1s*rho^2 + gamma3s*rho^2 + mu^2*nus + 3*mu*rho^2 + 3*mu^2*rho + nus*rho^2 + mu^2*sigmas + mu^2*taus + mu^2*thetas + rho^2*sigmas + rho^2*taus +...
%     rho^2*thetas + mu^3 + rho^3 + gamma1s*gamma3s*mu + gamma1s*gamma3s*rho + gamma1s*gamma3s*thetas +...
%     gamma1s*mu*nus + 2*gamma1s*mu*rho + 2*gamma3s*mu*rho + gamma1s*nus*rho + gamma3s*mu*sigmas +...
%     gamma3s*mu*taus + gamma1s*mu*thetas + gamma3s*mu*thetas + gamma1s*nus*thetas + gamma3s*rho*sigmas +...
%     gamma3s*rho*taus + gamma1s*rho*thetas + gamma3s*rho*thetas + 2*mu*nus*rho + gamma3s*sigmas*taus +...
%     gamma3s*sigmas*thetas + mu*nus*sigmas + gamma3s*taus*thetas + mu*nus*taus + mu*nus*thetas +...
%     2*mu*rho*sigmas + 2*mu*rho*taus + 2*mu*rho*thetas + nus*rho*sigmas + nus*rho*taus + mu*sigmas*taus + ...
%     nus*rho*thetas + mu*sigmas*thetas + mu*taus*thetas + nus*sigmas*taus + nus*sigmas*thetas +...
%     nus*taus*thetas + rho*sigmas*taus + rho*sigmas*thetas + rho*taus*thetas + sigmas*taus*thetas;
% lambdaS = max(numerateur/denominateur,0);
% 
% Ss  = b/(mu*max(RSyph,1)); Es = b*lambdaS/(mu*RSyph*(sigmas+rho+mu));
% I1s = sigmas*Es/(taus+gamma1s+rho+mu);
% I2s = taus*I1s/(thetas+rho+mu);
% I3s = thetas*I2s/(nus+gamma3s+rho+mu);
% [Ss,Es,I1s,I2s,I3s]
% [ts,YSyph] = ode45(@(t,Y) SEIIIS(t,Y,b,betaSyph,sigmas, gamma1s, gamma3s, taus, thetas, nus, mu,rho),tspan,[1,1,1,1,1], options);
% YSyph(end,:)
