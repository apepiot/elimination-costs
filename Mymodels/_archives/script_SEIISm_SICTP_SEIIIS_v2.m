%% This code generates the matrix and the ODE system associated to a SICR^n x SEIIIS^p x SEIIS^m model
% combined and targeted testing are mixed
clear all;

%% Parameter initialization
% Nothing to change below, this code generates a 4-disease model.
n = 1; %number of SICTPrEP (min=?, max=1)
p = 1; %number of SEIIIS (min=?,max=1)
m = 2; %number of SEIIS (min=?, max m=4)

syms mu b;
syms LambdaCt; syms nuCt; syms epsCt; syms sigmaCt; syms gammaCt;
syms LambdaNg; syms nuNg; syms epsNg; syms sigmaNg; syms gammaNg;

syms Lambdah thetah sigmah ph etah zetah
syms Lambdas sigmas taus thetas gamma1s gamma3s nus
% targeted testing rate:
syms rho_h rho_s rho_c rho_g
syms eta_c eta_g eta_s
syms eta_c_art eta_g_art eta_s_art
%combined testing of 2 infections
syms rho_hs rho_hc rho_hg
syms rho_sc rho_sg
syms rho_cg
syms rho_hsc rho_hsg rho_hcg rho_scg
syms rho_hscg

%(the ODE systems in the matlab functions have been created without considering the
% simplifications below)
gamma1s = 0;
nus     = 0;

nDis = m+n+p;       %number of infections in the model
no_dis = 1:nDis;    %numbering the infections
dis = ["HIV","syph","Ct","Ng"];
%dis = ["HIV","Ct"];
nbBoxesSICTP  = 7;  %number of boxes in the baseline SICTP model
nbBoxesSEIIS  = 4;  %number of boxes in the baseline SEIIS model
nbBoxesSEIIIS = 5;  %number of boxes in the baseline SEIIIS model

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
        elseif contains(inf,'Ct') || contains(inf,'Ng')
            boxesInf = boxesSEIIS;
        end
        
        otherDis = dis(dis~=inf);
        idx  = strfind(otherDis(:,:),'Ct','ForceCellOutput',true); mSTI=max([find([idx{:}]),0]);
        idx  = strfind(otherDis(:,:),'Ng','ForceCellOutput',true); mSTI=mSTI+max([find([idx{:}]),0]);
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
                Mh = M_SICTP(Lambdah,thetah,sigmah,ph,zetah,etah,0,rho_h,0);
            elseif contains(inf,'syph')
                Mh = M_SEIIIS(Lambdas,sigmas,taus,thetas,gamma1s,gamma3s,nus,0,rho_s,0);
            elseif contains(inf,'Ct')
                Mh = M_SEIIS(LambdaCt,epsCt,nuCt,gammaCt,sigmaCt,0,rho_c,0);
            elseif contains(inf,'Ng')
                Mh = M_SEIIS(LambdaNg,epsNg,nuNg,gammaNg,sigmaNg,0,rho_g,0);
            end
            M(newT.no,newT.no) = M(newT.no,newT.no) + Mh;
        end
    end
end

%% Adds additionnal rate: pi (input)
syms B [nbCompartments 1]; %contains other rate, not variable dependent (e.g. pi)
B(:) = 0;

%Adds input parameters
B(sum(table2array(tabComp(:,1:(n+m+p)))==repmat("S",1,n+m+p),2)==n+m+p) = (1-ph)*b;
B(sum(table2array(tabComp(:,1:(n+m+p)))==["P",repmat("S",1,m+p)],2)==n+m+p) = ph*b;

%% To adapt the code for a k-disease model, k<n+p+m
% i add columns for the other infections (the ones not considered in the model) and make it susceptible 
if m==1 && n==1 && p==0
    tabComp.Ng = repmat("S",nbCompartments,1);
    tabComp.syph = repmat("S",nbCompartments,1);
end

%% Adds combined voluntary testing rate rho_xxx
%Symptomatic individuals for STIs test because of symptoms and not with
%voluntary testing.


% For all combination of infections, we need to add combined voluntary
% testing rate
% e.g. for a given box *I(H)E(s)IA(C)E(G)*
% HIV and syph tested together, then combined testing leads to S(H)S(syph),
% Ct tested alone, targeted testing goes to S(C) (already accounted in the matrix M)
% Ng ---------------------------------------S(G) (idem)

CVTcomb = table({["HIV","syph"];["HIV","Ct"];["HIV","Ng"];["syph","Ct"];["syph","Ng"];["Ct","Ng"];...
    ["HIV","syph","Ct"];["HIV","syph","Ng"];["HIV","Ct","Ng"];["syph","Ct","Ng"];...
    ["HIV","syph","Ct","Ng"]},...
    [1;1;1;0;0;0;1;1;1;0;1],...
    [1;0;0;1;1;0;1;1;0;1;1],...
    [0;1;0;1;0;1;1;0;1;1;1],...
    [0;0;1;0;1;1;0;1;1;1;1],...
    [rho_hs;rho_hc;rho_hg;rho_sc;rho_sg;rho_cg;rho_hsc;rho_hsg;rho_hcg;rho_scg;rho_hscg],...
    'VariableNames',{'kit','isHIV','issyph','isCt','isNg','rho'});

CVTcombHIV = CVTcomb(logical(CVTcomb.isHIV),:);
CVTcombS = CVTcomb(logical(CVTcomb.issyph),:);
CVTcombCt = CVTcomb(logical(CVTcomb.isCt),:);
CVTcombNg = CVTcomb(logical(CVTcomb.isNg),:);

HIVflow = table( ["S";"I";"C";"P";"Ip";"Cp";"T"], ["S";"S";"S";"P";"P";"P";"T"], 'VariableNames',{'outState','inState'});

compSSSS = tabComp(sum(table2array(tabComp(:,1:n+m+p))==repmat("S",1,n+m+p),2)==n+m+p,:);
compPSSS = tabComp(sum(table2array(tabComp(:,1:n+m+p))==["P",repmat("S",1,m+p)],2)==n+m+p,:);
compTSSS = tabComp(sum(table2array(tabComp(:,1:n+m+p))==["T",repmat("S",1,m+p)],2)==n+m+p,:);

for k=1:nbCompartments
    currentComp = tabComp(k,:);
    if sum(currentComp.no==[compSSSS.no,compPSSS.no,compTSSS.no])==0 %if not fully susceptible
        if currentComp.Ct ~= "IS" && currentComp.Ng ~= "IS"
            if currentComp.Ng == "E" || currentComp.Ng == "IA"
                %M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-rho_hg-rho_sg-rho_cg-rho_hsg-rho_hcg-rho_hscg;
                for kitNo=1:size(CVTcombNg.kit,1)
                    %Departure compartment is the current one
                    M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-CVTcombNg.rho(kitNo);
                    
                    % Finding the reception compartment
                    currentKit = CVTcombNg.kit(kitNo);
                    otherFixedStates = currentComp(:,setdiff(dis,currentKit{:}));
                    receptComp = otherFixedStates;
                    for inf=currentKit{:}
                        if inf ~= "HIV"
                            receptComp.(inf) = "S";
                        else
                            receptComp.HIV = HIVflow(HIVflow.outState==currentComp.HIV,:).inState;
                        end
                    end
                    receiverComp = innerjoin(tabComp,receptComp);
                    M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+CVTcombNg.rho(kitNo);
                    %disp([join([table2array(currentComp(:,1:n+m+p)),"-",string(CVTcombNg.rho(kitNo)),"->",table2array(receiverComp(:,1:n+m+p))])])
                end
            else
                if currentComp.Ct == "E" || currentComp.Ct == "IA"
                    %copy paste with some adaptation from the block before
                    for kitNo=1:size(CVTcombCt.kit,1)
                        %Departure compartment is the current one
                        M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-CVTcombCt.rho(kitNo);
                        
                        % Finding the reception compartment
                        currentKit = CVTcombCt.kit(kitNo);
                        otherFixedStates = currentComp(:,setdiff(dis,currentKit{:}));
                        receptComp = otherFixedStates;
                        for inf=currentKit{:}
                            if inf ~= "HIV"
                                receptComp.(inf) = "S";
                            else
                                receptComp.HIV = HIVflow(HIVflow.outState==currentComp.HIV,:).inState;
                            end
                        end
                        receiverComp = innerjoin(tabComp,receptComp);
                        M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+CVTcombCt.rho(kitNo);
                        %disp([join([table2array(currentComp(:,1:n+m+p)),"-",string(CVTcombCt.rho(kitNo)),"->",table2array(receiverComp(:,1:n+m+p))])])
                    end
                else
                    if currentComp.syph == "E" || contains(currentComp.syph,"I")
                        for kitNo=1:size(CVTcombS.kit,1)
                            %Departure compartment is the current one
                            M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-CVTcombS.rho(kitNo);
                            
                            % Finding the reception compartment
                            currentKit = CVTcombS.kit(kitNo);
                            otherFixedStates = currentComp(:,setdiff(dis,currentKit{:}));
                            receptComp = otherFixedStates;
                            for inf=currentKit{:}
                                if inf ~= "HIV"
                                    receptComp.(inf) = "S";
                                else
                                    receptComp.HIV = HIVflow(HIVflow.outState==currentComp.HIV,:).inState;
                                end
                            end
                            receiverComp = innerjoin(tabComp,receptComp);
                            M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+CVTcombS.rho(kitNo);
                            %disp([join([table2array(currentComp(:,1:n+m+p)),"-",string(CVTcombS.rho(kitNo)),"->",table2array(receiverComp(:,1:n+m+p))])])
                        end
                    else
                        %in this case, the only boxes left are
                        %isss,csss,ipsss,cpsss
                        disp(join(['hiv infected only:',table2array(currentComp(:,1:n+m+p))]))
                        for kitNo=1:size(CVTcombHIV.kit,1)
                            %Departure compartment is the current one
                            M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-CVTcombHIV.rho(kitNo);
                            
                            % Finding the reception compartment
                            currentKit = CVTcombHIV.kit(kitNo);
                            otherFixedStates = currentComp(:,setdiff(dis,currentKit{:}));
                            receptComp = otherFixedStates;
                            for inf=currentKit{:}
                                receptComp.(inf) = "S";
                            end
                            receiverComp = innerjoin(tabComp,receptComp);
                            M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+CVTcombHIV.rho(kitNo);
                            %disp([join([table2array(currentComp(:,1:n+m+p)),"-",string(CVTcombHIV.rho(kitNo)),"->",table2array(receiverComp(:,1:n+m+p))])])
                        end
                    end
                end
            end
        else
            %no voluntary testing
            disp(join(['no VT:',table2array(currentComp(:,1:n+m+p))]))
        end
    end
end

%% Adding mandatory voluntary testing for STIs under PrEP

for k=1:nbCompartments
    currentComp = tabComp(k,:);
    if sum(currentComp.no==[compSSSS.no,compPSSS.no,compTSSS.no])==0
        if 1%currentComp.Ct ~= "IS" && currentComp.Ng ~= "IS"
            if currentComp.HIV=="P" || currentComp.HIV=="Ip" || currentComp.HIV=="Cp"
                M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-eta_c;
                recepComp = currentComp(:,1:4); recepComp.Ct ="S";
                receiverComp = innerjoin(tabComp,recepComp);
                M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+eta_c;
                
                M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-eta_g;
                recepComp = currentComp(:,1:4); recepComp.Ng ="S";
                receiverComp = innerjoin(tabComp,recepComp);
                M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+eta_g;
                
                M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-eta_s;
                recepComp = currentComp(:,1:4); recepComp.syph ="S";
                receiverComp = innerjoin(tabComp,recepComp);
                M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+eta_s;
                
            elseif currentComp.HIV=="T"
                M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-eta_c_art;
                recepComp = currentComp(:,1:4); recepComp.Ct ="S";
                receiverComp = innerjoin(tabComp,recepComp);
                M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+eta_c_art;
                
                M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-eta_g_art;
                recepComp = currentComp(:,1:4); recepComp.Ng ="S";
                receiverComp = innerjoin(tabComp,recepComp);
                M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+eta_g_art;
                
                M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-eta_s_art;
                recepComp = currentComp(:,1:4); recepComp.syph ="S";
                receiverComp = innerjoin(tabComp,recepComp);
                M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+eta_s_art;
            end
        end
    end
end


%% Converting the matrix product to ODE system
[C,dC,eqn] = matToODE_v2(nbCompartments,M,B);
eqn.'
tabComp.X = C;

%% Write the ODE system in a text file
pathW = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Mymodels\';
fileID = fopen([pathW,'ODE_SICTPSEIIISSEIIS_v3.txt'],'w');
for k=1:size(tabComp,1)
    fprintf(fileID,'%12s%1s\r\n ',eqn(k),';');    
end
fclose(fileID);


%% Test du systeme d'ODE (avec ode45)
clearvars -except tabComp M; close all;
b=2;
[paramTab,mu,vecAlphas] = sampleParameters_v3(true,true,true,true,b);   %Ct,Ng,HIV,syph

%paramTab{3}.modelType='SICTP';
paramTab{3}.p = 0.6; paramTab{3}.eta = 4; paramTab{3}.zeta=randPERT(46,60,71,1)/100;
%paramTab{3}.alpha0 = paramTab{3}.alpha;
paramTab{3}.mu = mu;
[paramTab{3}.RSICTP,~,paramTab{3}.alpha] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
    paramTab{3}.theta,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
    paramTab{3}.p,paramTab{3}.mu,b,paramTab{3}.rhob);
disp(paramTab{3}.RSICTP)
%vecAlphas(3) = paramTab{3}.alpha;

%%
Y0 = ones(560,1); tspan=[1,100];
[res] = ode45( @(t,Y) ODE_SICTPSEIIISSEIIS2_v2(t,Y,paramTab{3}.betaI,paramTab{3}.betaC,paramTab{3}.sigma,paramTab{3}.theta,paramTab{3}.zeta,paramTab{3}.eta,paramTab{3}.p,...
    paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.gamma3,paramTab{4}.tau,paramTab{4}.theta,...
    paramTab{1}.beta,paramTab{1}.gamma,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,...
    paramTab{2}.beta,paramTab{2}.gamma,paramTab{2}.nu,paramTab{2}.eps,paramTab{2}.sigma,...
    NaN,mu,b,...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),...
    tspan, Y0);

sum(res.y(:,end))
b/mu

% Comparaison avec le modèle SICTP
tic
Y0 = ones(560,1)/560; tspan=1:10;
tspan=[1,100];
[res4dis] = ode45( @(t,Y) ODE_SICTPSEIIISSEIIS2_v2(t,Y,paramTab{3}.betaI,paramTab{3}.betaC,paramTab{3}.sigma,paramTab{3}.theta,paramTab{3}.zeta,paramTab{3}.eta,paramTab{3}.p,...
    paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.gamma3,paramTab{4}.tau,paramTab{4}.theta,...
    0,paramTab{1}.gamma,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,...
    0,paramTab{2}.gamma,paramTab{2}.nu,paramTab{2}.eps,paramTab{2}.sigma,...
    NaN,mu,b,...
    0,0,0,0,...
    0,0,0, 0,0, 0.1,0,0,0,0,0,0,0,0,0,0,0),...
    tspan, Y0);

toc
%sum(res4dis.y(:,end))
%b/mu
%only HIV:
[sum(res4dis.y(tabComp(tabComp.HIV=="S",:).no,end));...
    sum(res4dis.y(tabComp(tabComp.HIV=="I",:).no,end));...
    sum(res4dis.y(tabComp(tabComp.HIV=="C",:).no,end));...
    sum(res4dis.y(tabComp(tabComp.HIV=="P",:).no,end));...
    sum(res4dis.y(tabComp(tabComp.HIV=="Ip",:).no,end));...
    sum(res4dis.y(tabComp(tabComp.HIV=="Cp",:).no,end));...
    sum(res4dis.y(tabComp(tabComp.HIV=="T",:).no,end))]'

%only syph:
[sum(res4dis.y(tabComp(tabComp.syph=="S",:).no,end));...
    sum(res4dis.y(tabComp(tabComp.syph=="E",:).no,end));...
    sum(res4dis.y(tabComp(tabComp.syph=="I1",:).no,end));...
    sum(res4dis.y(tabComp(tabComp.syph=="I2",:).no,end));...
    sum(res4dis.y(tabComp(tabComp.syph=="I3",:).no,end))]'

%% Test du systeme d'ODE (avec solve)
clearvars -except tabComp M; close all;
b=2;
[paramTab,mu,vecAlphas] = sampleParameters_v3(true,true,true,true,b);   %Ct,Ng,HIV,syph
%
%paramTab{3}.modelType='SICTP';
paramTab{3}.p = 0.6; paramTab{3}.eta = 4; paramTab{3}.zeta=randPERT(46,60,71,1)/100;
%paramTab{3}.alpha0 = paramTab{3}.alpha;
paramTab{3}.mu = mu;
[paramTab{3}.RSICTP,~,paramTab{3}.alpha] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
    paramTab{3}.theta,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
    paramTab{3}.p,paramTab{3}.mu,b,paramTab{3}.rhob);
disp(paramTab{3}.RSICTP)

%Parameters
betaIh = paramTab{3}.betaI;
betaCh = paramTab{3}.betaC;
sigmah = paramTab{3}.sigma;
thetah = paramTab{3}.theta;
etah = paramTab{3}.eta;
zetah = paramTab{3}.zeta;
ph = paramTab{3}.p;

betaS = paramTab{4}.beta;
sigmaS = paramTab{4}.sigma;
gamma1S = 0;
gamma3S = paramTab{4}.gamma3;
nuS = paramTab{4}.nu;
tauS = paramTab{4}.tau;
thetaS = paramTab{4}.theta;
betaC = paramTab{1}.beta;
gammaC = paramTab{1}.gamma;
nuC = paramTab{1}.nu;
epsC = paramTab{1}.eps;
sigmaC = paramTab{1}.sigma;

betaG = paramTab{2}.beta;
gammaG = paramTab{2}.gamma;
nuG = paramTab{2}.nu;
epsG = paramTab{2}.eps;
sigmaG = paramTab{2}.sigma;

rho_h = 0;%0.1;
rho_s = 0;%paramTab{4}.rhob;
rho_c = 0;%paramTab{1}.rhob;
rho_g = 0;%paramTab{2}.rhob;

rho_hs = 0;%0.1;
rho_hc = 0;%.2;
rho_hg = 0;%0.05;
rho_sc = 0;%0.05;
rho_sg = 0;
rho_cg = 0.;
rho_hsc = 0.0;
rho_hsg = 0.1;
rho_hcg = 0;
rho_scg = 0;
rho_hscg = 0;
eta_s = 1;
eta_c = 0%4;
eta_g = 0%4;
eta_s_art = 1;%0;%1;
eta_c_art = 1;%0;%4;
eta_g_art = 1;%0;%4;


%a priori on the populations dynamics, that should help the algorithm to
%converge faster
% Y0 = findApriori(betaIh,betaCh,sigmah,thetah,zetah,etah,ph,...
%                 betaS,sigmaS,gamma1S,gamma3S,tauS,thetaS,nuS,...
%                 betaC,gammaC,nuC,epsC,sigmaC,...
%                 betaG,gammaG,nuG,epsG,sigmaG,...
%                 mu,b,...
%                 rho_h, rho_s,rho_c,rho_g);

restart = true; iterNo=0;
Y0 = ones(560,1)/560;

solveur="fsolve"; tspan=[0,50];
while restart && iterNo<100
    iterNo=iterNo+1
    if iterNo>5
        solveur="ode45";
    end
    if solveur=="fsolve"
        tic
        options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-5,...
            'Algorithm','trust-region','SubproblemAlgorithm','cg');
        [ES,autre] = fsolve(@(Y) ODE_SICTPSEIIISSEIIS2_v2(0,Y,...
            betaIh,betaCh,sigmah,thetah,zetah,etah,ph,...
            betaS,sigmaS,gamma3S,tauS,thetaS,...
            betaC,gammaC,nuC,epsC,sigmaC,...
            betaG,gammaG,nuG,epsG,sigmaG,...
            0,mu,b,...
            rho_h,rho_s,rho_c,rho_g,...
            rho_hs,rho_hc,rho_hg,rho_sc,rho_sg,rho_cg,...
            rho_hsc,rho_hsg, rho_hcg, rho_scg,...
            rho_hscg,...
            eta_s,eta_c,eta_g,...
            eta_s_art,eta_c_art,eta_g_art),...
            Y0,options);
        Y0 = rand(561,1); %pour le coup d'apres
        toc
    elseif solveur=="ode45"
        disp('passe dans la boucle ode45')
        tspan=[tspan(end),tspan(end)+50];
        [res] = ode45(@(t,Y) ODE_SICTPSEIIISSEIIS2_v2(t,Y,...
            betaIh,betaCh,sigmah,thetah,zetah,etah,ph,...
            betaS,sigmaS,gamma3S,tauS,thetaS,...
            betaC,gammaC,nuC,epsC,sigmaC,...
            betaG,gammaG,nuG,epsG,sigmaG,...
            0,mu,b,...
            rho_h,rho_s,rho_c,rho_g,...
            rho_hs,rho_hc,rho_hg,rho_sc,rho_sg,rho_cg,...
            rho_hsc,rho_hsg, rho_hcg, rho_scg,...
            rho_hscg,...
            eta_s,eta_c,eta_g,...
            eta_s_art,eta_c_art,eta_g_art),...
            tspan, Y0);
        ES = res.y(:,end);
        Y0 = ES;
    end
    
    popTot = sum(ES(1:560))
    b/mu
    
    %only HIV:
    % res_4_h=[sum(ES(tabComp(tabComp.HIV=="S",:).no,end));...
    % sum(ES(tabComp(tabComp.HIV=="I",:).no,end));...
    % sum(ES(tabComp(tabComp.HIV=="C",:).no,end));...
    % sum(ES(tabComp(tabComp.HIV=="P",:).no,end));...
    % sum(ES(tabComp(tabComp.HIV=="Ip",:).no,end));...
    % sum(ES(tabComp(tabComp.HIV=="Cp",:).no,end));...
    % sum(ES(tabComp(tabComp.HIV=="T",:).no,end))];
    res_4_h=[sum(ES(1:7:554)), sum(ES(2:7:555)), sum(ES(3:7:556)), sum(ES(4:7:557)),...
        sum(ES(5:7:558)), sum(ES(6:7:559)), sum(ES(7:7:560))];
    disp(['4dis,HIV: ' ,num2str(res_4_h)])
    
    %sictp model
    X0 = ones(7,1);
    [res_sictp] = ode45(@(t,Y) ODE_SICTPrEP(t,Y,b,betaIh,betaCh,sigmah,...
        thetah,ph,zetah,etah,mu,rho_h,'frequency'),...
        [1,500], X0);
    disp(['SICTP: ' ,num2str(res_sictp.y(:,end)')])
    
    
    %only syph:
    % res_4_s = [sum(ES(tabComp(tabComp.syph=="S",:).no,end));...
    % sum(ES(tabComp(tabComp.syph=="E",:).no,end));...
    % sum(ES(tabComp(tabComp.syph=="I1",:).no,end));...
    % sum(ES(tabComp(tabComp.syph=="I2",:).no,end));...
    % sum(ES(tabComp(tabComp.syph=="I3",:).no,end))];
    res_4_s = [sum(ES(reshape(repmat((1:35:554),7,1)+[0:6]',[1,560/5]))),...
        sum(ES(reshape(repmat((8:35:561),7,1)+[0:6]',[1,560/5]))),...
        sum(ES(reshape(repmat((15:35:568),7,1)+[0:6]',[1,560/5]))),...
        sum(ES(reshape(repmat((22:35:569),7,1)+[0:6]',[1,560/5]))),...
        sum(ES(reshape(repmat((29:35:576),7,1)+[0:6]',[1,560/5])))];
    disp(['4dis,s: ' ,num2str(res_4_s)])
    
    %seiiis model
    X0 = ones(5,1);
    [res_seiiis] = ode45(@(t,Y) ODE_SEIIIS_v4(t,Y,betaS,sigmaS,tauS,...
        thetaS,gamma1S,gamma3S,nuS,rho_s,mu,b),...
        [1 200], X0);
    disp(['SEIIIS: ' ,num2str(res_seiiis.y(:,end)')])
    [Rp,~,alphaS] = Rp_SEIIIS_v4(betaS,sigmaS,tauS,...
        nuS,gamma1S,thetaS,gamma3S,mu,b,0);
    
    
    % seiis model Ct
    %only Ct:
    % res_4_Ct = [sum(ES(tabComp(tabComp.Ct=="S",:).no,end));...
    % sum(ES(tabComp(tabComp.Ct=="E",:).no,end));...
    % sum(ES(tabComp(tabComp.Ct=="IA",:).no,end));...
    % sum(ES(tabComp(tabComp.Ct=="IS",:).no,end))];
    
    res_4_Ct = [sum(ES(reshape(repmat((1:140:421),35,1)+[0:34]',[1,560/4]))),...
        sum(ES(reshape(repmat((36:140:456),35,1)+[0:34]',[1,560/4]))),...
        sum(ES(reshape(repmat((71:140:491),35,1)+[0:34]',[1,560/4]))),...
        sum(ES(reshape(repmat((106:140:527),35,1)+[0:34]',[1,560/4])))];
    disp(['4dis,Ct: ' ,num2str(res_4_Ct)])
    
    X0 = ones(4,1);
    [res_seiis] = ode45(@(t,Y) ODE_SEIIS_v4(t,Y,betaC,nuC,gammaC,...
        sigmaC,epsC,rho_c,mu,b),...
        [1 200], X0);
    disp(['SEIIS,Ct: ' ,num2str(res_seiis.y(:,end)')])
    
    % seiis model Ng
    % res_4_Ng = [sum(ES(tabComp(tabComp.Ng=="S",:).no,end));...
    % sum(ES(tabComp(tabComp.Ng=="E",:).no,end));...
    % sum(ES(tabComp(tabComp.Ng=="IA",:).no,end));...
    % sum(ES(tabComp(tabComp.Ng=="IS",:).no,end))];
    res_4_Ng = [sum(ES(1:140)),...
        sum(ES(141:280)),...
        sum(ES(281:420)),...
        sum(ES(421:560))];
    disp(['4dis,Ng: ' ,num2str(res_4_Ng)])
    
    X0 = ones(4,1);
    [res_seiis] = ode45(@(t,Y) ODE_SEIIS_v4(t,Y,betaG,nuG,gammaG,...
        sigmaG,epsG,rho_g,mu,b),...
        [1 200], X0);
    disp(['SEIIS,Ng: ' ,num2str(res_seiis.y(:,end)')])
    
    N=b/mu; myTol=1e-3;
    %alors que les infecitons ne doivent pas disparaitre:
    
    %abs(res_4_h(1)+res_4_h(4)-N)<myTol || abs(res_4_s(1)-N)<myTol ||...
    %abs(res_4_Ct(1)-N)<myTol || abs(res_4_Ng(1)-N)<myTol ||...
    if (abs(popTot-N)>myTol || sum(ES<-myTol)~=0)
        %nothing
    else
        restart=false;
        
    end
end
















if (0)
    %si : not SSSS et TSSS et not PSSS
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
    
    
    
    %% Verifications of the calculations
    clear all; close all;
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
    
    PrevC = sum(popEnd(table2array(tabComp(tabComp.Ct=="IA" | tabComp.Ct=="IS" | tabComp.Ct=="E" ,"no"))))./N
    1-1./RC
    
    PrevN = sum(popEnd(table2array(tabComp(tabComp.Ng=="IA" | tabComp.Ng=="IS" | tabComp.Ng=="E" ,"no"))))./N
    1-1./RN
    
    
    
    %%
    %
    % tspan = 0:.1:1000;
    % Y0 = (b/mu)*ones(80,1)';
    % options = odeset('AbsTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
    % % odefun = @(t,Y) ODE_SICRSEIIISSEIIS_3(t,Y,betaIHIV,betaCHIV,thetah,gammah,sigmah,...
    % %                                       betaSyph,sigmas,taus,thetas,gamma1s,gamma3s,nus,...
    % %                                       beta1,nuCt,epsCt,sigmaCt,gammaCt,mu,b,rho,indexInfSyph,indexInfSEIIS,indexIHIV,indexCHIV);
    % odefun = @(t,Y) ODE_SEIISSICRSEIIIS_3(t,Y,beta1,gammaCt,nuCt,sigmaCt,epsCt,...
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
    % R1 = (beta1*sigmaCt*(gammaCt*(1-epsCt) + mu + nuCt + epsCt*rho))./((mu + sigmaCt+rho).*(gammaCt + mu + nuCt).*(mu + nuCt + rho)); %ok
    % S1  = b/(mu*max(R1,1));
    % lambda1 = max(beta1*(R1-1)*(sigmaCt+rho+mu)/(beta1+(sigmaCt+rho+mu)*R1),0);
    % E1  = lambda1*S1/(sigmaCt+mu+rho);
    % I1  = (1-epsCt)*sigmaCt*E1/(nuCt+rho+mu);
    % J1  = epsCt*sigmaCt*E1/(gammaCt+nuCt+mu);
    % [S1,E1,I1,J1]
    % [ts,Ys1] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta1,nuCt,epsCt,sigmaCt,gammaCt,rho,mu),tspan,[1,1,1,1], options);
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
    % [beta1,beta2,gammaCt,gammaNg,nuCt,nuNg,sigmaCt,sigmaNg,rho] = random_parameters(true, true);
    % epsCt=0.2;epsNg=0.3;
    % RHIV = (betaIHIV*(thetah+mu) + betaCHIV*sigmah)/((mu+thetah)*(sigmah+gammah+mu));
    % R1 = (beta1*sigmaCt*(gammaCt*(1-epsCt) + mu + nuCt + epsCt*rho))./((mu + sigmaCt+rho).*(gammaCt + mu + nuCt).*(mu + nuCt + rho)); %ok
    % R2 = (beta2*sigmaNg*(gammaNg*(1-epsNg) + mu + nuNg + epsNg*rho))./((mu + sigmaNg+rho).*(gammaNg + mu + nuNg).*(mu + nuNg + rho)); %ok
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
    % %                                       beta1,nuCt,epsCt,sigmaCt,gammaCt,...
    % %                                       beta2,nuNg,epsNg,sigmaNg,gammaNg,mu,b,rho,...
    % %                                       indexInfSyph,indexInf1,indexInf2,indexIHIV,indexCHIV);
    %
    % odefun = @(t,Y) ODE_SEIIS2SICRSEIIIS_3(t,Y,beta1,gammaCt,nuCt,sigmaCt,epsCt,...
    %                                          beta2,gammaNg,nuNg,sigmaNg,epsNg,...
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
    % lambda1 = max(beta1*(R1-1)*(sigmaCt+rho+mu)/(beta1+(sigmaCt+rho+mu)*R1),0);
    % E1  = lambda1*S1/(sigmaCt+mu+rho);
    % I1  = (1-epsCt)*sigmaCt*E1/(nuCt+rho+mu);
    % J1  = epsCt*sigmaCt*E1/(gammaCt+nuCt+mu);
    % [S1,E1,I1,J1]
    % [ts,Ys1] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta1,nuCt,epsCt,sigmaCt,gammaCt,rho,mu),tspan,[1,1,1,1], options);
    % Ys1(end,:)
    %
    % %STI 2
    % indexS2 = find(compartments(:,4)==1); indexE2 = find(compartments(:,4)==2); indexI2 = find(compartments(:,4)==3); indexJ2 = find(compartments(:,4)==4);
    % [sum(T(indexS2)),sum(T(indexE2)),sum(T(indexI2)),sum(T(indexJ2))]
    % S2  = b/(mu*max(R2,1));
    % lambda2 = max(beta2*(R2-1)*(sigmaNg+rho+mu)/(beta2+(sigmaNg+rho+mu)*R2),0);
    % E2  = lambda2*S2/(sigmaNg+mu+rho);
    % I2  = (1-epsNg)*sigmaNg*E2/(nuNg+rho+mu);
    % J2  = epsNg*sigmaNg*E2/(gammaNg+nuNg+mu);
    % [S2,E2,I2,J2]
    % [ts,Ys2] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta2,nuNg,epsNg,sigmaNg,gammaNg,rho,mu),tspan,[1,1,1,1], options);
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
end







%% Test du systeme d'ODE (avec solve)
%tabRecap=[];
%clearvars -except tabRecap; 
clear all; close all;
tabRecap=[];
b=2; N=100;
P_inf_tot = zeros(N,4);
vecpHIV=[0.3,0.4,0.5,0.6,0.7];
k=1;
for pHIV = vecpHIV
for i=1:N
    disp([i,pHIV])
    tic
    [paramTab,mu,vecAlphas] = sampleParameters_v3(true,true,true,true,b);   %Ct,Ng,HIV,syph
    %
    paramTab{3}.modelType='SICTP';
    paramTab{3}.p = pHIV; paramTab{3}.eta = 4; paramTab{3}.zeta=randPERT(46,60,71,1)/100;
    paramTab{3}.alpha0 = paramTab{3}.alpha;
    paramTab{3}.mu = mu;
    [paramTab{3}.RSICTP,~,paramTab{3}.alpha] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
        paramTab{3}.theta,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
        paramTab{3}.p,paramTab{3}.mu,b,paramTab{3}.rhob);
    disp(paramTab{3}.RSICTP)
    
    paramRho.rho_h = paramTab{3}.rhob;
    paramRho.rho_s = paramTab{4}.rhob;
    paramRho.rho_c = paramTab{1}.rhob;
    paramRho.rho_g = paramTab{2}.rhob;
    
    paramRho.rho_hs = 0;
    paramRho.rho_hc = 0;
    paramRho.rho_hg = 0;
    paramRho.rho_sc = 0;
    paramRho.rho_sg = 0;
    paramRho.rho_cg = 0;
    paramRho.rho_hsc = 0;
    paramRho.rho_hsg = 0;
    paramRho.rho_hcg = 0;
    paramRho.rho_scg = 0;
    paramRho.rho_hscg = 0;
    paramRho.eta_s_prep = 1;
    paramRho.eta_c_prep = 4;
    paramRho.eta_g_prep = 4;
    paramRho.eta_s_art = 0;
    paramRho.eta_c_art = 0;
    paramRho.eta_g_art = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    
    c=0;f=1;
    %paramTab{2}.beta= 13.035779548202
    [U12,dU,P,P_inf] = U1234_SICTPSEIIISSEIIS2_v2(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,c,f);
    
    P_inf_tot(i,1)=P_inf(1);
    P_inf_tot(i,2)=P_inf(2);
    P_inf_tot(i,3)=P_inf(3);
    P_inf_tot(i,4)=P_inf(4);
    toc    
end

tabRecap(k,1) = pHIV;
tabRecap(k,2) = paramRho.eta_s_art;
tabRecap(k,3) = paramRho.eta_c_art;
tabRecap(k,4) = paramRho.eta_g_art;
tabRecap(k,5) = sum(P_inf_tot(:,1)>0.005)/N;
tabRecap(k,6) = sum(P_inf_tot(:,2)>0.005)/N;
tabRecap(k,7) = sum(P_inf_tot(:,3)>0.005)/N;
tabRecap(k,8) = sum(P_inf_tot(:,4)>0.005)/N;
k=k+1;
end
load handel
sound(y,Fs)
%pour resoudre le probleme d'optim selon une variable rho particuliere :...
%faire une fonction auxiliaire qui prend en entree cette variable et qui
%appelle le system d'ODE
tabRecap





%%
clear all;close all;
vecPCt=[]; vecPNg=[]; vecPHIV=[]; vecPS=[]; vecPprevagay=[];
b=2; 
for i=1:500
    [paramTab,mu,vecAlphas] = sampleParameters_v3(true,true,true,true,b);
    vecPCt(i)= paramTab{1}.P;
    vecPNg(i)= paramTab{2}.P;
    vecPHIV(i)= paramTab{3}.PSICT;
    vecPS(i)= paramTab{4}.P;
    vecPprevagay(i)= paramTab{3}.Pprev;
end


close all;
fig = figure()
BinWidth=0.01%(max(alphaTot(:))-min(alphaTot(:)))/60;
histogram(vecPHIV,'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(vecPS,'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(vecPCt,'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(vecPNg,'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
legend('$\Pi_h(\rho_0^h))$ (prevenir)',...
    '$\Pi_s(\rho_0^s)$ (adapt. prevenir)',...
    '$\Pi_c(\rho_0^c)$ (adapt. prevenir)',...
    '$\Pi_g(\rho_0^g)$ (adapt. prevenir)',...
    'Interpreter','latex','FontSize',12,'Box','off')


mean(vecPHIV)/mean(vecPprevagay)


