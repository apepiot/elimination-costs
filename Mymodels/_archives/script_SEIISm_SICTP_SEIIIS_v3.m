%% This code generates the matrix and the ODE system associated to a SICR^n x SEIIIS^p x SEIIS^m model
% combined and targeted testing are mixed

%there is a problem in the algo of assigning voluntary testing rates, see
%the following version (and notes in the blue notebook  pages 129-130)
clear all;

%% Parameter initialization
% Nothing to change below, this code generates a 4-disease model.
n = 1; %number of SICTPrEP (min=?, max=1)
p = 1; %number of SEIIIS (min=?,max=1)
m = 2; %number of SEIIS (min=?, max m=4)

syms mu b;
syms LambdaCt; syms nuCt; syms epsCt; syms sigmaCt; syms gammaCt;
syms LambdaNg; syms nuNg; syms epsNg; syms sigmaNg; syms gammaNg;

syms Lambdah thetah sigmah ph eta_h_prep zetah
syms Lambdas sigmas taus thetas gamma1s gamma3s nus
% targeted testing rate:
syms rho_h rho_s rho_c rho_g
syms eta_c_prep eta_g_prep eta_s_prep
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

boxesSICTP  = ["S","I","C","P","Ip","Cp","T"]; %should be the same order than in the ODE_SICTPrEP.m
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
        idx  = strfind(otherDis(:,:),'Ct','ForceCellOutput',true);   mSTI=max([find([idx{:}]),0]);
        idx  = strfind(otherDis(:,:),'Ng','ForceCellOutput',true);   mSTI=mSTI+max([find([idx{:}]),0]);
        idx  = strfind(otherDis(:,:),'HIV','ForceCellOutput',true);  nHIV=max([find([idx{:}]),0]);
        idx  = strfind(otherDis(:,:),'syph','ForceCellOutput',true); pSYPH=max([find([idx{:}]),0]);
        otherDiseaseStatesConstant = createTableComp(mSTI,nHIV,pSYPH,boxesSEIIS,boxesSICTP,boxesSEIIIS,otherDis);
        for i=1:size(otherDiseaseStatesConstant,1)
            %for every combination of disease states (other than HIV), applies the matrix
            %of the ODE system of HIV
            otherStates = otherDiseaseStatesConstant(i,:);
            T = innerjoin(otherStates,tabComp); %sub-table with otherStates and all the states of HIV
            
            %makes sure that the order of the states of inf is good (i.e. same that in the ODE system)
            [found, idx] = ismember(table2array(T(:,inf)), boxesInf');
            [~, sortorder] = sort(idx); newT = T(sortorder,:);
            
            %matrix of the ODE system for HIV
            %we don't include: rho, mu, b.
            if contains(inf,'HIV')
                Mh = M_SICTP(Lambdah,thetah,sigmah,ph,zetah,eta_h_prep,0,0,0);
            elseif contains(inf,'syph')
                Mh = M_SEIIIS(Lambdas,sigmas,taus,thetas,gamma1s,gamma3s,nus,0,0,0);
            elseif contains(inf,'Ct')
                Mh = M_SEIIS(LambdaCt,epsCt,nuCt,gammaCt,sigmaCt,0,0,0);
            elseif contains(inf,'Ng')
                Mh = M_SEIIS(LambdaNg,epsNg,nuNg,gammaNg,sigmaNg,0,0,0);
            end
            M(newT.no,newT.no) = M(newT.no,newT.no) + Mh;
        end
    end
end

%% Adds additionnal rate: pi (new individuals in the model)
syms B [nbCompartments 1]; %contains other rate, not variable dependent (e.g. pi)
B(:) = 0;

%Adds input parameters
if n==1
    B(sum(table2array(tabComp(:,1:(n+m+p)))==repmat("S",1,n+m+p),2)==n+m+p) = (1-ph)*b;
    B(sum(table2array(tabComp(:,1:(n+m+p)))==["P",repmat("S",1,m+p)],2)==n+m+p) = ph*b;
else
    error('need to add b')
end

%% To adapt the code for a k-disease model, k<n+p+m
% this adds columns for the other infections (the ones not considered in the model) and make them susceptible 
if m==1 && n==1 && p==0
    tabComp.Ng = repmat("S",nbCompartments,1);
    tabComp.syph = repmat("S",nbCompartments,1);
end

%% Adds combined voluntary testing rate rho_xxx
%Symptomatic individuals for STIs test because of symptoms and not with
%voluntary testing.
%Individuals on PrEP (P,Ip,Cp) test only routiely, not voluntarily.

% For all combination of infections *except if individuals are on PrEP*,
% we need to add combined voluntary testing rate
% e.g. for a given box *I(H)E(s)IA(C)E(G)*
% HIV and syph tested together, then combined testing leads to S(H)S(syph),
% Ct tested alone, targeted testing goes to S(C) (already accounted in the matrix M)
% Ng ---------------------------------------S(G) (idem)

%The following table shows 
CVTcomb = table({["HIV"];["syph"];["Ct"];["Ng"];["HIV","syph"];["HIV","Ct"];["HIV","Ng"];["syph","Ct"];["syph","Ng"];["Ct","Ng"];...
    ["HIV","syph","Ct"];["HIV","syph","Ng"];["HIV","Ct","Ng"];["syph","Ct","Ng"];...
    ["HIV","syph","Ct","Ng"]},...
    [1;0;0;0;1;1;1;0;0;0;1;1;1;0;1],...
    [0;1;0;0;1;0;0;1;1;0;1;1;0;1;1],...
    [0;0;1;0;0;1;0;1;0;1;1;0;1;1;1],...
    [0;0;0;1;0;0;1;0;1;1;0;1;1;1;1],...
    [rho_h;rho_s;rho_c;rho_g;rho_hs;rho_hc;rho_hg;rho_sc;rho_sg;rho_cg;rho_hsc;rho_hsg;rho_hcg;rho_scg;rho_hscg],...
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
    if sum(currentComp.no==[compSSSS.no,compPSSS.no,compTSSS.no])==0 %if not fully (non infected/infectious)
        if currentComp.HIV ~= "P" && currentComp.HIV ~= "Ip" && currentComp.HIV ~= "Cp" %if not under PrEP
            if currentComp.Ct ~= "IS" && currentComp.Ng ~= "IS" %if not symptomatics
                if currentComp.Ng == "E" || currentComp.Ng == "IA"
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
            else %symptmatics for Ct or Ng -> no voluntary testing at all, even for HIV and syphilis
                disp(join(['no VT:',table2array(currentComp(:,1:n+m+p))]))
            end
        else %if under PrEP
        end
    end
end


%% Adding mandatory/recommanded testing for STIs under PrEP
for k=1:nbCompartments
    currentComp = tabComp(k,:);
    if sum(currentComp.no==[compSSSS.no,compPSSS.no,compTSSS.no])==0
        if 1%currentComp.Ct ~= "IS" && currentComp.Ng ~= "IS"
            if currentComp.HIV=="P" || currentComp.HIV=="Ip" || currentComp.HIV=="Cp"
                M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-eta_c_prep;
                recepComp = currentComp(:,1:4); recepComp.Ct ="S";
                receiverComp = innerjoin(tabComp,recepComp);
                M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+eta_c_prep;
                
                M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-eta_g_prep;
                recepComp = currentComp(:,1:4); recepComp.Ng ="S";
                receiverComp = innerjoin(tabComp,recepComp);
                M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+eta_g_prep;
                
                M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no)-eta_s_prep;
                recepComp = currentComp(:,1:4); recepComp.syph ="S";
                receiverComp = innerjoin(tabComp,recepComp);
                M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+eta_s_prep;
                
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
[X,dX,eqn,F] = matToODE_v2(nbCompartments,M,B);
eqn.'
tabComp.X = X;

%% Write the ODE system in a text file
pathW = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Mymodels\';
fileID = fopen([pathW,'ODE_SICTPSEIIISSEIIS_v4.txt'],'w');
for k=1:size(tabComp,1)
    fprintf(fileID,'%12s%1s\r\n ',eqn(k),';');    
end
fclose(fileID);


%% Computation of the gradient
%syms dXkdXi [nbCompartments nbCompartments]

dXkdXi = []%zeros(nbCompartments,nbCompartments);

clear all;
pathW = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Mymodels\';
fileID = fopen([pathW,'gradF_ODE_SICTPSEIIISSEIIS_v4_2.txt'],'w');

F_ODESICTPSEIIISSEIIS2_v4_detailed
%syms dXkdXi [1 1]
for k=1:560
    for i=1:560
        %dXkdXi(k,i)=diff(F(k),X(i));
        res=diff(F(k),X(i));
        fprintf(fileID,'%10s%100s%1s\r\n ',['dF(',num2str(k),',',num2str(i),') = '],res,';'); 
    end
    disp([k,i]);
end
fclose(fileID);


%% Test du fichier du gradient:
tic
for i=1:10
grad = gradF_ODE_SICTPSEIIISSEIIS2_v4_fun(NaN,rand(560,1),paramTab{3}.betaI,paramTab{3}.betaC,paramTab{3}.sigma,paramTab{3}.theta,paramTab{3}.zeta,paramTab{3}.eta,paramTab{3}.p,...
    paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.gamma3,paramTab{4}.tau,paramTab{4}.theta,...
    0,paramTab{1}.gamma,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,...
    0,paramTab{2}.gamma,paramTab{2}.nu,paramTab{2}.eps,paramTab{2}.sigma,...
    mu,b,...
    0,0,0,0,...
    0,0,0,0,0,0.1,0,0,0,0,0,0,0,0,0,0,0);
disp(i);
toc
end
toc






















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
Y0 = ones(560,1); tspan=[1,500];
[res] = ode45( @(t,Y) ODE_SICTPSEIIISSEIIS2_v4(t,Y,paramTab{3}.betaI,paramTab{3}.betaC,paramTab{3}.sigma,paramTab{3}.theta,paramTab{3}.zeta,paramTab{3}.eta,paramTab{3}.p,...
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
Y0 = ones(560,1)/560; tspan=1:100;
tspan=[1,100];
[res4dis] = ode45(@(t,Y) ODE_SICTPSEIIISSEIIS2_v4(t,Y,paramTab{3}.betaI,paramTab{3}.betaC,paramTab{3}.sigma,paramTab{3}.theta,paramTab{3}.zeta,paramTab{3}.eta,paramTab{3}.p,...
    paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.gamma3,paramTab{4}.tau,paramTab{4}.theta,...
    0,paramTab{1}.gamma,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,...
    0,paramTab{2}.gamma,paramTab{2}.nu,paramTab{2}.eps,paramTab{2}.sigma,...
    NaN,mu,b,...
    0,0,0,0,...
    0,0,0,0,0,0.1,0,0,0,0,0,0,0,0,0,0,0),...
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
% clearvars -except tabComp M; close all;
% b=2;
% [paramTab,mu,vecAlphas] = sampleParameters_v3(true,true,true,true,b);   %Ct,Ng,HIV,syph
% %
% %paramTab{3}.modelType='SICTP';
% paramTab{3}.p = 0.6; paramTab{3}.eta = 4; paramTab{3}.zeta=randPERT(46,60,71,1)/100;
% %paramTab{3}.alpha0 = paramTab{3}.alpha;
% paramTab{3}.mu = mu;
% [paramTab{3}.RSICTP,~,paramTab{3}.alpha] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
%     paramTab{3}.theta,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
%     paramTab{3}.p,paramTab{3}.mu,b,paramTab{3}.rhob);
% disp(paramTab{3}.RSICTP)

%Parameters
betaIh = paramTab{3}.betaI;
betaCh = paramTab{3}.betaC;
sigmah = paramTab{3}.sigma;
thetah = paramTab{3}.theta;
eta_h_prep = paramTab{3}.eta;
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

rho_h = paramTab{3}.rhob;
rho_s = paramTab{4}.rhob;
rho_c = paramTab{1}.rhob;
rho_g = paramTab{2}.rhob;

rho_hs = 0;%0.1;
rho_hc = 0;%.2;
rho_hg = 0;%0.05;
rho_sc = 0;%0.05;
rho_sg = 0;
rho_cg = 0.;
rho_hsc = 0.0;
rho_hsg = 0.;
rho_hcg = 0;
rho_scg = 0;
rho_hscg = 0;
eta_s_prep = 0;
eta_c_prep = 0%4;
eta_g_prep = 0%4;
eta_s_art = 0;%0;%1;
eta_c_art = 0;%0;%4;
eta_g_art = 0;%0;%4;


%a priori on the populations dynamics, that should help the algorithm to
%converge faster
% Y0 = findApriori(betaIh,betaCh,sigmah,thetah,zetah,eta_h_prep,ph,...
%                 betaS,sigmaS,gamma1S,gamma3S,tauS,thetaS,nuS,...
%                 betaC,gammaC,nuC,epsC,sigmaC,...
%                 betaG,gammaG,nuG,epsG,sigmaG,...
%                 mu,b,...
%                 rho_h, rho_s,rho_c,rho_g);

% Finding the initial condition of the ODE system by solving the following system:

% One-disease models
[res_sictp] = ode45(@(t,Y) ODE_SICTPrEP(t,Y,b,betaIh,betaCh,sigmah,...
    thetah,ph,zetah,eta_h_prep,mu,rho_h,'frequency'),...
    [1,500], ones(7,1)); ES_sictp = res_sictp.y(:,end)';

[res_seiiis] = ode45(@(t,Y) ODE_SEIIIS_v4(t,Y,betaS,sigmaS,tauS,...
    thetaS,gamma1S,gamma3S,nuS,rho_s,mu,b),...
    [1 200],ones(5,1)); ES_seiiis = res_seiiis.y(:,end)';

[res_seiisCt] = ode45(@(t,Y) ODE_SEIIS_v4(t,Y,betaC,nuC,gammaC,...
    sigmaC,epsC,rho_c,mu,b),...
    [1 200], ones(4,1)); ES_ct = res_seiisCt.y(:,end)';

[res_seiisNg] = ode45(@(t,Y) ODE_SEIIS_v4(t,Y,betaG,nuG,gammaG,...
    sigmaG,epsG,rho_g,mu,b),...
    [1 200], ones(4,1)); ES_ng = res_seiisNg.y(:,end)';

Y0 = ones(561,1)/560;
options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-6,'MaxFunctionEvaluations',100);
% [ES,fval,exitflag,output] = fsolve(@(Y) CI_sol(Y, ES_sictp(1),ES_sictp(2),ES_sictp(3),ES_sictp(4),ES_sictp(5),ES_sictp(6),ES_sictp(7),...
%                                                ES_seiiis(1),ES_seiiis(2),ES_seiiis(3),ES_seiiis(4),ES_seiiis(5),...
%                                                ES_ct(1),ES_ct(2),ES_ct(3),ES_ct(4),...
%                                                ES_ng(1),ES_ng(2),ES_ng(3),ES_ng(4),b/mu),...
%                               Y0,options);
% 
% [A,B] = CI_sol_mat(ES_sictp(1),ES_sictp(2),ES_sictp(3),ES_sictp(4),ES_sictp(5),ES_sictp(6),ES_sictp(7),...
%                    ES_seiiis(1),ES_seiiis(2),ES_seiiis(3),ES_seiiis(4),ES_seiiis(5),...
%                    ES_ct(1),ES_ct(2),ES_ct(3),ES_ct(4),...
%                    ES_ng(1),ES_ng(2),ES_ng(3),ES_ng(4),...
%                    b/mu); 
 
%X = linsolve(A,B')
% XCI = lsqlin(A,B,[],[],[],[],zeros(560,1),b/mu*ones(560,1)) %lsqlin(C,d,A,b,Aeq,beq,lb,ub)
% 
% %verif de XCI
% F = CI_sol(XCI, ES_sictp(1),ES_sictp(2),ES_sictp(3),ES_sictp(4),ES_sictp(5),ES_sictp(6),ES_sictp(7),...
%                                                ES_seiiis(1),ES_seiiis(2),ES_seiiis(3),ES_seiiis(4),ES_seiiis(5),...
%                                                ES_ct(1),ES_ct(2),ES_ct(3),ES_ct(4),...
%                                                ES_ng(1),ES_ng(2),ES_ng(3),ES_ng(4),b/mu)
% Y0=XCI;

restart=true; iterNo=0;
solveur="ode45"; tspan=[0,50];
while restart && iterNo<100
    iterNo=iterNo+1
    if iterNo>10
        solveur="ode45";
    end
    if solveur=="fsolve"
        tic
        options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-6,'MaxFunctionEvaluations',100000,...
            'Algorithm','trust-region','SubproblemAlgorithm','cg');
        [ES,fval,exitflag,output] = fsolve(@(Y) ODE_SICTPSEIIISSEIIS2_v4(0,Y,...
                        betaIh,betaCh,sigmah,thetah,zetah,eta_h_prep,ph,...
                        betaS,sigmaS,gamma3S,tauS,thetaS,...
                        betaC,gammaC,nuC,epsC,sigmaC,...
                        betaG,gammaG,nuG,epsG,sigmaG,...
                        0,mu,b,...
                        rho_h,rho_s,rho_c,rho_g,...
                        rho_hs,rho_hc,rho_hg,rho_sc,rho_sg,rho_cg,...
                        rho_hsc,rho_hsg, rho_hcg, rho_scg,...
                        rho_hscg,...
                        eta_s_prep,eta_c_prep,eta_g_prep,...
                        eta_s_art,eta_c_art,eta_g_art),...
                        Y0,options);
              
        Y0 = rand(560,1)*b/mu; %pour le coup d'apres
        toc
    elseif solveur=="knitro"
        my_f = @(Y) ODE_SICTPSEIIISSEIIS2_v4(0,Y,...
                        betaIh,betaCh,sigmah,thetah,zetah,eta_h_prep,ph,...
                        betaS,sigmaS,gamma3S,tauS,thetaS,...
                        betaC,gammaC,nuC,epsC,sigmaC,...
                        betaG,gammaG,nuG,epsG,sigmaG,...
                        0,mu,b,...
                        rho_h,rho_s,rho_c,rho_g,...
                        rho_hs,rho_hc,rho_hg,rho_sc,rho_sg,rho_cg,...
                        rho_hsc,rho_hsg, rho_hcg, rho_scg,...
                        rho_hscg,...
                        eta_s_prep,eta_c_prep,eta_g_prep,...
                        eta_s_art,eta_c_art,eta_g_art);
        options = knitro_options('maxtime_real',900.0);
        %options = knitro_options('maxtime_real',10.0, 'ms_enable', 1, 'ms_maxsolves', 5, 'ms_maxtime_real', 500);
        [ES,fval,exitflag,output] = knitro_nlneqs(myf, Y0, {}, options);
        objfunc = @(Y) ([1]);
        confunc = @(Y) (disperse2([], myf(Y)));
        %objfunc = @(Y) ([norm(myf(Y))]);
        %confunc = @(Y) (disperse2([], []));
        lb=zeros(size(Y0));
        [ES,fval,exitflag,output,lambda,grad,hessian] = knitro_nlp(objfunc, Y0,...
                                                                   [], [], [], [],...
                                                                   lb, [], confunc, [], options);
        
        disp(['error=', num2str(max(abs(myf(ES))))])
        Y0 = rand(561,1)*b/mu; %pour le coup d'apres
    elseif solveur=="ode45"
        disp('passe dans la boucle ode45')
        tspan=[tspan(end),tspan(end)+50];
        [res] = ode45(@(t,Y) ODE_SICTPSEIIISSEIIS2_v4(t,Y,...
                            betaIh,betaCh,sigmah,thetah,zetah,eta_h_prep,ph,...
                            betaS,sigmaS,gamma3S,tauS,thetaS,...
                            betaC,gammaC,nuC,epsC,sigmaC,...
                            betaG,gammaG,nuG,epsG,sigmaG,...
                            0,mu,b,...
                            rho_h,rho_s,rho_c,rho_g,...
                            rho_hs,rho_hc,rho_hg,rho_sc,rho_sg,rho_cg,...
                            rho_hsc,rho_hsg, rho_hcg, rho_scg,...
                            rho_hscg,...
                            eta_s_prep,eta_c_prep,eta_g_prep,...
                            eta_s_art,eta_c_art,eta_g_art),...
                            tspan, Y0);
        ES = res.y(:,end);
        Y0 = ES;
    end
    
    popTot = sum(ES(1:560))
    b/mu
    
    %only HIV:
    res_4_h=[sum(ES(1:7:554)), sum(ES(2:7:555)), sum(ES(3:7:556)), sum(ES(4:7:557)),...
        sum(ES(5:7:558)), sum(ES(6:7:559)), sum(ES(7:7:560))];
    disp(['4dis,HIV: ' ,num2str(res_4_h)])  
    %from the sictp model
    disp(['SICTP: ' ,num2str(ES_sictp)])
    
    
    %only syph:
    res_4_s = [sum(ES(reshape(repmat((1:35:554),7,1)+[0:6]',[1,560/5]))),...
        sum(ES(reshape(repmat((8:35:561),7,1)+[0:6]',[1,560/5]))),...
        sum(ES(reshape(repmat((15:35:568),7,1)+[0:6]',[1,560/5]))),...
        sum(ES(reshape(repmat((22:35:569),7,1)+[0:6]',[1,560/5]))),...
        sum(ES(reshape(repmat((29:35:576),7,1)+[0:6]',[1,560/5])))];
    disp(['4dis,s: ' ,num2str(res_4_s)])
    
    %from the seiiis model
    disp(['SEIIIS: ' ,num2str(ES_seiiis)])
    [Rp,~,alphaS] = Rp_SEIIIS_v4(betaS,sigmaS,tauS,...
        nuS,gamma1S,thetaS,gamma3S,mu,b,0);
    
    
    % seiis model Ct
    %only Ct:
    res_4_Ct = [sum(ES(reshape(repmat((1:140:421),35,1)+[0:34]',[1,560/4]))),...
        sum(ES(reshape(repmat((36:140:456),35,1)+[0:34]',[1,560/4]))),...
        sum(ES(reshape(repmat((71:140:491),35,1)+[0:34]',[1,560/4]))),...
        sum(ES(reshape(repmat((106:140:527),35,1)+[0:34]',[1,560/4])))];
    disp(['4dis,Ct: ' ,num2str(res_4_Ct)]) 
    %from the Ct model
    disp(['SEIIS,Ct: ' ,num2str(ES_ct)])
    
    % seiis model Ng
    res_4_Ng = [sum(ES(1:140)),...
        sum(ES(141:280)),...
        sum(ES(281:420)),...
        sum(ES(421:560))];
    disp(['4dis,Ng: ' ,num2str(res_4_Ng)])    
    %from the Ng model
    disp(['SEIIS,Ng: ' ,num2str(ES_ng)])
    
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


