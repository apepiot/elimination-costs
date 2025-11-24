function [nbCompartments,M,B,tabComp] = createODEsystem_v7(n,p,m)
%from the previous previous version, we add a boolean variable VTunderART.
%from the previous version (v_6), I add possibility to use a kit if
%symptomatic for an infection that is not in the kit.
%create the ode system v_7

%n = number of SICTPrEP (min=?, max=1)
%p = number of SEIIIS (min=?,max=1)
%m = number of SEIIS (min=?, max m=4)

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

syms VTunderART %boolean, 1 if VT can be practice under ART, 0 otherwise

if (m==0 && n==1 && p==1)
    warning('attention a certains rho mis a 0')
    rho_sc=0;rho_scg=0;rho_hscg=0;
    rho_hg=0; rho_hc=0;
    rho_hsc=0; rho_hsg=0;
    rho_sg=0; rho_hcg=0;
    rho_g=0;
    rho_c=0;
end
if (m==1 && n==1 && p==0)
    warning('attention a certains rho mis a 0')
    rho_sc=0;rho_scg=0;rho_hscg=0;
    rho_hg=0;
    rho_hsc=0; rho_hsg=0;
    rho_sg=0; rho_hcg=0;
    rho_hs=0;
    rho_s=0;
    rho_g=0;
end


%(the ODE systems in the matlab functions have been created without considering the
% simplifications below)
gamma1s = 0;
nus     = 0;

nDis = m+n+p;       %number of infections in the model
dis = ["HIV","syph","Ct","Ng"];
dis = dis([n==1,p==1,m>=1,m==2]); warning('verifier ici')
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
%% To adapt the code for a k-disease model, k<4
% this adds columns for the other infections (the ones not considered in the model) and make them susceptible
if m==1 && n==1 && p==0
    tabComp.Ng = repmat("S",nbCompartments,1);
    tabComp.syph = repmat("S",nbCompartments,1);
elseif m==0 && p==1 && n==1
    tabComp.Ct = repmat("S",nbCompartments,1);
    tabComp.Ng = repmat("S",nbCompartments,1);
elseif m==1 && p==0 && n==1
    tabComp.syph = repmat("S",nbCompartments,1);
    tabComp.Ng = repmat("S",nbCompartments,1);
end

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
            [~, idx] = ismember(table2array(T(:,inf)), boxesInf');
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
    warning('tester ici, line 118, B, sans HIV')
    B(sum(table2array(tabComp(:,1:(m+p)))==repmat("S",1,n+m+p),2)==n+m+p) = b;
end


%% Adds voluntary testing rates rho_xxx
%Symptomatic individuals for STIs test because of symptoms and not with
%voluntary testing.
%Individuals on PrEP (P,Ip,Cp) test only routiely, not voluntarily.

% For all combination of infections *except if individuals are on PrEP*,
% we need to add combined voluntary testing rate
% e.g. for a given box *I(H)E(s)IA(C)E(G)*
% HIV and syph tested together, then combined testing leads to S(H)S(syph),
% Ct tested alone, targeted testing goes to S(C) (already accounted in the matrix M)
% Ng ---------------------------------------S(G) (idem)

%The following table shows which kit test for which infection
CVTcomb = table({["HIV"];["syph"];["Ct"];["Ng"];["HIV","syph"];["HIV","Ct"];["HIV","Ng"];["syph","Ct"];["syph","Ng"];["Ct","Ng"];...
    ["HIV","syph","Ct"];["HIV","syph","Ng"];["HIV","Ct","Ng"];["syph","Ct","Ng"];...
    ["HIV","syph","Ct","Ng"]},...
    [1;0;0;0;1;1;1;0;0;0;1;1;1;0;1],...
    [0;1;0;0;1;0;0;1;1;0;1;1;0;1;1],...
    [0;0;1;0;0;1;0;1;0;1;1;0;1;1;1],...
    [0;0;0;1;0;0;1;0;1;1;0;1;1;1;1],...
    [rho_h;rho_s;rho_c;rho_g;rho_hs;rho_hc;rho_hg;rho_sc;rho_sg;rho_cg;rho_hsc;rho_hsg;rho_hcg;rho_scg;rho_hscg],...
    [1:15]',...
    'VariableNames',{'kit','isHIV','issyph','isCt','isNg','rho','kitNo'});

CVTcombHIV = CVTcomb(logical(CVTcomb.isHIV),:);
CVTcombS = CVTcomb(logical(CVTcomb.issyph),:);
CVTcombCt = CVTcomb(logical(CVTcomb.isCt),:);
CVTcombNg = CVTcomb(logical(CVTcomb.isNg),:);

CVTcombHIV_without_Ct = CVTcombHIV(CVTcombHIV.isCt==0,:);
CVTcombHIV_without_Ng = CVTcombHIV(CVTcombHIV.isNg==0,:);
CVTcombHIV_without_Ct_and_Ng = CVTcombHIV_without_Ct(CVTcombHIV_without_Ct.isNg==0,:);
CVTcombS_without_Ct = CVTcombS(CVTcombS.isCt==0,:);
CVTcombS_without_Ng = CVTcombS(CVTcombS.isNg==0,:);
CVTcombS_without_Ct_and_Ng = CVTcombS_without_Ct(CVTcombS_without_Ct.isNg==0,:);
CVTcombCt_without_Ng = CVTcombCt(CVTcombCt.isNg==0,:);
CVTcombNg_without_Ct = CVTcombNg(CVTcombNg.isCt==0,:);


if n==1
    compSSSS = tabComp(sum(table2array(tabComp(:,1:n+m+p))==repmat("S",1,n+m+p),2)==n+m+p,:);
    compPSSS = tabComp(sum(table2array(tabComp(:,1:n+m+p))==["P",repmat("S",1,m+p)],2)==n+m+p,:);
    compTSSS = tabComp(sum(table2array(tabComp(:,1:n+m+p))==["T",repmat("S",1,m+p)],2)==n+m+p,:);
else
    error('a faire')
end

VT_table = zeros(nbCompartments,size(CVTcomb,1));
%VT_table is a 560x15 matrix with VT_table(i,j)>0 if voluntary testing
%applied in the compartment with the kit j.
if(1)
for i=1:nbCompartments
    currentComp = tabComp(i,:);
    if sum(currentComp.no==[compSSSS.no,compPSSS.no,compTSSS.no])>0 %not fully (non infected/infectious)
        %do nothing, not VT testing in these compartments
    else
        if currentComp.HIV == "P" || currentComp.HIV == "Ip" || currentComp.HIV == "Cp"
            % no voluntary testing under PrEP
        else
            if currentComp.Ct == "IS" || currentComp.Ng == "IS"        
                %voluntary testing only if individuals are asymptomatic for
                %the kit considered
                if currentComp.Ct == "IS" && currentComp.Ng == "IS"
                    if any(currentComp.HIV==["I","C"])    
                        VT_table(i,CVTcombHIV_without_Ct_and_Ng.kitNo) = VT_table(i,CVTcombHIV_without_Ct_and_Ng.kitNo)+1;
                    end
                    if any(currentComp.syph==["E","I1","I2","I3"])
                        VT_table(i,CVTcombS_without_Ct_and_Ng.kitNo) = VT_table(i,CVTcombS_without_Ct_and_Ng.kitNo)+1;
                    end
                elseif currentComp.Ct == "IS" && ~(currentComp.Ng == "IS")
                    if any(currentComp.HIV==["I","C"])    
                        VT_table(i,CVTcombHIV_without_Ct.kitNo) = VT_table(i,CVTcombHIV_without_Ct.kitNo)+1;
                    end
                    if any(currentComp.syph==["E","I1","I2","I3"])
                        VT_table(i,CVTcombS_without_Ct.kitNo) = VT_table(i,CVTcombS_without_Ct.kitNo)+1;
                    end
                    if any(currentComp.Ng==["E","IA"])
                        VT_table(i,CVTcombNg_without_Ct.kitNo) = VT_table(i,CVTcombNg_without_Ct.kitNo)+1;
                    end        
                elseif ~(currentComp.Ct == "IS") && currentComp.Ng == "IS"
                    if any(currentComp.HIV==["I","C"])    
                        VT_table(i,CVTcombHIV_without_Ng.kitNo) = VT_table(i,CVTcombHIV_without_Ng.kitNo)+1;
                    end
                    if any(currentComp.syph==["E","I1","I2","I3"])
                        VT_table(i,CVTcombS_without_Ng.kitNo) = VT_table(i,CVTcombS_without_Ng.kitNo)+1;
                    end
                    if any(currentComp.Ct==["E","IA"])
                        VT_table(i,CVTcombCt_without_Ng.kitNo) = VT_table(i,CVTcombCt_without_Ng.kitNo)+1;
                    end 
                end
            else
                if any(currentComp.HIV==["I","C"])
                    VT_table(i,CVTcombHIV.kitNo) = VT_table(i,CVTcombHIV.kitNo)+1;
                end
                if any(currentComp.syph==["E","I1","I2","I3"])
                    VT_table(i,CVTcombS.kitNo) = VT_table(i,CVTcombS.kitNo)+1;
                end
                if any(currentComp.Ct==["E","IA"])
                    VT_table(i,CVTcombCt.kitNo) = VT_table(i,CVTcombCt.kitNo)+1;
                end
                if any(currentComp.Ng==["E","IA"])
                    VT_table(i,CVTcombNg.kitNo) = VT_table(i,CVTcombNg.kitNo)+1;
                end
                if any(currentComp.HIV==["T"])
                    VT_table(i,CVTcombHIV.kitNo) = VT_table(i,CVTcombHIV.kitNo)+1000;
                end
            end
        end
    end
end

HIVflow_hiv = table( ["S";"I";"C";"P";"Ip";"Cp";"T"], ["S";"T";"T";"P";"T";"T";"T"], 'VariableNames',{'outState','inState'});

for i=1:nbCompartments
    currentComp = tabComp(i,:);
    for k=1:15
        if VT_table(i,k)>=1
            % Finding the reception compartment of the rate rho
            currentKit = CVTcomb.kit(k);
            otherFixedStates = currentComp(:,setdiff(dis,currentKit{:}));
            receptComp = otherFixedStates;
            for inf=currentKit{:}
                if inf ~= "HIV"
                    receptComp.(inf) = "S";
                else  %si l'infection est HIV
                    receptComp.HIV = HIVflow_hiv(HIVflow_hiv.outState==currentComp.HIV,:).inState;
                end
            end
            receiverComp = innerjoin(tabComp,receptComp);
            
            if VT_table(i,k)<1000
                M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no) - min(VT_table(i,k),1)*CVTcomb(k,:).rho ; %min(...,1) to not count twice VT rates
                M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+min(VT_table(i,k),1)*CVTcomb(k,:).rho;
            else
                M(currentComp.no,currentComp.no) = M(currentComp.no,currentComp.no) - VTunderART*min(VT_table(i,k),1)*CVTcomb(k,:).rho ; %min(...,1) to not count twice VT rates
                M(receiverComp.no,currentComp.no) = M(receiverComp.no,currentComp.no)+VTunderART*min(VT_table(i,k),1)*CVTcomb(k,:).rho;
            end
        end
    end
end
end
%% Adding mandatory/recommended testing for STIs under PrEP and ART
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

end

