%clear all; 
close all;
addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
addpath('MAIN')
setupOnce;

%couleurs def
rougeHIV=[215, 0, 0]/255;
jauneS  =[250, 215, 0]/255;
bleuCt  =[56, 57, 186]/255;
vertNg  =[43, 152, 38]/255;

%----------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
%----------------%


% Lecture des paramètres
log_path = 'C:/Users/Moi/Desktop/Temporaire/tests';
ampl_models_dir = [pwd,'\MAIN\AMPL_models\'];
pathParam = ['.\ParameterAnalysis\paramSets_',num2str(paramNo),'\round_',num2str(roundNo),'\'];
paramH = readtable([pathParam, 'allParametersSets_HIV.txt']);
paramS = readtable([pathParam, 'allParametersSets_syphilis.txt']);
paramC = readtable([pathParam, 'allParametersSets_Ct.txt']);
paramG = readtable([pathParam, 'allParametersSets_Ng.txt']);
ID_ech = paramH(paramH.p==pHIV & paramH.nbEch==nbEch,:).IDech_id;
paramTab{1} = table2struct(paramC(paramC.IDech_id==ID_ech,:));
paramTab{2} = table2struct(paramG(paramG.IDech_id==ID_ech,:));
paramTab{3} = table2struct(paramH(paramH.IDech_id==ID_ech & paramH.p==pHIV,:));
paramTab{4} = table2struct(paramS(paramS.IDech_id==ID_ech,:));
mu=paramTab{3}.mu; b = paramTab{3}.pi;

verbose=0;
paramSolver.tolP0 = 0.5e-4;
paramSolver.maxBndAlpha=20;
paramSolver.nbRelanceMax=5;
paramSolver.timeLimit = 20; %seconds
paramSolver.iterMaxDicho = 20;
paramSolver.tolAlpha = 1e-4;
paramSolver.method_alpha = 'dicho';
paramSolver.timeSolver = 20;
%opt.TolP0=0.5e-4;

pathRes = ['.\ParameterAnalysis\results_',num2str(paramNo),'\_round_',num2str(roundNo),'\'];
paramRho = table2struct(readtable([pathRes ,'paramRho.txt']));

paramTab{1}
paramTab{2}
%%
elimCost = readtable([pathRes,'elimCosts_f_',num2str(f),'.txt']);

opts = detectImportOptions([pathRes,'tabAlpha.txt']);
opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
tabAlpha  = readtable([pathRes,'tabAlpha.txt'],opts);
alphas = tabAlpha(tabAlpha.IDech==ID_ech & tabAlpha.p==pHIV,:);

costsOfElim = elimCost(elimCost.IDech==ID_ech & elimCost.p==pHIV,:);

%%
kit = {'HIV','syphilis','Ct','Ng'};%,'Ct','Ng'};
mod = 'hscg_hscg';
costsOfElim = elimCost(elimCost.IDech==ID_ech & ismember(elimCost.kit,indexKit(kit)),:);

%% ALPHA 
addpath('MAIN')
tic
paramSolver.method_alpha = 'dicho';
[alpha,Palpha,ES,msg,elim_i] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir);
toc

paramSolver.method_alpha = 'pas_dicho';
tic
[alpha,Palpha,ES,msg,elim_i] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir);
toc

%% Evolution de la prevalence / Calcul de la prevalence (reelle, percue)
f   = 3;
mod = 'hscg_hscg';
kit = {'HIV','syphilis'};%,'Ct','Ng'};
paramSolver.varToChange = 'rho_sc';
paramRho_bis = paramRho;
IA=[];
IS=[];
vecRho = [0,0.05,0.1:0.005:0.2,0.4:0.005:0.5,1];
Pk_r = zeros(1,length(vecRho)); Pk_p = Pk_r;
Ph = [];Pc = [];Pg = [];Ps = [];Ph_p = [];
ampl=0;
i=1;
for rho = vecRho
    paramRho_bis.(paramSolver.varToChange) = rho;
    [P,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl);
    IA(i) = sum(ES([71:105,211:245,351:385,491:525]))/(b/mu); 
    IS(i) = sum(ES([106:140,246:280,386:420,526:560]))/(b/mu);    
    [Pk_r(i),Pk_p(i),Ph_p(i)] = P_kit_v3(ES,kit,f);
    Ph(i) = P.h;
    Pc(i) = P.c;
    Ps(i) = P.s;
    Pg(i) = P.g;
    disp([num2str(i),'/',num2str(length(vecRho))])
    i=i+1;
end

%% showing prevalences (computed in the last section)
close all;
fig = figure()
plot(vecRho,Pk_r,'k-','LineWidth',2)
hold on
plot(vecRho,Pk_p,'k--','LineWidth',2)
plot(vecRho,Ph,'Color',rougeHIV,'LineWidth',2)
plot(vecRho,Ph_p,'--','Color',rougeHIV,'LineWidth',2)
plot(vecRho,Ps,'Color',jauneS,'LineWidth',2)
plot(vecRho,Pc,'Color',bleuCt,'LineWidth',2)
plot(vecRho,Pg,'Color',vertNg,'LineWidth',2)
lgd = legend('$\Pi_{k}$','$\Pi_{k}^b$','$\Pi_h$','$\Pi_{h}^b$','$\Pi_s$','$\Pi_c$','$\Pi_g$','Interpreter','latex');
%lgd = legend(['$\Pi_{',indexKit(kit),'}$'],'$\Pi_h$','$\Pi_s$','$\Pi_c$','$\Pi_g$','Interpreter','latex','EdgeColor','w');
set(lgd,'FontSize',18);
%plot(vecRho([1,end]),paramSolver.tolP0*ones(2,1),'k','DisplayName','tolP')
%title('Prevalence','Interpreter','latex')
%xlabel(['\',paramSolver.varToChange])
xlabel(['$\rho_{',indexKit(kit),'}$'],'FontSize',18,'Interpreter','latex')
ylabel('Prevalence','Interpreter','latex','FontSize',18)
xlim([0,max(vecRho)])
set(fig,'position',[100,100,370,300])
%set(fig,'PaperSize',[20 10]); %set the paper size to what you want  
%print(fig,[log_path,'/prevalence_',indexKit(kit),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'.pdf'],'-dpdf') 
saveas(fig,[log_path,'/prevalence_',indexKit(kit),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'.pdf'])


%%
%utility en fonction de rho et c
figure() 
hold on

for c=[-1,0,0.5]%[-0.5,-0.6,-0.7,-0.8]           % [-0.2,-0.19, -0.11,-0.05,0]
    U = vecRho.*(Pk_r-c);   
    plot(vecRho,U,'-','DisplayName',['c=',num2str(c)],'LineWidth',2);
    
    U = vecRho.*(Pk_p-c);    %quand f=1
    plot(vecRho,U,'--','DisplayName',['c=',num2str(c)],'LineWidth',2);
end
plot(vecRho, zeros(length(vecRho),1),'k','HandleVisibility','off')
legend()
title(['ID=',num2str(ID_ech),' p=', num2str(pHIV),' kit=',indexKit(kit)])
%xlabel(paramSolver.varToChange)
%xlabel(['\rho_{hscg}'])

%plot(0.2617*ones(2,1),[0,0.25],'k','HandleVisibility','off')
ylimMax = 1;
%plot(alphas(alphas.p==pHIV,:).hscg_hscg*ones(2,1),[0,ylimMax],'Color',[0.5 0.5 0.5],'HandleVisibility','off')
%plot(alphas(alphas.p==pHIV,:).hsc_hsc*ones(2,1),[0,ylimMax],'Color',[0.5 0.5 0.5],'HandleVisibility','off')
%plot(alphas(alphas.p==pHIV,:).hs_hs*ones(2,1),[0,ylimMax],'Color',[0.5 0.5 0.5],'HandleVisibility','off')
%plot(alphas(alphas.p==pHIV,:).h_h*ones(2,1),[0,ylimMax],'Color',[0.5 0.5 0.5],'HandleVisibility','off')

%plot([0,0.7],0.0104686*ones(2,1),'k--','HandleVisibility','off')
%plot([0,0.7],0.0103936*ones(2,1),'k-.','HandleVisibility','off')

%% utilité par morceaux
figure() 
Y = U
%vecC = [-0.19249,-0.10216,-0.05624];
vecC = [-0.3:0.1:0,0.05];
ey=0.01;
j=0; newmax=0;
for c = vecC
    hold on
    U = vecRho.*(Pk_r-c);
    Y=U;
    transp =  1 - (length(vecC)-j)*0/length(vecC);
    plot(vecRho,Y,'Color',[0.9290 0.6940 0.1250 transp],'LineWidth',3);
    plot(vecRho(Pg<=paramSolver.tolP0),Y(Pg<=paramSolver.tolP0),'Color',[0.4660 0.6740 0.1880 transp],'LineWidth',3);
    plot(vecRho(Pc<=paramSolver.tolP0),Y(Pc<=paramSolver.tolP0),'Color',[0 0.4470 0.7410 transp],'LineWidth',3);
    plot(vecRho(Ps<=paramSolver.tolP0),Y(Ps<=paramSolver.tolP0),'Color',[0.4940 0.1840 0.5560 transp],'LineWidth',3);
    plot(vecRho(Ph<=paramSolver.tolP0),Y(Ph<=paramSolver.tolP0),'Color',[0 0 0 transp],'LineWidth',3);
    text(0.95*vecRho(end), U(end)+ey,['c=',num2str(round(c,2))])
    
    %ad horizontal line for the max
    plot(vecRho(1):vecRho(end),max(Y(Ph>=paramSolver.tolP0))*ones(1,2),'--','Color',[0.8 0.8 0.8],'HandleVisibility','off','LineWidth',1);
    newmax = max([newmax,max(U)]);
    j=j+1;
end
plot(vecRho, zeros(length(vecRho),1),'k','HandleVisibility','off')
title(['ID=',num2str(ID_ech),' round=',num2str(roundNo),' p=', num2str(pHIV),' kit=',indexKit(kit),' b=',num2str(f)])
xlabel(['\rho_{hscg}'],'FontSize',12)

plot(alphas(alphas.p==pHIV,:).hscg_hscg*ones(2,1),[0,1.2*newmax],'Color',[0.8 0.8 0.8],'HandleVisibility','off','LineWidth',0.8)
plot(alphas(alphas.p==pHIV,:).hsc_hsc*ones(2,1),[0,1.2*newmax],'Color',[0.8 0.8 0.8],'HandleVisibility','off','LineWidth',0.8)
plot(alphas(alphas.p==pHIV,:).hs_hs*ones(2,1),[0,1.2*newmax],'Color',[0.8 0.8 0.8],'HandleVisibility','off','LineWidth',0.8)
plot(alphas(alphas.p==pHIV,:).h_h*ones(2,1),[0,1.2*newmax],'Color',[0.8 0.8 0.8],'HandleVisibility','off','LineWidth',0.8)

xticks(round([0,alphas(alphas.p==pHIV,:).hscg_hscg,alphas(alphas.p==pHIV,:).hsc_hsc,alphas(alphas.p==pHIV,:).hs_hs,alphas(alphas.p==pHIV,:).h_h,1],2))
%ylim([0,1.2*newmax])

%% Calcul de rhohat pour une valeur de c
c=0.1;
newmod=1;
mySeed=1;
ampl_c=0;
paramSolver.verbose=1;
[rhohat,Cval,ES,msg_tot,ampl_c] = findRhohat_kit_v9(paramTab,mu,b,f,paramRho,kit,c,alphas,paramSolver.verbose,mySeed,log_path,paramSolver,ampl_c,newmod);

rhohat
ampl_c.close();

paramRho_bis.(paramSolver.varToChange) = rhohat.tot;
ampl=0;
[P,ES,msgSolver,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,paramSolver.verbose,paramSolver,log_path,ampl_models_dir,ampl);

%(P.h-0.5e-4)

%% Cnn
paramsC         = paramSolver;
paramsC.sup     = 1;
paramsC.inf     = -5;
paramsC.tolC    = 1e-3;
paramsC.iterMax = 60;
afficherOutput  = 0;
paramSolver.verbose       = 0;
paramSolver.sup_bnd_alpha = 20;
paramSolver.timeSolver    = 20;
paramSolver.tolP0         = paramSolver.tolP0;

tic;
apriori = []; %[-0.1923, -0.1021, -0.0564];
[cnn,tabP,msg] = find_Cnn_kit_2(paramTab,mu,b,paramRho,f,kit,paramsC,alphas,afficherOutput,log_path,paramSolver,apriori);
toc


%% Representation de l'impact de rho_hg sur rho_sc'
kit={'HIV','Ng'};
mod='hscg';
vecAlpha=[];
paramRho_bis = paramRho;
vecRho_sc = 0:0.1:1;
for rho_sc=vecRho_sc
    disp(rho_sc)
    paramRho_bis.rho_sc =  rho_sc;
    [alpha,Palpha,ES,msg,elim_i] = findAlpha_v4(paramTab,paramRho_bis,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir);
    disp(alpha)
    vecAlpha = [vecAlpha,alpha];
end

figure()
%plot(vecRho_hs, vecAlpha, 'k--','LineWidth',0.5)
hold on
plot(vecRho_sc, round(vecAlpha,3), 'k', 'LineWidth',1.5)
%hold on
%plot(vecRho_hs(6:end), vecAlpha(6:end), 'k','LineWidth',1.5)
ylim([0,0.1])

title(['ID=',num2str(ID_ech),' p=', num2str(pHIV),' kit=',indexKit(kit)])
xlabel('$\rho_{sc}$','interpreter','latex','FontSize',15)
ylabel('$\rho_{hg}^\prime$','interpreter','latex','FontSize',15)

%% Representation de l'impact de rho_hg sur \hat\rho_sc et sur \rho_sc'
kit={'syphilis','Ct'}; i=1;
c=0; paramSolver.verbose=0;mySeed=1;ampl_c=0;newmod=1;
paramSolver.tolAlpha = 1e-3;
k_bis = 'hg';
mod='hscg';
paramRho_bis = paramRho;
vecRho = 0:0.1:1;
vecAlpha=zeros(length(vecRho),3); vecRhohat=zeros(length(vecRho),1);
for rho=vecRho
    disp(rho)
    paramRho_bis.(['rho_',k_bis]) =  rho;
    [infElim,~,~] = checkFeasability_v2(paramTab,paramRho_bis,b,mu,f,mod,verbose,paramSolver,log_path,ampl_models_dir);
    if infElim.h
        [alphas.sc_sc,~,~,~,~] = findAlpha_v4(paramTab,paramRho_bis,mu,b,f,kit,'hscg',verbose,log_path,paramSolver,ampl_models_dir);
        alphas.elim_h=1;
    else
        [alphas.sc_hsc,~,~,~,~] = findAlpha_v4(paramTab,paramRho_bis,mu,b,f,kit,'hscg',verbose,log_path,paramSolver,ampl_models_dir);
        alphas.elim_h=0;
    end
    [alphas.s_s,~,~,~,~] = findAlpha_v4(paramTab,paramRho_bis,mu,b,f,{'syphilis'},'hscg',verbose,log_path,paramSolver,ampl_models_dir);
    [alphas.c_c,~,~,~,~] = findAlpha_v4(paramTab,paramRho_bis,mu,b,f,{'Ct'},'hscg',verbose,log_path,paramSolver,ampl_models_dir);
    
    [rhohat,~,~,msg_tot,ampl_c] = findRhohat_kit_v9(paramTab,mu,b,f,paramRho_bis,kit,c,alphas,paramSolver.verbose,mySeed,log_path,paramSolver,ampl_c,newmod);

    vecAlpha(i,1) = alphas.sc_sc;
    vecAlpha(i,2) = alphas.s_s;
    vecAlpha(i,3) = alphas.c_c;
    vecRhohat(i)  = rhohat.tot;
    i=i+1;
end
ampl_c.close();


figure()
hold on
plot(vecRho, vecRhohat, 'k--', 'LineWidth',1.5)
plot(vecRho(1:6), vecRhohat(1:6), 'k', 'LineWidth',1.5)
plot(vecRho(7:end), vecRhohat(7:end), 'k', 'LineWidth',1.5)

title(['ID=',num2str(ID_ech),' p=', num2str(pHIV),' kit=',indexKit(kit), ' b=',num2str(f),' c=',num2str(c)],'interpreter','latex')
xlabel('$\rho_{hg}$','interpreter','latex','FontSize',15)
ylabel(['$\hat\rho_{',indexKit(kit),'}$'],'interpreter','latex','FontSize',15)

%%
figure()
plot(vecRho, round(vecAlpha(:,1),3), 'k--', 'LineWidth',1.5)
hold on
plot(vecRho, round(vecAlpha(:,1),3), 'k', 'LineWidth',1.5)
plot(vecRho, round(vecAlpha(:,1),3), 'k', 'LineWidth',1.5)

plot(vecRho, round(vecAlpha(:,2),3), '--','Color',jauneS, 'LineWidth',1.5)
plot(vecRho, round(vecAlpha(:,2),3), 'Color',jauneS, 'LineWidth',1.5)
plot(vecRho, round(vecAlpha(:,2),3), 'Color',jauneS, 'LineWidth',1.5)

plot(vecRho(:), round(vecAlpha(:,3),3), '--','Color',bleuCt, 'LineWidth',1.5)
plot(vecRho, round(vecAlpha(:,3),3), 'Color',bleuCt, 'LineWidth',1.5)
plot(vecRho, round(vecAlpha(:,3),3), 'Color',bleuCt, 'LineWidth',1.5)


%% U_sc as a function of rho_hg

vecRho_sc = 0:0.01:0.8; n = length(vecRho_sc);
vecRho_hg = 0.6; m = length(vecRho_hg);
Psc = zeros(n,m);
Phg = zeros(n,m);
Ph  = zeros(n,m);
Ps  = zeros(n,m);
Pc  = zeros(n,m);
Pg  = zeros(n,m);
for i=1:n
    rho_sc = vecRho_sc(i);
    disp([num2str(i),'/',num2str(n)])
    for j=1:m
        rho_hg = vecRho_hg(j);
        paramRho_bis.rho_sc = rho_sc;
        paramRho_bis.rho_hg = rho_hg;
        
        ampl=0;
        [~,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,paramSolver.verbose,paramSolver,log_path,ampl_models_dir,ampl);
        ampl.close();
        Psc(i,j) = P_kit_v3(ES,{'syphilis','Ct'},1);
        Phg(i,j) = P_kit_v3(ES,{'HIV','Ng'},1);
        Ps(i,j)  = P_kit_v3(ES,{'syphilis'},1);
        Pc(i,j)  = P_kit_v3(ES,{'Ct'},1);
        Ph(i,j)  = P_kit_v3(ES,{'HIV'},1);
        Pg(i,j)  = P_kit_v3(ES,{'Ng'},1);
    end
end
%% en 1D
figure()
plot(vecRho_sc,vecRho_sc'.*(Psc-c),'LineWidth',2)
xlabel('$\rho_{sc}$','interpreter','latex','FontSize',15)
title(['$U_{sc}$ with $\rho_{hg}=$',num2str(rho_hg)],'interpreter','latex','FontSize',15)
%% 3D
figure(1)
s=surf(vecRho_hg,vecRho_sc,Psc);
s.EdgeColor = 'none';
xlabel('$\rho_{hg}$','interpreter','latex','FontSize',16)
ylabel('$\rho_{sc}$','interpreter','latex','FontSize',16)
title(['$\Pi_{sc}$, ID',num2str(ID_ech),' p=',num2str(pHIV)],'interpreter','latex','FontSize',14)
colorbar;

figure(2)
s=surf(vecRho_hg,vecRho_sc,Phg);
s.EdgeColor = 'none';
xlabel('$\rho_{hg}$','interpreter','latex','FontSize',16)
ylabel('$\rho_{sc}$','interpreter','latex','FontSize',16)
title(['$\Pi_{hg}$, ID',num2str(ID_ech),' p=',num2str(pHIV)],'interpreter','latex','FontSize',14)
colorbar;

%% U 3D
c=0;
Usc = vecRho_sc'.*(Psc-c);
Uhg = vecRho_hg.*(Phg-c);

figure(3)
s=surf(vecRho_hg,vecRho_sc,Usc);
s.EdgeColor = 'none';
xlabel('$\rho_{hg}$','interpreter','latex','FontSize',16)
ylabel('$\rho_{sc}$','interpreter','latex','FontSize',16)
title(['$U_{sc}$, ID',num2str(ID_ech),' p=',num2str(pHIV), ' c=',num2str(c)],'interpreter','latex','FontSize',14)
colorbar;

figure(4)
s=surf(vecRho_hg,vecRho_sc,Uhg);
s.EdgeColor = 'none';
xlabel('$\rho_{hg}$','interpreter','latex','FontSize',16)
ylabel('$\rho_{sc}$','interpreter','latex','FontSize',16)
title(['$U_{hg}$, ID',num2str(ID_ech),' p=',num2str(pHIV),' c=',num2str(c)],'interpreter','latex','FontSize',14)
colorbar;

%%
Us = vecRho.*(Ps-c);
figure(3)
plot(vecRho,Us)

%%
vecRho = 0.:0.05:1; 
i=1; P=zeros(length(vecRho),1);
for rho=vecRho
    disp(rho)
    [P(i),P12,PHIV,PIST,ES] = P12_SICTPSEIIS_v7(param1,param2,mu,b,paramRho,rho,f,'knitro-ampl',opt);
    disp(P(i))
    i=i+1;
end
plot(vecRho,P)


%%
[alpha,P_alpha,ES,msg,~] = findAlpha_v3(paramTab,paramRho,mu,b,f,kit,mod,0,log_path,paramSolver,ampl_models_dir);

figure()
vecC = -1:0.1:1;
rhohat = zeros(length(vecC),1)
i=1;
% il faut borner rho dans le calcul de minU par le plus des alpha
% et le plus petit des alhas doit etre determiner avec findAlpha
% ou un truc comme ça
opt.tolP0=1e-5;
opt.verbose=0;
opt.up_bnd_alpha=alpha;
paramTab = {param1,param2};  
for c=vecC
    disp(['c=',num2str(c)])     
    [rhohat21,Cval,ES,msg] = fminU_knitro_v7(kit,paramTab,paramRho,mu,b,c,opt);
%         if rhohat21~=0
%             Pval(i) = c-Cval/rhohat21
%         else
%             Pval(i) = 1;
%         end
%         disp(rhohat21)
%         disp(c)
    rhohat(i)=rhohat21;
    disp(rhohat21)
    i=i+1;         
end

plot(vecC,rhohat,'DisplayName','simu ','LineWidth',1)
hold on
legend()


figure()
U = vecRho.*(P'-c)
plot(vecRho,U)

%%


kit={'HIV','Ng'};
mod = 'hg_hg';
[alpha,P,ES,msg] = findAlpha_v3(paramTab,paramRho,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir);
P


[~,~,alphaS] = Rp_SEIIIS_v4(paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.tau,...
    paramTab{4}.nu,paramTab{4}.gamma10,paramTab{4}.theta,paramTab{4}.gamma30,mu,b,0);

[~,~,alphaH,~,~] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
    paramTab{3}.theta0,paramTab{3}.gamma0,paramTab{3}.sigma,paramTab{3}.zeta,...
    paramRho.eta_h_prep,paramTab{3}.p,mu,b,0)



%%
opts = detectImportOptions([myPath,'tabAlpha_last.txt']);
opts.VariableTypes(19:26)={'char'};
tabAlpha   = readtable([myPath,'tabAlpha_last.txt'],opts);


id = tabAlpha.IDech(1);
p=0

alphas = tabAlpha(tabAlpha.IDech==id & tabAlpha.p==p,:);

paramTab{1} = table2struct(paramC(paramC.IDech_id==id,:));
paramTab{2} = table2struct(paramG(paramG.IDech_id==id,:));
paramTab{3} = table2struct(paramH(paramH.IDech_id==id & paramH.p==p,:));
paramTab{4} = table2struct(paramS(paramS.IDech_id==id,:));
mu = paramTab{1}.mu;
b = paramTab{1}.pi;
verbose = 0;
log_path = 'C:/Users/Moi/Desktop/Temporaire/tests';
paramSolver.tolP0  = 0.5e-4;%;
paramSolver.maxBndAlpha = 10;%;
paramSolver.nbRelanceMax = 5;%;
paramSolver.timeLimit = 20;%; 
ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';

f=1;
mod = 'scg_hscg';
kit={'syphilis','Ct','Ng'};
[alpha,Palpha,ES,msg,elim_i] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir)
paramRho.rho_hs = alpha;
[P,ES,msgSolver] = P_mod_v7(paramTab,paramRho,b,mu,f,'hscg',verbose,paramSolver,log_path,ampl_models_dir);
P_k = P_kit_v2(ES,kit);
P.h
P.c
P.s;
P.g;
    
[rhohat,cval,ES,status] = findRhohat_kit_v8(paramTab,mu,b,f,paramRho,kit,c,alphas,verbose,mySeed,log_path,paramSolver);

%%
%close all;
log_path = 'C:/Users/Moi/Desktop/Temporaire/tests';
ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
verbose=0;

paramSolver.tolP0  = 0.5e-4;%;
paramSolver.maxBndAlpha = 50;%;
paramSolver.nbRelanceMax = 5;%;
paramSolver.timeLimit = 20;%; 
paramSolver.verbose = verbose;%; 
paramSolver.method_alpha = '';
paramSolver.sup_bnd_alpha = paramSolver.maxBndAlpha;
f=1;

%attention le mod doit prendre en compte s'il y a HIV ou non
mod = 'hscg'; %mettre p à 0 devrait suffire.
kit = {'HIV','Ct','Ng'};
paramSolver.varToChange = 'rho_hcg';

paramRho_bis = paramRho;
IA=[];
IS=[];
vecRho = 0.:0.01:0.9;
P_k=zeros(1,length(vecRho));
Ph=[];Pc=[];Pg=[];Ps=[];
ampl=0;
%paramTab{3}.p = p;
i=1;
for rho = vecRho
    paramRho_bis.(paramSolver.varToChange) = rho;
    [P,ES,msgSolver,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,f,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl);
    IA(i) = sum(ES([71:105,211:245,351:385,491:525]))/(b/mu);
    IS(i) = sum(ES([106:140,246:280,386:420,526:560]))/(b/mu);
    
    P_k(i) = P_kit_v2(ES,kit);
    Ph(i) = P.h;
    Pc(i) = P.c;
    Ps(i) = P.s;
    Pg(i) = P.g;
    i=i+1
end

figure(4)
hold on
for c=-0.03:-0.005:-0.04
    U = vecRho.*(P_k-c);
    plot(vecRho,U,'DisplayName',num2str(c),'LineWidth',1.5);
end
plot(vecRho, zeros(length(vecRho),1),'k')
legend()

figure(2)
%plot(vecRho, P_k)
hold on
plot(vecRho,Pc,'b','LineWidth',1.5)
plot(vecRho,Ps,'y','LineWidth',1.5)
plot(vecRho,Pg,'g','LineWidth',1.5)
plot(vecRho,Ph,'r','LineWidth',1.5)
legend('Pc','Ps','Pg','Ph')
plot(vecRho([1,end]),paramSolver.tolP0*ones(2,1),'k')

Us = vecRho.*(Ps-c);
figure(3)
plot(vecRho,Us)


%%
[alpha,P_alpha,ES,msg,elim_i] = findAlpha_v3(paramTab,paramRho,mu,b,f,kit,mod,0,log_path,paramSolverAlpha,ampl_models_dir);

mySeed = 1;
[rhohat,cvals,ES,status] = findRhohat_kit_v8(paramTab,mu,b,f,paramRho,kit,c,alphas,verbose,mySeed,log_path,paramSolver);
%[rhohat,Cval,ES,msg_tot] = findRhohat_kit_v8_by_dis(paramTab,mu,b,f,paramRho,kit,'h',c,alphas,verbose,mySeed,log_path,paramSolver);

mod_obj = 'h_hg';
paramSolver.up_bnd_alpha = paramSolver.maxBndAlpha;
paramSolver.inf_bnd_alpha = 0;
c=0;
[rhohat,Cval,ES,P,msg] = fmaxU_knitro_v7_bis(kit,mod_obj,paramTab,paramRho,mu,b,c,paramSolver);

%% Recherche de c_elim
%1125183724000100e, p=0

paramsC = paramSolver;
paramsC.sup=1;
paramsC.inf=-5;
paramsC.tolC=1e-3;
paramsC.iterMax=30;
afficherOutput=0;

createParamRho;
[cnn,tabP,msg] = find_Cnn_kit(paramTab,mu,b,paramRho,f,kit,paramsC,alphas,afficherOutput,log_path,paramSolver);





%% Evolution des prevalences des infections quand rho_h varie (3D graph)
mod = 'hscg'; %mettre p à 0 devrait suffire.
PH = [];
PS = [];
PC = [];
PG = [];
vec_rho_s = 0:0.1:1;
vec_rho_c = 0:0.1:1;
RHOS = [];
RHOC = [];
paramRho_bis = paramRho;
for i=1:length(vec_rho_s)
    disp([num2str(i),'/',num2str(length(vec_rho_s))])
    rho_s = vec_rho_s(i);
    
    paramRho_bis = paramRho;
    paramRho_bis.rho_s = rho_s;

    for j=1:length(vec_rho_c)
        ampl = 0;

        rho_c = vec_rho_c(j);
        RHOS(i,j) = rho_s;
        RHOC(i,j) = rho_c;        
        paramRho_bis.rho_c = rho_c;
        [P,ES,status,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,f,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl);
        ampl.close();
        PH(i,j) = P.h;
        PS(i,j) = P.s;
        PC(i,j) = P.c;
        PG(i,j) = P.g;
    end   
end
%%
close all

figure(1)
s=surf(RHOS,RHOC,round(PH,4))
xlabel('$\rho_s$','Interpreter','latex','FontSize', 20)
ylabel('$\rho_c$','Interpreter','latex','FontSize', 20)
zlabel('$\Pi_h(\rho_s,\rho_c)$','Interpreter','latex','FontSize', 20)
title(['$\Pi_h(\rho_s,\rho_c)$, parameter set no ',num2str(ID_ech),' p=', num2str(pHIV)],'interpreter','latex','FontSize', 15)
s.EdgeColor = 'none';
colorbar

figure(2)
s=surf(RHOS,RHOC,round(PS,4))
xlabel('$\rho_s$','Interpreter','latex','FontSize', 20)
ylabel('$\rho_c$','Interpreter','latex','FontSize', 20)
zlabel('$\Pi_s(\rho_s,\rho_c)$','Interpreter','latex','FontSize', 20)
title(['$\Pi_s(\rho_s,\rho_c)$, parameter set no ',num2str(ID_ech),' p=', num2str(pHIV)],'interpreter','latex','FontSize', 15)
s.EdgeColor = 'none';
colorbar

figure(3)
s=surf(RHOS,RHOC,round(PC,4))
xlabel('$\rho_s$','Interpreter','latex','FontSize', 20)
ylabel('$\rho_c$','Interpreter','latex','FontSize', 20)
zlabel('$\Pi_c(\rho_s,\rho_c)$','Interpreter','latex','FontSize', 20)
title(['$\Pi_c(\rho_s,\rho_c)$, parameter set no ',num2str(ID_ech),' p=', num2str(pHIV)],'interpreter','latex','FontSize', 15)
s.EdgeColor = 'none';
colorbar

figure(4)
s=surf(RHOS,RHOC,round(PG,4))
xlabel('$\rho_s$','Interpreter','latex','FontSize', 20)
ylabel('$\rho_c$','Interpreter','latex','FontSize', 20)
zlabel('$\Pi_g(\rho_s,\rho_c)$','Interpreter','latex','FontSize', 20)
title(['$\Pi_g(\rho_s,\rho_c)$, parameter set no ',num2str(ID_ech),' p=', num2str(pHIV)],'interpreter','latex','FontSize', 15)
s.EdgeColor = 'none';
colorbar


%% 
clear all; 
close all;
addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
addpath('MAIN')
log_path = 'C:/Users/Moi/Desktop/Temporaire/tests/';
fig_path = 'C:\Users\Moi\Documents\IPLESP\These\Articles\Article 2\figures\';
ampl_models_dir = [pwd,'\MAIN\AMPL_models\'];
setupOnce;

%couleurs def
rougeHIV = [215, 0, 0]/255;
jauneS   = [250, 215, 0]/255;
bleuCt   = [56, 57, 186]/255;
vertNg   = [43, 152, 38]/255;
[WGB_log_cm] = subdivisedColormap([[1,1,1];[0.08,0.08,0.08];[0,0,0]], 100, 'log');
%% Comparaison findAlpha dichotomy ou methode directe (avec probleme d'optim)
%close all; 
%clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit     = {'HIV','syphilis','Ct','Ng'};
paramSolver.varToChange = 'rho_hscg';
vecRho  = 0:0.01:1;%[0,0.05,0.1:0.005:0.2,0.4:0.005:0.5,1];
fontSize = 18;
getParameters; disp(num2str(ID_ech))
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit     = {'HIV','syphilis','Ct','Ng'};
paramSolver.varToChange = 'rho_hscg';
fontSize = 18;
getParameters; disp(num2str(ID_ech))
%----------------------------------%
k_mod = 'hscg_hscg';

%%
%Bisection method
paramSolver.iterMaxDicho = 10;
paramSolver.method_alpha = 'dicho';
[a11,P11,~,~,elim11] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit,k_mod,verbose,log_path,paramSolver,ampl_models_dir);

%Optim pb method
paramSolver.iterMaxDicho = 0;
paramSolver.method_alpha = 'pas_dicho';
paramSolver.timeLimit = 30;
verbose=1;
[a12,P12,~,~,elim12] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit,k_mod,verbose,log_path,paramSolver,ampl_models_dir);

%% prevalence avec les parametres ci-dessus
vecRho = [0,0.1,0.12:0.001:0.135,0.14:0.1:0.8];
Ph = [];Pc = [];Pg = [];Ps = [];Ph_p = [];Pk_r = zeros(1,length(vecRho));
ampl = 0;
i=1; paramRho_bis = paramRho;

for rho = vecRho
    paramRho_bis.(paramSolver.varToChange) = rho;
    [P,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl); 
    Ph(i) = P.h;
    Pc(i) = P.c;
    Ps(i) = P.s;
    Pg(i) = P.g;
    [Pk_r(i),~,~] = P_kit_v3(ES,kit,f); %prevalence of asymptomatic and undiagnosed
    %[Pk_r(i),~,~] = P_kit_v3(ES,kit,f); %prevalence of asymptomatic and undiagnosed
    disp([num2str(i),'/',num2str(length(vecRho))])
    i=i+1;
end
ampl.close();

fig = figure();
plot(vecRho,Pk_r,'k-','LineWidth',2)
hold on
plot(vecRho,Ph,'Color',rougeHIV,'LineWidth',2)
plot(vecRho,Ps,'Color',jauneS,'LineWidth',2)
plot(vecRho,Pc,'Color',bleuCt,'LineWidth',2)
plot(vecRho,Pg,'Color',vertNg,'LineWidth',2)
lgd = legend(['$\Pi_{',indexKit(kit),'}$'],'$\Pi_h$','$\Pi_s$','$\Pi_c$','$\Pi_g$','Interpreter','latex','FontSize',fontSize);
xlabel(['$\rho_{',indexKit(kit),'}$'],'FontSize',fontSize,'Interpreter','latex')
ylabel('Prevalence','Interpreter','latex','FontSize',fontSize)
yline(paramSolver.tolP0,'LineWidth',0.5)
yline(0.59532,'LineWidth',0.5)
xline(0.11187,'LineWidth',0.5)

