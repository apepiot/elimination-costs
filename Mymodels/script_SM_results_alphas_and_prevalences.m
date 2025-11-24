clear all; close all;
parametrizationNo = 48;
vecRounds         = 1:176;
rougeHIV = [215, 0, 0]/255;     cols.h = rougeHIV;
jauneS   = [250, 215, 0]/255;   cols.s = jauneS;
bleuCt   = [56, 57, 186]/255;   cols.c = bleuCt;
vertNg   = [43, 152, 38]/255;   cols.g = vertNg;
fontSize = 18;
colMat = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];...
    [0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];...
    [0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];...
    [0.6350 0.0780 0.1840]];
%/ParameterAnalysis/createParameterSets_alea.m

backupPath = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis/results_',num2str(parametrizationNo),'/'];
pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\parameterAnalysis\';
%Reading parameterSets
setsPath = ['C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/ParameterAnalysis/paramSets_',num2str(parametrizationNo)];
fig_path_sm = 'C:\Users\Moi\Documents\IPLESP\These\Articles\Drafts\Article2_v5\supplementary_material\figures\parameterAnalysis\';
%
paramSetsCt_tot=[];paramSetsNg_tot=[];paramSetsHIV_tot=[];paramSetsS_tot=[];
for roundNo = vecRounds
    disp(roundNo)
    pathFile     = [setsPath,'/round_',num2str(roundNo)];
    paramSetsCt  = readtable([pathFile,'/allParametersSets_Ct.txt']);
    paramSetsNg  = readtable([pathFile,'/allParametersSets_Ng']);
    paramSetsHIV = readtable([pathFile,'/allParametersSets_HIV']);
    paramSetsS   = readtable([pathFile,'/allParametersSets_syphilis']);
    
    paramSetsCt_tot = [paramSetsCt_tot;paramSetsCt];
    paramSetsNg_tot = [paramSetsNg_tot;paramSetsNg];
    paramSetsHIV_tot = [paramSetsHIV_tot;paramSetsHIV];
    paramSetsS_tot = [paramSetsS_tot;paramSetsS];
end

%% 0.
ech = paramSetsS_tot.beta;
round(mean(ech),2)
disp(['[',num2str(round(min(ech),2)),',',num2str(round(max(ech),2)),']'])


%% 1.a Histogram des 4 alpha (sans la prep)
%pas dans le SM
fig = figure(1);
BinWidth=0.01;
histogram(paramSetsHIV_tot.alpha_prev(paramSetsHIV_tot.p==0,:),'FaceColor',cols.h,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(paramSetsS_tot.alpha,'FaceColor',cols.s,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(paramSetsCt_tot.alpha,'FaceColor',cols.c,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(paramSetsNg_tot.alpha,'FaceColor',cols.g,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
legend('$\rho_{h,p=0}^\prime$',...
    '$\rho_{s}^\prime$',...
    '$\rho_{c}^\prime$',...
    '$\rho_{g}^\prime$',...
    'Interpreter','latex','FontSize',12,'Box','off')
xlim([0,1])
xlabel('1/year','Interpreter','latex')


%% 1.b Histogram des 4 alpha (sans la prep) avec comparaison des rhob
fig = figure(2);
BinWidth=0.01;
%adding rhob
histogram(paramSetsHIV_tot.rhob(paramSetsHIV_tot.p==0,:),'FaceColor','none','EdgeColor',cols.h,'Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(paramSetsS_tot.rhob,'FaceColor','none','EdgeColor',cols.s,'Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(paramSetsCt_tot.rhob,'FaceColor','none','EdgeColor',cols.c,'Normalization','probability','BinWidth',BinWidth);  %Ct
histogram(paramSetsNg_tot.rhob,'FaceColor','none','EdgeColor',cols.g,'Normalization','probability','BinWidth',BinWidth); %Ng

histogram(paramSetsHIV_tot.alpha_prev(paramSetsHIV_tot.p==0,:),'FaceColor',cols.h,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
histogram(paramSetsS_tot.alpha,'FaceColor',cols.s,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(paramSetsCt_tot.alpha,'FaceColor',cols.c,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth);  %Ct
histogram(paramSetsNg_tot.alpha,'FaceColor',cols.g,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng

legend('$\rho_{h}^0$',...
    '$\rho_{s}^0$',...
    '$\rho_{c}^0$',...
    '$\rho_{g}^0$',...
    '$\rho_{h,p=0}^\prime$',...
    '$\rho_{s}^\prime$',...
    '$\rho_{c}^\prime$',...
    '$\rho_{g}^\prime$',...
    'Interpreter','latex','FontSize',12,'Box','off')
xlim([0,max(paramSetsHIV_tot.alpha_prev)])
xlabel('$\rho$','Interpreter','latex','FontSize',fontSize)
ylabel('Frequency','Interpreter','latex','FontSize',fontSize)
yticks([0 0.1 0.2 0.3 0.4])
lgd = legend('boxoff','FontSize',fontSize);
lgd.NumColumns = 2;
x1=22 ; x2=18;
lgd.ItemTokenSize = [x1,x2];
lgd.FontSize = fontSize;
set(fig,'position',[100,100,480,350])
%thisFigPath = [fig_path_sm,'hist_alpha0_1dis-mod_comp_rhob.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)


%% 1.b Histogram des 4 1/alpha (sans la prep) avec comparaison des 1/rhob 
%pas dans le SM
fig = figure(3);
BinWidth=0.1;
%adding rhob
histogram(1./paramSetsHIV_tot.rhob(paramSetsHIV_tot.p==0,:),'FaceColor','none','EdgeColor',cols.h,'Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(1./paramSetsS_tot.rhob,'FaceColor','none','EdgeColor',cols.s,'Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(1./paramSetsCt_tot.rhob,'FaceColor','none','EdgeColor',cols.c,'Normalization','probability','BinWidth',BinWidth); %Ct
histogram(1./paramSetsNg_tot.rhob,'FaceColor','none','EdgeColor',cols.g,'Normalization','probability','BinWidth',BinWidth); %Ng

histogram(1./paramSetsHIV_tot.alpha_prev(paramSetsHIV_tot.p==0,:),'FaceColor',cols.h,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
histogram(1./paramSetsS_tot.alpha,'FaceColor',cols.s,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(1./paramSetsCt_tot.alpha,'FaceColor',cols.c,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(1./paramSetsNg_tot.alpha,'FaceColor',cols.g,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng

legend('$1/\rho_{h}^0$',...
    '$1/\rho_{s}^0$',...
    '$1/\rho_{c}^0$',...
    '$1/\rho_{g}^0$',...
    '$1/\rho_{h,p=0}^\prime$',...
    '$1/\rho_{s}^\prime$',...
    '$1/\rho_{c}^\prime$',...
    '$1/\rho_{g}^\prime$',...
    'Interpreter','latex','FontSize',12,'Box','off')
xlim([0,10])
xlabel('year','Interpreter','latex')


%% 2. Histogram des alpha_HIV (en faisant varier p)
vecP = sort(unique(paramSetsHIV_tot.p));
BinWidth=0.02;

fig=figure(4);
i=4;
histogram(paramSetsHIV_tot.rhob(paramSetsHIV_tot.p==0,:),'FaceColor','none','EdgeColor',cols.h,'Normalization','probability','BinWidth',BinWidth); %HIV
var{1} = ['$\rho_{h}^0$'];
hold on
var{2} = [' '];
histogram(0.1:0.01:1,'EdgeColor','none','Normalization','probability','FaceColor','none');
var{3} = [' '];
histogram(0.1:0.01:1,'EdgeColor','none','Normalization','probability','FaceColor','none');

histogram(paramSetsHIV_tot.alpha_prep(paramSetsHIV_tot.p==0,:),'FaceColor',cols.h,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
var{4} = ['$\rho_{h,p=',num2str(0),'}^\prime$'];
hold on
for p=vecP(2:end)'
    i=i+1;
    histogram(paramSetsHIV_tot.alpha_prep(paramSetsHIV_tot.p==p,:),'FaceColor',colMat(i-4,:),'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
    var{i} = ['$\rho_{h,p=',num2str(p),'}^\prime$'];
end
legend(var,'Interpreter','latex','FontSize',fontSize,'Box','off')
xlabel('$\rho$','Interpreter','latex','FontSize',fontSize);
xlim([0,1.05*max(paramSetsHIV_tot.alpha_prep(paramSetsHIV_tot.p==0,:))])
ylim([0,0.4])
lgd = legend('boxoff','FontSize',fontSize);
lgd.NumColumns = 3;
x1=22 ; x2=18;
lgd.ItemTokenSize = [x1,x2];
lgd.FontSize = fontSize;
ylabel('Frequency','Interpreter','latex','FontSize',fontSize)
yticks([0 0.1 0.2 0.3 0.4])
set(fig,'position',[100,100,480,350])
thisFigPath = [fig_path_sm,'hist_alpha_prep_hiv.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

%% 3. Histogram des 4 Pi(0)^all sans le depistage volontaire (nombre de personnes infectees)
close all
fig = figure(5);
BinWidth=0.01;
histogram(paramSetsHIV_tot.Ptot_prev_base(paramSetsHIV_tot.p==0,:),'FaceColor','none','EdgeColor',cols.h,'Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(paramSetsS_tot.P_base,'FaceColor','none','EdgeColor',cols.s,'Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(paramSetsCt_tot.P_base,'FaceColor','none','EdgeColor',cols.c,'Normalization','probability','BinWidth',BinWidth); %Ct
histogram(paramSetsNg_tot.P_base,'FaceColor','none','EdgeColor',cols.g,'Normalization','probability','BinWidth',BinWidth); %Ng

histogram(paramSetsHIV_tot.Ptot_prev_0(paramSetsHIV_tot.p==0,:),'FaceColor',cols.h,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
histogram(paramSetsS_tot.P0,'FaceColor',cols.s,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(paramSetsCt_tot.P0,'FaceColor',cols.c,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(paramSetsNg_tot.P0,'FaceColor',cols.g,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng

lgd = legend('$\Pi_{h,p=0}^{\rm all}(\rho_h^0)$',...
    '$\Pi_{s}^{\rm all}(\rho_s^0)$',...
    '$\Pi_{c}^{\rm all}(\rho_c^0)$',...
    '$\Pi_{g}^{\rm all}(\rho_g^0)$',...
    '$\Pi_{h,p=0}^{\rm all}(0)$',...
    '$\Pi_{s}^{\rm all}(0)$',...
    '$\Pi_{c}^{\rm all}(0)$',...
    '$\Pi_{g}^{\rm all}(0)$',...
    'Interpreter','latex','FontSize',0.95*fontSize,'Box','off','Location','NorthEast')
x1=22 ; x2=18;
lgd.ItemTokenSize = [x1,x2];
lgd.NumColumns = 2;

xlim([0,0.9])
ylim([0,0.65])
ylabel('Frequency','Interpreter','latex','FontSize',fontSize)
xlabel('Prevalence','Interpreter','latex','FontSize',fontSize)
yticks([0 0.2 0.4 0.6])
xticks([0 0.2 0.4 0.6 0.8 1])
set(fig,'position',[100,100,480,350])

%thisFigPath = [fig_path_sm,'hist_prevalence0_1dis-mod.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

%% 4. Histogram des Ph(0)
close all
vecP = sort(unique(paramSetsHIV_tot.p));
BinWidth=0.5/100;

fig=figure(6);
i=7;
histogram(paramSetsHIV_tot.Pun_prep_base(paramSetsHIV_tot.p==0,:),'FaceColor','none','EdgeColor',cols.h,'Normalization','probability','BinWidth',BinWidth); %HIV
var{1} = ['$\Pi_{h,p=0}^{\rm und}(\rho^0_h)$'];
hold on
histogram((0.1:0.01:1)/100,'EdgeColor','none','Normalization','probability','FaceColor','none');
var{2} = [' '];
histogram((0.1:0.01:1)/100,'EdgeColor','none','Normalization','probability','FaceColor','none');
var{3} = [' '];
histogram((0.1:0.01:1)/100,'EdgeColor','none','Normalization','probability','FaceColor','none');
var{4} = [' '];
histogram((0.1:0.01:1)/100,'EdgeColor','none','Normalization','probability','FaceColor','none');
var{5} = [' '];
histogram((0.1:0.01:1)/100,'EdgeColor','none','Normalization','probability','FaceColor','none');
var{6} = [' '];
histogram(paramSetsHIV_tot.Pun_prep_0(paramSetsHIV_tot.p==0,:),'FaceColor',cols.h,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
var{7} = ['$\Pi_{h,p=0}^{\rm und}(0)$'];
for p=vecP(2:end)'
    i=i+1;
    histogram(paramSetsHIV_tot.Pun_prep_0(paramSetsHIV_tot.p==p,:),'FaceColor',colMat(i-7,:),'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
    var{i} = ['$\Pi_{h,p=',num2str(p),'}^{\rm und}(0)$'];
end
xlabel('$\Pi_{h}(0), \%$','Interpreter','latex')
xlim([0,16/100])
lgd = legend(var,'Interpreter','latex','FontSize',0.95*fontSize,'Box','off');
lgd.NumColumns = 2;
xticks([0 0.05 0.1 0.15])
x1=22 ; x2=18;
lgd.ItemTokenSize = [x1,x2];
lgd.FontSize = fontSize;
xlabel('Prevalence','FontSize',fontSize,'Interpreter','latex')
ylabel('Frequency','FontSize',fontSize,'Interpreter','latex')
set(fig,'position',[100,100,480,350])

thisFigPath = [fig_path_sm,'hist_P_und_hiv.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

%% 6. Impact of P_prev_base on P_und_prep_base
close all
fontSize = 18;
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\INITIALISATION\parametersSets')
nbSim = 500;
b = 100;
vecvalPHIV = [0.143;5/100;10/100;15/100;20/100;25/100;30/100];
BinWidth   = 0.0005*50/100;
colFace = [1,1,1;colMat];
colEdge = {rougeHIV,'none','none','none','none','none','none'};
fig = figure(7);
k=1;
for valPHIV=vecvalPHIV'
    Pund_base=[];
    for i = 1:nbSim
        if valPHIV==0.143
        	PHIV_prevagay_VT = randPERT(12.0,14.3,16.9,1)/100;
        else
            PHIV_prevagay_VT = valPHIV;    
        end
        mu              = 1/randPERT(27.2,30.6,33.7,1); 
        R_sict_VT       = 1/(1-PHIV_prevagay_VT);
        sigma_sict      = 52/randPERT(6.7,8.2,9.8,1);
        rhobh           = 1/randPERT(1.6,1.9,2.2,1);
        theta_sict0     = 1/randPERT(4,4.4,10,1);
        theta_sict      = theta_sict0+rhobh;
        gamma_sict0     = 0;
        gamma_sict      = rhobh+gamma_sict0;
        ratioBeta       = randPERT(8.4,9.1,9.6,1);
        betaC_sict      = R_sict_VT*(sigma_sict+gamma_sict+mu)*(theta_sict+mu)/(ratioBeta*(theta_sict+mu)+sigma_sict); %check
        betaI_sict      = ratioBeta*betaC_sict; %check
        [R0_sict,~,alpha_sict] = Rp_SICR_v4(betaI_sict,betaC_sict,theta_sict0,sigma_sict,gamma_sict0,mu,b,0); %check
        %P0_sict = 1-1/R0_sict;
        %undiagnosed prevalence :
        Pund_base(i) = mu*(theta_sict+sigma_sict+mu)/(betaI_sict*(theta_sict+mu)+betaC_sict*sigma_sict)*(R_sict_VT-1);
        %Pund_0    = mu*(theta_sict0+sigma_sict+mu)/(betaI_sict*(theta_sict0+mu)+betaC_sict*sigma_sict)*(R0_sict-1);
    end
    histogram(Pund_base,'EdgeColor',colEdge{k},'FaceColor',colFace(k,:),'Normalization','probability','BinWidth',BinWidth); %HIV
    hold on
    k=k+1;
end
%
xlim([0,0.02])
lgd = legend(cellstr(['ref.';num2str(vecvalPHIV(2:end))]),'Interpreter','latex','FontSize',fontSize,'Box','off')
lgd.Title.String = '$\Pi^{\rm all}_h(\rho^0_h)$';
xlabel('$\Pi_{h,p=0}^{\rm und}(\rho^0_h)$','Interpreter','latex','FontSize',fontSize)
ylabel('Frequency','Interpreter','latex','FontSize',fontSize)
xticks([0 0.005 0.01 0.015 0.02])
yticks([0 0.1 0.2 0.3 0.4])
x1=22 ; x2=18;
lgd.ItemTokenSize = [x1,x2];
lgd.NumColumns = 2;

set(fig,'position',[100,100,480,350])

thisFigPath = [fig_path_sm,'hist_Ph_und_base_acc_to_Ph_base.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

%% 7. Impact of P_prev_base on P_und_prep_base
close all
nbSim = 2000;
b=100;
vecRhobh = [1/1.9,1./[1,2,3,4,5]];
BinWidth = 0.001*50/100;
fig = figure(8);
k=1;
colFace = [1,1,1;colMat];
colEdge = {rougeHIV,'none','none','none','none','none','none'};
for valRhobh=vecRhobh
    Pund_base=[];
    for i = 1:nbSim
        mu              = 1/randPERT(27.2,30.6,33.7,1); 
        PHIV_prevagay_VT= randPERT(12.0,14.3,16.9,1)/100;
        R_sict_VT       = 1/(1-PHIV_prevagay_VT);
        sigma_sict      = 52/randPERT(6.7,8.2,9.8,1);
        if (valRhobh==1/1.9)
            rhobh           =  1/randPERT(1.6,1.9,2.2,1);
        else
            rhobh           = valRhobh;
        end
        theta_sict0     = 1/randPERT(4,4.4,10,1);
        theta_sict      = theta_sict0+rhobh;
        gamma_sict0     = 0;
        gamma_sict      = rhobh+gamma_sict0;
        ratioBeta       = randPERT(8.4,9.1,9.6,1);
        betaC_sict      = R_sict_VT*(sigma_sict+gamma_sict+mu)*(theta_sict+mu)/(ratioBeta*(theta_sict+mu)+sigma_sict); %check
        betaI_sict      = ratioBeta*betaC_sict; %check
        [R0_sict,~,alpha_sict] = Rp_SICR_v4(betaI_sict,betaC_sict,theta_sict0,sigma_sict,gamma_sict0,mu,b,0); %check
        %P0_sict = 1-1/R0_sict;
        %undiagnosed prevalence :
        Pund_base(i) = mu*(theta_sict+sigma_sict+mu)/(betaI_sict*(theta_sict+mu)+betaC_sict*sigma_sict)*(R_sict_VT-1);
        %Pund_0    = mu*(theta_sict0+sigma_sict+mu)/(betaI_sict*(theta_sict0+mu)+betaC_sict*sigma_sict)*(R0_sict-1);
    end
    histogram(Pund_base,'EdgeColor',colEdge{k},'FaceColor',colFace(k,:),'Normalization','probability','BinWidth',BinWidth); %HIV
    hold on
    k=k+1;
end
%
xlim([0,0.018])
lgd = legend({'ref.','1','1/2','1/3','1/4','1/5'},'Interpreter','latex','FontSize',fontSize,'Box','off')
lgd.Title.String = '$\rho^0_h$';
xlabel('$\Pi_{h,p=0}^{\rm und}(\rho^0_h)$','Interpreter','latex','FontSize',fontSize)
ylabel('Frequency','Interpreter','latex','FontSize',fontSize)
xticks([0 0.005 0.01 0.015 0.02])
x1=22 ; x2=18;
lgd.ItemTokenSize = [x1,x2];
set(fig,'position',[100,100,480,350])

thisFigPath = [fig_path_sm,'hist_Ph_und_base_acc_to_rho_base.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)