clear all; close all;
parametrizationNo = 4;
vecRounds         = 1:160;
rougeHIV = [215, 0, 0]/255;     cols.h = rougeHIV;
jauneS   = [250, 215, 0]/255;   cols.s = jauneS;
bleuCt   = [56, 57, 186]/255;   cols.c = bleuCt;
vertNg   = [43, 152, 38]/255;   cols.g = vertNg;

%/ParameterAnalysis/createParameterSets_alea.m

backupPath = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis/results_',num2str(parametrizationNo),'/'];

%Reading parameterSets
setsPath = ['C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/ParameterAnalysis/paramSets_',num2str(parametrizationNo)];

%
paramSetsCt_tot=[];paramSetsNg_tot=[];paramSetsHIV_tot=[];paramSetsS_tot=[];
for roundNo = vecRounds
    pathFile = [setsPath,'/round_',num2str(roundNo)];
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
saveas(fig,[backupPath,'hist_alpha0_1dis-mod.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_alpha0_1dis-mod.pdf'])
saveas(fig,['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\parameterAnalysis\hist_alpha0_1dis-mod.png'])


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
    '$\rho_{h,(p=0)}^\prime$',...
    '$\rho_{s}^\prime$',...
    '$\rho_{c}^\prime$',...
    '$\rho_{g}^\prime$',...
    'Interpreter','latex','FontSize',12,'Box','off')
xlim([0,1])
xlabel('$1/year$','Interpreter','latex')

saveas(fig,[backupPath,'hist_alpha0_1dis-mod_comp_rhob.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_alpha0_1dis-mod_comp_rhob.pdf'])
saveas(fig,['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\parameterAnalysis\hist_alpha0_1dis-mod_comp_rhob.png'])


%% 1.b Histogram des 4 1/alpha (sans la prep) avec comparaison des 1/rhob
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
saveas(fig,[backupPath,'hist_duration_alpha0_1dis-mod_comp_rhob.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_duration_alpha0_1dis-mod_comp_rhob.pdf'])
saveas(fig,['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\parameterAnalysis\hist_duration_alpha0_1dis-mod_comp_rhob.png'])


%% 2. Histogram des alpha_HIV (en faisant varier p)
vecP = sort(unique(paramSetsHIV_tot.p));
BinWidth=0.02;

fig=figure(4);
i=2;
histogram(paramSetsHIV_tot.rhob(paramSetsHIV_tot.p==0,:),'FaceColor','none','EdgeColor',cols.h,'Normalization','probability','BinWidth',BinWidth); %HIV
var{1} = ['$\rho_{h}^0$'];
hold on
histogram(paramSetsHIV_tot.alpha_prep(paramSetsHIV_tot.p==0,:),'FaceColor',cols.h,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
var{2} = ['$\rho_{h,(p=',num2str(0),')}^\prime$'];
hold on
for p=vecP(2:end)'
    i=i+1;
    histogram(paramSetsHIV_tot.alpha_prep(paramSetsHIV_tot.p==p,:),'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
    var{i} = ['$\rho_{h,(p=',num2str(p),')}^\prime$'];
end
legend(var,'Interpreter','latex','FontSize',12,'Box','off')
xlabel('1/year','Interpreter','latex')
saveas(fig,[backupPath,'hist_alpha_prep_hiv.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_alpha_prep_hiv.pdf'])
saveas(fig,['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\parameterAnalysis\hist_alpha_prep_hiv.png'])

%% 3. Histogram des 4 Pi(0) sans le depistage volontaire (nombre de personnes infectees)
fig = figure(5);
BinWidth=0.01;
histogram(paramSetsHIV_tot.Ptot_prev_0(paramSetsHIV_tot.p==0,:),'FaceColor',cols.h,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(paramSetsS_tot.P0,'FaceColor',cols.s,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(paramSetsCt_tot.P0,'FaceColor',cols.c,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(paramSetsNg_tot.P0,'FaceColor',cols.g,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
legend('$\Pi_{h,p=0}(0)$',...
    '$\Pi_{s}(0)$',...
    '$\Pi_{c}(0)$',...
    '$\Pi_{g}(0)$',...
    'Interpreter','latex','FontSize',12,'Box','off')
xlim([0,1])
saveas(fig,[backupPath,'hist_prevalence0_1dis-mod.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_prevalence0_1dis-mod.pdf'])
saveas(fig,['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\parameterAnalysis\hist_prevalence0_1dis-mod.png'])


%% 4. Histogram des Ph(0)
vecP = sort(unique(paramSetsHIV_tot.p));
BinWidth=0.5;

fig=figure(6);
i=2;
histogram(paramSetsHIV_tot.Pun_prep_base(paramSetsHIV_tot.p==0,:)*100,'FaceColor','none','EdgeColor',cols.h,'Normalization','probability','BinWidth',BinWidth); %HIV
var{1} = ['$\Pi_{h,p=0}^{und}(\rho^0_h)$'];
hold on
histogram(paramSetsHIV_tot.Pun_prep_0(paramSetsHIV_tot.p==0,:)*100,'FaceColor',cols.h,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
var{2} = ['$\Pi_{h,p=0}^{und}(0)$'];
hold on
for p=vecP(2:end)'
    i=i+1;
    histogram(paramSetsHIV_tot.Pun_prep_0(paramSetsHIV_tot.p==p,:)*100,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
    var{i} = ['$\Pi_{h,p=',num2str(p),'}^{und}(0)$'];
end
xlabel('$\Pi_{h}(0), \%$','Interpreter','latex')
xlim([0,16])
legend(var,'Interpreter','latex','FontSize',12,'Box','off')
saveas(fig,[backupPath,'hist_P_und_hiv.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_P_und_hiv.pdf'])
saveas(fig,['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\parameterAnalysis\hist_P_und_hiv.png'])


%% 6. Impact of P_prev_base on P_und_prep_base
close all
nbSim = 1500;
b = 100;
vecvalPHIV = [5/100,10/100,15/100,20/100,25/100,30/100,35/100];
BinWidth   = 0.0005*50;
fig = figure(7);
for valPHIV=vecvalPHIV
    Pund_base=[];
    for i = 1:nbSim
        PHIV_prevagay_VT = valPHIV;    %randPERT(12.0,14.3,16.9,1)/100;
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
    histogram(Pund_base*100,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
    hold on
end
xlim([0,0.025*100])
legend({'$\Pi(\rho^0_h)=5\%$','$\Pi(\rho^0_h)=10\%$','$\Pi(\rho^0_h)=15\%$','$\Pi(\rho^0_h)=20\%$',...
    '$\Pi(\rho^0_h)=25\%$','$\Pi(\rho^0_h)=30\%$','$\Pi(\rho^0_h)=35\%$' },'Interpreter','latex','FontSize',12)
xlabel('$\Pi_{h,p=0}^{und}(\rho^0_h), \%$','Interpreter','latex')
saveas(fig,[backupPath,'hist_Ph_und_base_acc_to_Ph_base.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_Ph_und_base_acc_to_Ph_base.pdf'])
saveas(fig,['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\parameterAnalysis\hist_Ph_und_base_acc_to_Ph_base.png'])


%% 7. Impact of P_prev_base on P_und_prep_base
close all
nbSim = 1000;
b=100;
vecRhobh = 1./[1,2,3,4,5];
BinWidth = 0.001*50;
fig = figure(8);
for valRhobh=vecRhobh
    Pund_base=[];
    for i = 1:nbSim
        mu              = 1/randPERT(27.2,30.6,33.7,1); 
        PHIV_prevagay_VT= randPERT(12.0,14.3,16.9,1)/100;
        R_sict_VT       = 1/(1-PHIV_prevagay_VT);
        sigma_sict      = 52/randPERT(6.7,8.2,9.8,1);
        rhobh           = valRhobh;         %1/randPERT(1.6,1.9,2.2,1);
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
    histogram(Pund_base*100,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
    hold on
end
xlim([0,0.02*100])
legend({'$\rho^0_h=1$','$\rho^0_h=1/2$','$\rho^0_h=1/3$',...
    '$\rho^0_h=1/4$','$\rho^0_h=1/5$'},'Interpreter','latex','FontSize',12)
xlabel('$\Pi_{h,p=0}^{und}(\rho^0_h), \%$','Interpreter','latex')
saveas(fig,[backupPath,'hist_Ph_und_base_acc_to_rho_base.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_Ph_und_base_acc_to_rho_base.pdf'])
saveas(fig,['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\parameterAnalysis\hist_Ph_und_base_acc_to_rho_base.png'])
%%prevalence sans les traités ? utiliser la fonction scitp qui rend les
%%prevalences, en mettant p a 0
[Rp,Lambdap,alpha,Ptot,Pun] = Rp_SICTP(betaI_sict,betaC_sict,theta_sict0,gamma_sict0,sigma_sict,0,0,0,mu,b,rhobh);


