
% difference with v_4 : the algorithm to find alpha has changed
clear all; close all;
parametrizationNo = 3;
vecRounds = 1:240;  %[50:61,40:44];%[15,17];%[14,15,17,18];%14:15;

%/ParameterAnalysis/createParameterSets.m

backupPath = ['C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/ParameterAnalysis/results_',num2str(parametrizationNo),'/'];

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


%% 1.a Histogram des 4 alpha (sans la prep)
fig = figure(1);
BinWidth=0.01;
histogram(paramSetsHIV_tot.alpha_prev(paramSetsHIV_tot.p==0,:),'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(paramSetsS_tot.alpha,'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(paramSetsCt_tot.alpha,'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(paramSetsNg_tot.alpha,'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
legend('$\rho_{h}^\prime(p=0)$',...
    '$\rho_{s}^\prime$',...
    '$\rho_{c}^\prime$',...
    '$\rho_{g}^\prime$',...
    'Interpreter','latex','FontSize',12,'Box','off')
xlim([0,1])
xlabel('1/year')
saveas(fig,[backupPath,'hist_alpha0_1dis-mod.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_alpha0_1dis-mod.pdf'])

%% 1.b Histogram des 4 alpha (sans la prep) avec comparaison des rhob
fig = figure(2);
BinWidth=0.01;
%adding rhob
histogram(paramSetsHIV_tot.rhob(paramSetsHIV_tot.p==0,:),'FaceColor','none','EdgeColor','red','Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(paramSetsS_tot.rhob,'FaceColor','none','EdgeColor','yellow','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(paramSetsCt_tot.rhob,'FaceColor','none','EdgeColor','blue','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(paramSetsNg_tot.rhob,'FaceColor','none','EdgeColor','green','Normalization','probability','BinWidth',BinWidth); %Ng

histogram(paramSetsHIV_tot.alpha_prev(paramSetsHIV_tot.p==0,:),'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
histogram(paramSetsS_tot.alpha,'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(paramSetsCt_tot.alpha,'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(paramSetsNg_tot.alpha,'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng

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
xlabel('1/year')

saveas(fig,[backupPath,'hist_alpha0_1dis-mod_comp_rhob.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_alpha0_1dis-mod_comp_rhob.pdf'])

%% 1.b Histogram des 4 1/alpha (sans la prep) avec comparaison des 1/rhob
fig = figure(3);
BinWidth=0.1;
%adding rhob
histogram(1./paramSetsHIV_tot.rhob(paramSetsHIV_tot.p==0,:),'FaceColor','none','EdgeColor','red','Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(1./paramSetsS_tot.rhob,'FaceColor','none','EdgeColor','yellow','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(1./paramSetsCt_tot.rhob,'FaceColor','none','EdgeColor','blue','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(1./paramSetsNg_tot.rhob,'FaceColor','none','EdgeColor','green','Normalization','probability','BinWidth',BinWidth); %Ng

histogram(1./paramSetsHIV_tot.alpha_prev(paramSetsHIV_tot.p==0,:),'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
histogram(1./paramSetsS_tot.alpha,'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(1./paramSetsCt_tot.alpha,'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(1./paramSetsNg_tot.alpha,'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng

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
xlabel('year')
saveas(fig,[backupPath,'hist_duration_alpha0_1dis-mod_comp_rhob.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_duration_alpha0_1dis-mod_comp_rhob.pdf'])


%% 2. Histogram des alpha_HIV (en faisant varier p)
vecP = sort(unique(paramSetsHIV_tot.p));
BinWidth=0.02;

fig=figure(4);
i=2;
histogram(paramSetsHIV_tot.rhob(paramSetsHIV_tot.p==0,:),'FaceColor','none','EdgeColor','red','Normalization','probability','BinWidth',BinWidth); %HIV
var{1} = ['$\rho_{h}^0$'];
hold on
histogram(paramSetsHIV_tot.alpha_prep(paramSetsHIV_tot.p==0,:),'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
var{2} = ['$\rho_{h,(p=',num2str(0),')}^\prime$'];
hold on
for p=vecP(2:end)'
    i=i+1;
    histogram(paramSetsHIV_tot.alpha_prep(paramSetsHIV_tot.p==p,:),'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
    var{i} = ['$\rho_{h,(p=',num2str(p),')}^\prime$'];
end
legend(var,'Interpreter','latex','FontSize',12,'Box','off')
xlabel('1/year')
saveas(fig,[backupPath,'hist_alpha_prep_hiv.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_alpha_prep_hiv.pdf'])

%% 3. Histogram des 4 Pi(0) sans le depistage volontaire (nombre de personnes infectees)
fig = figure(5);
BinWidth=0.01;
histogram(paramSetsHIV_tot.Ptot_prev_0(paramSetsHIV_tot.p==0,:),'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(paramSetsS_tot.P0,'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(paramSetsCt_tot.P0,'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(paramSetsNg_tot.P0,'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
legend('$\Pi_{h,p=0}(0)$',...
    '$\Pi_{s}(0)$',...
    '$\Pi_{c}(0)$',...
    '$\Pi_{g}(0)$',...
    'Interpreter','latex','FontSize',12,'Box','off')
xlim([0,1])
saveas(fig,[backupPath,'hist_prevalence0_1dis-mod.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_prevalence0_1dis-mod.pdf'])

%% 4. Histogram des Ph(0)
vecP = sort(unique(paramSetsHIV_tot.p));
BinWidth=0.5;

fig=figure(6);
i=2;
histogram(paramSetsHIV_tot.Pun_prep_base(paramSetsHIV_tot.p==0,:)*100,'FaceColor','none','EdgeColor','red','Normalization','probability','BinWidth',BinWidth); %HIV
var{1} = ['$\Pi_{h,p=0}^{und}(\rho^0_h)$'];
hold on
histogram(paramSetsHIV_tot.Pun_prep_0(paramSetsHIV_tot.p==0,:)*100,'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
var{2} = ['$\Pi_{h,p=0}^{und}(0)$'];
hold on
for p=vecP(2:end)'
    i=i+1;
    histogram(paramSetsHIV_tot.Pun_prep_0(paramSetsHIV_tot.p==p,:)*100,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
    var{i} = ['$\Pi_{h,p=',num2str(p),'}^{und}(0)$'];
end
xlabel('\Pi_{h}(0), %')
xlim([0,16])
legend(var,'Interpreter','latex','FontSize',12,'Box','off')
saveas(fig,[backupPath,'hist_P_und_hiv.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_P_und_hiv.pdf'])


%% 5. Trouver alpha avec le modèle à 4 infections
ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
log_path = ['C:/Users/Moi/Documents/Projets/Resultats_knitro/parametrization_',num2str(parametrizationNo),'/'];
backupPath = ['C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/ParameterAnalysis/results_',num2str(parametrizationNo),'/'];
paramSolverAlpha.tolP0  = 0.5e-4;
paramSolverAlpha.maxBndAlpha = 50;
paramSolverAlpha.nbRelanceMax = 2;
paramSolverAlpha.timeLimit = 30;
paramSolverAlpha.method_alpha = 'dicho';
paramSolverAlpha.tolAlpha = 1e-3;
paramSolverAlpha.iterMaxDicho = 20;


f=1;
%sauvegarder les parametres ci dessus dans un txt.

afficherOutput=0;

addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/';
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/Initialisation';
addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models'
setupOnce;

%attention a ne pas ecraser le fichier tabAlpha

list_kits = {{'HIV'},{'syphilis'},{'Ct'},{'Ng'},...
    {'HIV','syphilis'},{'HIV','Ct'},{'HIV','Ng'}, {'syphilis','Ct'},...
    {'syphilis','Ng'},{'Ct','Ng'}, {'HIV','syphilis','Ct'},{'HIV','syphilis','Ng'},...
    {'HIV','Ct','Ng'},{'syphilis','Ct','Ng'}, {'HIV','syphilis','Ct','Ng'}};

%list_kits_to_consider = list_kits([1,5,6,7,11,12,13,15]); 
list_kits_to_consider = list_kits(15);%
createParamRho;
paramRho.rho_h = 0;
paramRho.rho_s = 0;
paramRho.rho_c = 0;
paramRho.rho_g = 0;
paramRho.VTunderART = 1;

for roundNo=vecRounds
    log_path_round = [log_path, '_round_', num2str(roundNo),'/'];
    mkdir([backupPath,'_round_',num2str(roundNo),'/'])
    writetable(struct2table(paramRho),[backupPath,'_round_',num2str(roundNo),'/paramRho.txt'],'WriteVariableNames',true)
    writetable(struct2table(paramSolverAlpha),[backupPath,'_round_',num2str(roundNo),'/paramSolverAlpha.txt'],'WriteVariableNames',true)

    paramCt_round  = paramSetsCt_tot(paramSetsCt_tot.roundNo==roundNo,:);
    paramNg_round  = paramSetsNg_tot(paramSetsNg_tot.roundNo==roundNo,:);
    paramHIV_round = paramSetsHIV_tot(paramSetsHIV_tot.roundNo==roundNo,:);
    paramS_round   = paramSetsS_tot(paramSetsS_tot.roundNo==roundNo,:);
    numSim=0;
    elimCostbyModel = createelimCostbyModel(list_kits_to_consider);
    tabAlpha        = createTabAlpha(list_kits_to_consider,0);
    
    vecID_ech = unique(paramCt_round.IDech_id);
    for id_ech=vecID_ech'
        paramCt_ech  = paramCt_round(paramCt_round.IDech_id==id_ech,:);
        paramNg_ech  = paramNg_round(paramNg_round.IDech_id==id_ech,:);
        paramS_ech   = paramS_round(paramS_round.IDech_id==id_ech,:);
        paramHIV_ech = paramHIV_round(paramHIV_round.IDech_id==id_ech,:);
        vecP = paramHIV_ech.p;
        b  = paramS_ech.pi;
        mu = paramS_ech.mu;
        
        for p=vecP'
            disp(['ID_ech=',num2str(id_ech),', p=',num2str(p),', nEch=',num2str(paramCt_ech.nbEch)]);
            
            if  (1)%id_ech>=11251838530003 && p>0.2) || (id_ech>11251838530003)
            paramTab{1} = table2struct(paramCt_ech);
            paramTab{2} = table2struct(paramNg_ech);
            paramTab{3} = table2struct(paramHIV_ech(paramHIV_ech.p==p,:));
            paramTab{4} = table2struct(paramS_ech);
            nbEch = paramTab{1}.nbEch;
            numSim = numSim+1;
            
            % sauvegarder paramRho dans une table aussi
            % et modifier le code MAIN pour lire le paramRho
            % ajouter un time compil a alpha
            
            clear alphas;
            tabAlpha = [tabAlpha;array2table(zeros(1,size(tabAlpha,2)),'VariableNames',tabAlpha.Properties.VariableNames)];
            tabAlpha(numSim,:).IDech = id_ech;
            tabAlpha(numSim,:).roundNo = roundNo;
            tabAlpha(numSim,:).nbEch = nbEch;
            tabAlpha(numSim,:).p = p;
            msg = [];
            
            %check Feasability
            [infElim,msgSolver,msgSol] = checkFeasability_v2(paramTab,paramRho,b,mu,f,'hscg',afficherOutput,paramSolverAlpha,log_path_round,ampl_models_dir);
            
            tabAlpha(numSim,["elim_h","elim_s","elim_c","elim_g"]) = struct2table(infElim);
            
            tic;
            for numKit=1:length(list_kits_to_consider)
                kit = list_kits_to_consider{numKit}; k = indexKit(kit);
                if ~contains(k,'h')
                    mods = {[k,'_h',k],[k,'_',k]};
                else
                    mods = {[k,'_',k]};
                end
                
                for mod=mods
                    m=mod{:};
                    doNotRunModel=0;
                    for inf='hscg'
                        if contains(m,inf) && infElim.(inf)
                            doNotRunModel = 1;
                        end
                    end
          
                    if ~doNotRunModel
                        [alpha,P,~,msgAlpha_mod,elim_i] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit,m,afficherOutput,...
                            log_path_round,paramSolverAlpha,ampl_models_dir);
                    else
                        alpha=0;
                        msgAlpha_mod = '1';
                        elim_i = '0';
                    end
                    %alphas.(m) = alpha;   msgAlphas.(m) = msgAlpha_mod;
                    tabAlpha(numSim,[m,'_elim'])= table({elim_i}, 'VariableNames',{[m,'_elim']});
                    elimCostbyModel(numKit,:).alphaK    = alpha;
                    elimCostbyModel(numKit,:).msgAlphaK = msgAlpha_mod;
                    tabAlpha(numSim,:).(m) = alpha;
                    msg = [msg, ' ',msgAlpha_mod];
                end
            end
            
            tabAlpha(numSim,:).timeCompil = toc;
            tabAlpha(numSim,:).msg = msg;
            
            tabAlpha
            end
            writetable(tabAlpha,[backupPath,'_round_',num2str(roundNo),'/tabAlpha'],'WriteVariableNames',true)
        end
    end
end


%% Merging tabAlpha_i with tabAlpha_ijk...
roundNo=18;
backupPath = ['C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/ParameterAnalysis/results_1/_round_',num2str(roundNo),'/'];

opts1 = detectImportOptions([backupPath,'tabAlpha_all.txt']);
opts1.VariableTypes(32:52)={'char'}; %decalé
tab1 = readtable([backupPath,'tabAlpha_all.txt'],opts1);

opts2 = detectImportOptions([backupPath,'tabAlpha_hscg.txt']);
opts2.VariableTypes([7,12])={'char'}; %decalé
tab2 = readtable([backupPath,'tabAlpha_hscg.txt'],opts2);

tabTot = join(tab1,tab2,'Keys',{'IDech','p','nbEch','roundNo'} );
tabTot.timeCompil = tabTot.timeCompil_tab1 + tabTot.timeCompil_tab2;
tabTot.msg = strcat(tabTot.msg_tab1,strcat({' '},tabTot.msg_tab2));

tabTot.elim_h = tabTot.elim_h_tab1;
tabTot.elim_s = tabTot.elim_s_tab1;
tabTot.elim_c = tabTot.elim_c_tab1;
tabTot.elim_g = tabTot.elim_g_tab1;

writetable(tabTot,[backupPath,'tabAlpha_last.txt']);

% tab2 = readtable([backupPath,'tabAlpha_ij_v1.txt']);
% tab3 = readtable([backupPath,'tabAlpha_ijk_v1.txt']);
%
% tabTot = join(join(tab1,tab2),tab3);
%
% writetable(tabTot,[backupPath,'tabAlpha_all_v1.txt']);

%% 6. Porportion of infections eliminated without introducing VT (because of a high value of p)


