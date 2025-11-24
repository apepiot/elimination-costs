
% difference with v_4 : the algorithm to find alpha has changed
clear all; close all;
parametrizationNo = 1;
vecRounds = 1;

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

list_kits_to_consider = list_kits([1,5,6,7,11,12,13,15]);
createParamRho;
paramRho.rho_h = 0;
paramRho.rho_s = 0;
paramRho.rho_c = 0;
paramRho.rho_g = 0;
paramRho.VTunderART = 1;

for roundNo=1
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
        
        for p=[0,0.1,0.2,0.3]
            disp(['ID_ech=',num2str(id_ech),', p=',num2str(p),', nEch=',num2str(paramCt_ech.nbEch)]);
            
            if  (1)
            paramTab{1} = table2struct(paramCt_ech);
            paramTab{2} = table2struct(paramNg_ech);
            paramTab{3} = table2struct(paramHIV_ech(paramHIV_ech.p==p,:));
            paramTab{4} = table2struct(paramS_ech);
            nbEch = paramTab{1}.nbEch;
            numSim = numSim+1;
            
            
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
            
            tStart = tic;
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
            
            tabAlpha(numSim,:).timeCompil = toc(tStart);
            tabAlpha(numSim,:).msg = msg;
            
            tabAlpha
            end
        end
    end
end




