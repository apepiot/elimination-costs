% difference with MAIN_8 and MAIN_8bis: we focus only on HIV models (not strategies)

clear all;
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/';
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/Initialisation/';
setupOnce;

f=1;
vecRound = 16;
backUpFolder='Runtest_1'; %mettre un nom de dossier qui n'existe pas
afficherOutput=0;
byModels = 1;
paramSolver.tolP0 = 0.5e-4;
paramsC.sup=1;
paramsC.inf=-15;
paramsC.tolC=1e-3;
paramsC.iterMax=20;
paramSolver.verbose = afficherOutput;
paramSolver.sup_bnd_alpha = 20;
parametrizationNo = 1;
paramSolver.timeSolver = 30;
%---------------------------------------------------------------%
pathBackup = ['C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/ParameterAnalysis/results_',num2str(parametrizationNo),'/'];
pathParams = ['C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/ParameterAnalysis/paramSets_',num2str(parametrizationNo),'/'];
%%
pathKnitroLog = ['C:/Users/Moi/Documents/Projets/Resultats_knitro/knitro_out_rounds/',backUpFolder,'/'];
ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';

elimCost_roundTemplate = table('Size',[1,10],...
    'VariableTypes',{'double','double','double','double','double','double','string','double','double','double'},...
    'VariableNames',{'IDech','strat','HIV','syphilis','Ct','Ng','code_err','timeCompil','f','p'});
elimCostbyModelTemplate = table('Size',[1,10],...
    'VariableTypes',{'double','string','double','double','double','double','string','double','double','double'},...
    'VariableNames',{'IDech','kit','HIV','syphilis','Ct','Ng','msgCnnK','timeCompil','f','p'});
%%
list_kits = {{'HIV'},{'syphilis'},{'Ct'},{'Ng'},...
    {'HIV','syphilis'},{'HIV','Ct'},{'HIV','Ng'}, {'syphilis','Ct'},...
    {'syphilis','Ng'},{'Ct','Ng'}, {'HIV','syphilis','Ct'},{'HIV','syphilis','Ng'},...
    {'HIV','Ct','Ng'},{'syphilis','Ct','Ng'}, {'HIV','syphilis','Ct','Ng'}};

list_kits_to_consider = list_kits(15);%[1,5,6,7,11,12,13]); %15

%%
for roundNo=vecRound
    resultsBackupPath = [pathBackup, '_round_',num2str(roundNo),'/'];
    pathElimCost_old = [resultsBackupPath,'elimCosts.txt'];
    if isfile(pathElimCost_old)
        elimCostbyRound = readtable(pathElimCost_old);
    else
        elimCostbyRound=[];
    end
    
    paramBackupPath = [pathParams, 'round_', num2str(roundNo),'/'];
    paramsCt   = readtable([paramBackupPath,'allParametersSets_Ct.txt']);
    paramsNg   = readtable([paramBackupPath,'allParametersSets_Ng.txt']);
    paramsHIV  = readtable([paramBackupPath,'allParametersSets_HIV.txt']);
    paramsS    = readtable([paramBackupPath,'allParametersSets_syphilis.txt']);
    opts       = detectImportOptions([resultsBackupPath,'tabAlpha_last.txt']);
    %opts.VariableTypes(19:26)={'char'}; %decalé
    %opts.VariableTypes(32:52)={'char'}; %decalé
    opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
    tabAlpha   = readtable([resultsBackupPath,'tabAlpha_last.txt'],opts);
    paramRho   = readtable([resultsBackupPath,'paramRho.txt']);
    
    for Id_ech=unique(tabAlpha.IDech)'
        paramHIV_ech = paramsHIV(paramsHIV.IDech_id==Id_ech ,:);
        tabAlpha_ech = tabAlpha(tabAlpha.IDech==Id_ech ,:);
        for pHIV=unique(tabAlpha_ech.p)'
            
            if (1)%Id_ech>=11251838360004 && pHIV>0.1)|| Id_ech>11251838360004
                %             if (any(elimCostbyRound(elimCostbyRound.IDech==Id_ech,:).p==pHIV))
                %                 continue
                %             end
                log_path_id = [pathKnitroLog,'_round_',num2str(roundNo),'/sim_',num2str(Id_ech),'_',num2str(pHIV*100),'/'];
                mkdir(log_path_id)
                paramHIV_sim = paramsHIV(paramsHIV.IDech_id==Id_ech & paramsHIV.p==pHIV,:);
                if size(paramHIV_sim,1)>1
                    error('Size issue')
                end
                noEch = paramHIV_sim.nbEch;
                
                %-- Récupération des jeux de paramètres et des alphas calculés --%
                disp(['round=',num2str(roundNo),', b=',num2str(f),', p=',num2str(pHIV),', noEch=',num2str(noEch),', IDech=',num2str(Id_ech)])
                paramTab{1} = table2struct(paramsCt(paramsCt.IDech_id==Id_ech,:));
                paramTab{2} = table2struct(paramsNg(paramsNg.IDech_id==Id_ech,:));
                paramTab{3} = table2struct(paramHIV_sim);
                paramTab{4} = table2struct(paramsS(paramsS.IDech_id==Id_ech,:));
                mu = paramTab{3}.mu;
                b  = paramTab{3}.pi;
                
                alphas = tabAlpha(tabAlpha.IDech==Id_ech & tabAlpha.p==pHIV,:);
                
                elimCostbyModel = elimCostbyModelTemplate;
                
                %------------- Calcul des cnn par echantillon -------------%
                for numKit=1:length(list_kits_to_consider)
                    tic;
                    kit = list_kits_to_consider{numKit}; k = indexKit(kit); disp(['kit=',k])
                    
                    disp([' Calcul des c_i avec i dans le kit {',k,'}, alpha_k=',num2str(alphas.([k,'_',k]))])
                    doNotComputeCnn = 0;
                    for dis=k
                        if alphas.(['elim_',dis])
                            doNotComputeCnn = 1;
                        end
                    end
                    
                    if ~doNotComputeCnn
                        [cnn,~,msgCnn] = find_Cnn_kit(paramTab,mu,b,paramRho,f,kit,paramsC,alphas,afficherOutput,log_path_id,paramSolver);
                    else
                        for dis=kit
                            cnn.(dis{:}) = NaN;
                        end
                    end
                    elimCostbyModel(numKit,:).kit = indexKit(kit);
                    for dis=kit
                        elimCostbyModel(numKit,:).(dis{:}) = cnn.(dis{:});
                    end
                    
                    elimCostbyModel(numKit,:).msgCnnK    = msgCnn;
                    elimCostbyModel(numKit,:).timeCompil = round(toc,0);
                    elimCostbyModel(numKit,:).IDech = Id_ech;
                    elimCostbyModel(numKit,:).p = pHIV;
                    elimCostbyModel(numKit,:).f = f;
                    elimCostbyModel(numKit,:)
                    if numKit<length(list_kits_to_consider)
                        elimCostbyModel = [elimCostbyModel;elimCostbyModelTemplate];
                    end
                end
                %------------- fin du calcul des cnn -------------%
                
                elimCostbyRound = [elimCostbyRound; elimCostbyModel];
                %---------- end of assigning by strategies ---------%
                writetable(elimCostbyRound,[resultsBackupPath,'elimCosts.txt'],'WriteVariableNames',true)
            end%
        end
    end
end
