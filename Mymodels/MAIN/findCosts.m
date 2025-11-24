%=MAIN_9 in my own files
% difference with MAIN_8 and MAIN_8bis: we focus only on HIV models (not strategies)

% clear all;
% %addpath '';
% addpath './AMPL_models/';
% addpath './Initialisation/';
% setupOnce;

%f=1;

% backUpFolder='Runtest_1'; %mettre un nom de dossier qui n'existe pas
% afficherOutput=0;
% byModels = 1;
% paramSolver.tolP0 = 0.5e-4;
% paramsC.sup=1;
% paramsC.inf=-15;
% paramsC.tolC=1e-3;
% paramsC.iterMax=20;
% paramSolver.verbose = afficherOutput;
% paramSolver.sup_bnd_alpha = 20;
% parametrizationNo = 1;
% paramSolver.timeSolver = 30;
% %---------------------------------------------------------------%
% pathBackup = ['./ParameterAnalysis/results_',num2str(parametrizationNo),'/'];
% pathParams = ['./ParameterAnalysis/paramSets_',num2str(parametrizationNo),'/'];
%
% ampl_models_dir = [pwd,'/AMPL_models/'];
% pathKnitroLog = ['./Resultats_knitro/parametrization_',num2str(parametrizationNo),'/knitro_out_rounds/',backUpFolder,'/'];

%elimCost_roundTemplate = table('Size',[1,10],...
%    'VariableTypes',{'double','double','double','double','double','double','string','double','double','double'},...
%    'VariableNames',{'IDech','strat','HIV','syphilis','Ct','Ng','code_err','timeCompil','f','p'});

elimCostbyModelTemplate = table('Size',[1,12],...
    'VariableTypes',{'double','string','double','double','double','double','string','double','double','double','double','double'},...
    'VariableNames',{'IDech','kit','HIV','syphilis','Ct','Ng','msgCnnK','timeCompil','f','p','roundNo','noEch'});
elimCostbyModelTemplate.HIV      = NaN;
elimCostbyModelTemplate.syphilis = NaN;
elimCostbyModelTemplate.Ct       = NaN;
elimCostbyModelTemplate.Ng       = NaN;


%%
for roundNo=vecRounds
    resultsBackupPath = [pathBackup, '_round_',num2str(roundNo),'/'];
    pathElimCost_old = [resultsBackupPath,'elimCosts_f_',num2str(f),'.txt'];
    if isfile(pathElimCost_old) && ~ecraser_elimCosts
        elimCostbyRound = readtable(pathElimCost_old);
    else
        elimCostbyRound=[];
    end
    
    paramBackupPath = [pathParams, 'round_', num2str(roundNo),'/'];
    paramsCt   = readtable([paramBackupPath,'allParametersSets_Ct.txt']);
    paramsNg   = readtable([paramBackupPath,'allParametersSets_Ng.txt']);
    paramsHIV  = readtable([paramBackupPath,'allParametersSets_HIV.txt']);
    paramsS    = readtable([paramBackupPath,'allParametersSets_syphilis.txt']);
    opts       = detectImportOptions([resultsBackupPath,'tabAlpha.txt']);
    opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
    tabAlpha   = readtable([resultsBackupPath,'tabAlpha.txt'],opts);
    paramRho   = readtable([resultsBackupPath,'paramRho.txt']);
    
    for Id_ech=unique(tabAlpha.IDech)'
        paramHIV_ech = paramsHIV(paramsHIV.IDech_id==Id_ech ,:);
        tabAlpha_ech = tabAlpha(tabAlpha.IDech==Id_ech ,:);
        for pHIV=unique(tabAlpha_ech.p)'
            paramHIV_sim = paramsHIV(paramsHIV.IDech_id==Id_ech & paramsHIV.p==pHIV,:);
            if size(paramHIV_sim,1)>1
                error('Size issue')
            end
            noEch = paramHIV_sim.nbEch;
            
            if (eval(cond_part_Cnn))
                %             if (any(elimCostbyRound(elimCostbyRound.IDech==Id_ech,:).p==pHIV))
                %                 continue
                %             end
                log_path_id = [pathKnitroLog,'_round_',num2str(roundNo),'/sim_',num2str(Id_ech),'_',num2str(pHIV*100),'/'];
                mkdir(log_path_id)

                
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
                    ticS = tic;
                    kit = list_kits_to_consider{numKit}; k = indexKit(kit); disp(['kit=',k])
                    
                    disp([' Calcul des c_i avec i dans le kit {',k,'}, alpha_k=',num2str(alphas.([k,'_',k]))])
                    doNotComputeCnn = 0;
                    for dis=k
                        if alphas.(['elim_',dis])
                            doNotComputeCnn = 1;
                        end
                    end
                    
                    if ~doNotComputeCnn
                        %aprioriC = [];
                        aprioriC = retrieveC(elimCostbyModel(elimCostbyModel.IDech==Id_ech & elimCostbyModel.p==pHIV & elimCostbyModel.f==f,:),alphas,k);
                        %[cnn,~,msgCnn] = find_Cnn_kit(paramTab,mu,b,paramRho,f,kit,paramsC,alphas,afficherOutput,log_path_id,paramSolver);
                        [cnn,~,msgCnn] = find_Cnn_kit_2(paramTab,mu,b,paramRho,f,kit,paramsC,alphas,afficherOutput,log_path_id,paramSolver,aprioriC,ampl_models_dir);
                    else
                        for dis=kit
                            cnn.(dis{:}) = NaN;
                        end
                        msgCnn='NaN';
                    end
                    elimCostbyModel(numKit,:).kit = indexKit(kit);
                    for dis=kit
                        elimCostbyModel(numKit,:).(dis{:}) = round(cnn.(dis{:}),5);
                    end
                    
                    elimCostbyModel(numKit,:).msgCnnK    = msgCnn;
                    elimCostbyModel(numKit,:).timeCompil = round(toc(ticS),1);
                    elimCostbyModel(numKit,:).IDech = Id_ech;
                    elimCostbyModel(numKit,:).p = pHIV;
                    elimCostbyModel(numKit,:).f = f;
                    elimCostbyModel(numKit,:).roundNo = roundNo;
                    elimCostbyModel(numKit,:).noEch = noEch;
                    elimCostbyModel(numKit,:)
                    if numKit<length(list_kits_to_consider)
                        elimCostbyModel = [elimCostbyModel;elimCostbyModelTemplate];
                    end
                end
                disp('- - - - - - - - - - - -')
                %------------- fin du calcul des cnn -------------%
                
                elimCostbyRound = [elimCostbyRound; elimCostbyModel];
                %---------- end of assigning by strategies ---------%
                writetable(elimCostbyRound,[resultsBackupPath,'elimCosts_f_',num2str(f),'.txt'],'WriteVariableNames',true)
            end%
        end
    end
end
