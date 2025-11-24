clear all;
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/';
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/Initialisation/';
setupOnce;

f=1;
vecRound = 1:2;
vecStratNos=1:3;
backUpFolder='Runtest_1'; %mettre un nom de dossier qui n'existe pas
afficherOutput=0;
byModels = 1;
aprioriBndsC.sup=1;
tolP0 = 1e-3;
aprioriBndsC.inf=-1;
parametrizationNo = 1;
alphaMax = 20;
%---------------------------------------------------------------%
pathBackup = ['C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/ParameterAnalysis/results_',num2str(parametrizationNo),'/'];
pathParams = ['C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/ParameterAnalysis/paramSets_',num2str(parametrizationNo),'/'];
%%
pathKnitroLog = ['C:/Users/Moi/Documents/Projets/Resultats_knitro/knitro_out_rounds/',backUpFolder,'/'];
ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';

elimCost_roundTemplate = table('Size',[1,10],...
    'VariableTypes',{'double','double','double','double','double','double','string','double','double','double'},...
    'VariableNames',{'IDech','strat','HIV','syphilis','Ct','Ng','code_err','timeCompil','f','p'});
elimCost_stratTemplate = table('Size',[1,7],...
    'VariableTypes',{'double','double','double','double','double','double','double'}, ...
    'VariableNames',{'HIV','syphilis','Ct','Ng','timeCompil','IDech','p'});
elimCostbyModelTemplate = table('Size',[1,10],...
    'VariableTypes',{'double','string','double','double','double','double','string','double','double','double'},...
    'VariableNames',{'IDech','kit','HIV','syphilis','Ct','Ng','msgCnnK','timeCompil','f','p'});
%%
strategies{1} = {{'HIV'},{'syphilis'},{'Ct'},{'Ng'}};
strategies{2} = {{'HIV'},{'syphilis'},{'Ct','Ng'}};
strategies{3} = {{'HIV','Ct'},{'syphilis'},{'Ng'}};
strategies{4} = {{'HIV'},{'syphilis','Ct'},{'Ng'}};
strategies{5} = {{'HIV','Ng'},{'syphilis'},{'Ct'}};
strategies{6} = {{'HIV'},{'syphilis','Ng'},{'Ct'}};
strategies{7} = {{'HIV','syphilis'},{'Ct'},{'Ng'}};
strategies{8} = {{'HIV','syphilis'},{'Ct','Ng'}};
strategies{9} = {{'HIV','Ct'},{'syphilis','Ng'}};
strategies{10} = {{'HIV','Ng'},{'syphilis','Ct'}};
strategies{11} = {{'HIV','Ct','Ng'},{'syphilis'}};
strategies{12} = {{'HIV'},{'syphilis','Ct','Ng'}};
strategies{13} = {{'HIV','syphilis','Ct'},{'Ng'}};
strategies{14} = {{'HIV','syphilis','Ng'},{'Ct'}};
strategies{15} = {{'HIV','syphilis','Ct','Ng'}};

list_kits = {{'HIV'},{'syphilis'},{'Ct'},{'Ng'},...
    {'HIV','syphilis'},{'HIV','Ct'},{'HIV','Ng'}, {'syphilis','Ct'},...
    {'syphilis','Ng'},{'Ct','Ng'}, {'HIV','syphilis','Ct'},{'HIV','syphilis','Ng'},...
    {'HIV','Ct','Ng'},{'syphilis','Ct','Ng'}, {'HIV','syphilis','Ct','Ng'}};
strategiesToConsiderTemp = {strategies{vecStratNos}};
list_kits_to_consider = uniqueInCell([strategiesToConsiderTemp{:}]);
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
    tabAlpha   = readtable([resultsBackupPath,'tabAlpha_all_v1.txt']);
        
    for Id_ech=unique(paramsHIV.IDech_id)'
        paramHIV_ech = paramsHIV(paramsHIV.IDech_id==Id_ech ,:);
        for pHIV=unique(paramHIV_ech.p)'
            if (any(elimCostbyRound(elimCostbyRound.IDech==Id_ech,:).p==pHIV))
                continue
            end
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
            mu=paramTab{3}.mu; b=paramTab{3}.pi;
            createParamRho; %tout à 0
            alphas      = tabAlpha(tabAlpha.IDech==Id_ech & tabAlpha.p==pHIV ,:);
            elimCostbyModel = elimCostbyModelTemplate;
            
            %------------- Calcul des cnn par echantillon -------------%
            for numKit=1:length(list_kits_to_consider)
                tic;
                kit = list_kits_to_consider{numKit}; k = indexKit(kit); disp(['kit=',k])
                alpha = alphas.(k); %msgAlpha_k = msgAlphas.(k);
                if abs(alpha-alphaMax)>=10^(-3)
                    %if ~isequal(msgAlpha_k(end-1:end),'-0')    %alpha n'a pas été trouvé
                    %    elimCostbyModel(numKit,:).kit = [kit{:}];
                    %else
                    disp([' Calcul des c_i avec i dans le kit {',k,'}, alpha_k=',num2str(alpha)])
                    [cnn,~,msgCnn] = find_Cnn_kit(paramTab,mu,b,paramRho,f,kit,aprioriBndsC,alpha,afficherOutput,log_path_id);
                    elimCostbyModel(numKit,:).kit = [kit{:}];
                    for dis=kit
                        elimCostbyModel(numKit,:).(dis{:}) = cnn.(dis{:});
                    end
                    %end
                else
                    disp([' On ne calcule pas les c_i dans le kit {',k,'}, car alpha_k=',num2str(alpha)])
                    elimCostbyModel(numKit,:).kit = [kit{:}];
                    for dis=kit
                        elimCostbyModel(numKit,:).(dis{:}) = Inf;
                    end
                end
                elimCostbyModel(numKit,:).msgCnnK    = msgCnn;
                elimCostbyModel(numKit,:).timeCompil = toc;
                elimCostbyModel(numKit,:).IDech = Id_ech;
                elimCostbyModel(numKit,:).p = pHIV;
                elimCostbyModel(numKit,:).f = f;
                elimCostbyModel(numKit,:)
                elimCostbyModel = [elimCostbyModel;elimCostbyModelTemplate];
                
            end
            %------------- fin du calcul des cnn -------------%
            
            %------------ Assigning by strategies ----------%
            elimCostsByStrat = table('Size',[length(vecStratNos),9],...
                'VariableNames',{'HIV','syphilis','Ct','Ng','timeCompil','IDech','p','errStrat','strat'},...
                'VariableTypes',{'double', 'double', 'double', 'double', 'double','double','double','string','double'});
            for noStrat=vecStratNos %strategies
                strat = strategies{noStrat}; timeStrat= 0 ;
                for kit=strat
                    kit=kit{:};
                    temp = elimCostbyModel(elimCostbyModel.kit==[kit{:}],:);
                    timeStrat = timeStrat+temp.timeCompil;
                    msgStrat = [];
                    for dis=kit
                        disp(dis)
                        elimCost_strat.(dis{:}) = temp.(dis{:});
                    end
                    msgStrat = [msgStrat,temp.msgCnnK];
                end
                elimCostsByStrat(noStrat,:).HIV      =  elimCost_strat.HIV;
                elimCostsByStrat(noStrat,:).syphilis =  elimCost_strat.syphilis;
                elimCostsByStrat(noStrat,:).Ct       =  elimCost_strat.Ct;
                elimCostsByStrat(noStrat,:).Ng       =  elimCost_strat.Ng;
                elimCostsByStrat(noStrat,:).timeCompil = timeStrat;
                elimCostsByStrat(noStrat,:).errStrat = msgStrat;
                elimCostsByStrat(noStrat,:).strat = noStrat;
                elimCostsByStrat(noStrat,:).IDech = Id_ech;
                elimCostsByStrat(noStrat,:).p = pHIV; 
                
                elimCostbyRound = [elimCostbyRound; ...
                    [elimCostsByStrat(noStrat,["HIV","syphilis","Ct","Ng","timeCompil","IDech","p","errStrat"]),...
                    table(noStrat,f,noEch,roundNo,'VariableNames',{'strat','f','noEch','roundNo'})]];
            end
            %---------- end of assigning by strategies ---------%
            writetable(elimCostbyRound,[resultsBackupPath,'elimCosts.txt'],'WriteVariableNames',true)
        end%
    end
end
