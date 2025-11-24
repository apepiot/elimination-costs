% New version : here we sample over each strategy %%%%PC%%%%%
% i.e. 1500 new set of parameters for each strategy
% strategies are independent
% we look for the cost of elimination of each kit
% difference with RHOHAT_v6: ?
% on parcourt strategie par strategie

% version of ODE_SICTPSEIISSEIIS2_v5 used

clear all;
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/';
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
addpath 'C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\Initialisation\';
setupOnce;

%% README
%---------------------------------------------------------------%
nbEchPerRound = 3; %1500 %echantillon de paramètres
nbRounds = 2;
vecPHIV=[0.1,0.2]; %ptage of people on PrEP
vecF = [1,2];
vecStratNos=1:15;
backUpFolder='Runtest10'; %mettre un nom de dossier qui n'existe pas
afficherOutput=0;
byModels = 1;
b=100; %input parameters
aprioriBndsC.sup=1;
tolP0 = 1e-3;
aprioriBndsC.inf=-1;
%---------------------------------------------------------------%
nbSimPerRound = length(vecF)*length(vecPHIV)*nbEchPerRound;
path2 = ['C:/Users/Moi/Documents/IPLESP/These/Codes/Simulations_StrategiesPC/',backUpFolder,'/'];
pathKnitroLog = ['C:/Users/Moi/Documents/Projets/Resultats_knitro/knitro_out_rounds/',backUpFolder,'/'];
ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
nbStrats = length(vecStratNos);

recapSimIDTemplate=NaN(nbSimPerRound,5);
tabElimTemplate = table('Size',[nbStrats*nbSimPerRound,10],...
    'VariableTypes',{'double','double','double','double','double','double','string','double','double','double'},...
    'VariableNames',{'IDech','strat','HIV','syphilis','Ct','Ng','code_err','timeCompil','f','p'});
elimCost_stratTemplate = table('Size',[1,5],...
    'VariableTypes',{'double','double','double','double','double'}, ...
    'VariableNames',{'HIV','syphilis','Ct','Ng','timeCompil'});
paramSolverAlpha.tolP0  = tolP0;
paramSolverAlpha.maxBndAlpha = 20;
paramSolverAlpha.nbRelanceMax = 5;
paramSolverAlpha.timeLimit = 20;  %for one start of the multistart
%%
%Liste des strategies
%1:Ct, 2:Ng, 3:HIV, 4:syph
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
list_kits_to_considerTemp = uniqueInCell([strategiesToConsiderTemp{:}]);

[status, msg, msgID] = mkdir([path2]); %where to store the results
if 1%isempty(msg)
    for roundNo=1:nbRounds
        fprintf(['---------------------------------------- \n ------------ Round number ',num2str(roundNo),' ------------ \n ---------------------------------------- \n'])
        
        %Initialisation des paramètres, tableaux de résulat du round
        simRndNo  = 0; %simulation number
        pathRound = [path2,'_round_',num2str(roundNo)]; [status, msg, msgID] = mkdir(pathRound);
        log_path  = [pathKnitroLog,'_round_',num2str(roundNo),'/knitro_out/'];
        recapSimID = recapSimIDTemplate;
        tabElim = tabElimTemplate;
        list_kits_to_consider = list_kits_to_considerTemp;
        tabAlpha = createTabAlpha(list_kits_to_consider,nbSimPerRound);
        
        for f=vecF
            for pHIV=vecPHIV
                disp(['round=',num2str(roundNo),', b=',num2str(f),', p=',num2str(pHIV)]) 
                for noSim=1:nbEchPerRound
                    simRndNo = simRndNo+1;
                    IDech_id = str2double([char(datetime('now','TimeZone','local','Format','yyMMddHHmmss')),num2str(simRndNo)]); %unique ID based on the current time
                    log_path_id = [log_path,num2str(IDech_id),'/'];
                    recapSimID(simRndNo,1:5) = [IDech_id,roundNo,f,pHIV,simRndNo];
                    
                    if byModels
                        elimCostbyModel = createelimCostbyModel(list_kits_to_consider);
                    end
                    disp('------------------------------------------')
                    disp(['New set of parameters: ',num2str(simRndNo),'/', num2str(nbSimPerRound),' with b=',num2str(f),' and p=',num2str(pHIV)])
                    [paramTab,mu,vecAlphas] = sampleParameters_v3_extent(true,true,true,true,b,pHIV);   %Ct,Ng,HIV,syph
                    createParamRho;
                    
                    %------------------------------- begin save parameters -----------------------------------%
                    paramSetsCt(simRndNo,:) = [struct2table(paramTab{1}),table(IDech_id)];
                    paramSetsNg(simRndNo,:) = [struct2table(paramTab{2}),table(IDech_id)];
                    paramSetsHIV(simRndNo,:) = [struct2table(paramTab{3}),table(IDech_id)];
                    paramSetsS(simRndNo,:) = [struct2table(paramTab{4}),table(IDech_id)];
                    writetable(paramSetsCt,[pathRound,'/allParametersSets_Ct'],'WriteVariableNames',true)
                    writetable(paramSetsNg,[pathRound,'/allParametersSets_Ng'],'WriteVariableNames',true)
                    writetable(paramSetsHIV,[pathRound,'/allParametersSets_HIV'],'WriteVariableNames',true)
                    writetable(paramSetsS,[pathRound,'/allParametersSets_syphilis'],'WriteVariableNames',true)
                    %----------------------------- end save parameters ----------------------------------------%
                    
                    %----------------------------- begin feasability check ----------------------------------------%
                    disp('  Checking the feasability of the parameter set')
                    [msgSolver,msgSol] = checkFeasability(paramTab,paramRho,b,mu,f,{'HIV','syphilis','Ct','Ng'},afficherOutput,tolP0);
                    disp(['  err_solver=',msgSolver])
                    %We delete the list of kits/models for which one of the diseases are already eliminated without VT.
                    %list_kits_to_consider = delKitsFeas(list_kits_to_consider,msgSol);
                    %------------------------------ end feasability check -----------------------------------------%
                    
                    disp('Calcul des alpha pour chaque kit')
                    if ~isempty(list_kits_to_consider) %au moins une des maladies est endémique
                        tic;
                        if byModels
                            tabAlpha(simRndNo,:).('IDech') = IDech_id;
                            
                            %------------- Calcul des alphas -------------%
                            for numKit=1:length(list_kits_to_consider)
                                kit = list_kits_to_consider{numKit}; k = indexKit(kit);
                                [alpha,~,msgAlpha_k] = findAlpha(paramTab,paramRho,mu,b,f,kit,afficherOutput,...
                                    log_path_id,paramSolverAlpha,ampl_models_dir);
                                %alpha=1; msgAlpha_k='-0';
                                alphas.(k) = alpha; msgAlphas.(k) = msgAlpha_k;
                                elimCostbyModel(numKit,:).alphaK    = alpha;
                                elimCostbyModel(numKit,:).msgAlphaK = msgAlpha_k;
                                tabAlpha(simRndNo,:).(k) = round(alpha,5);
                            end
                            writetable(tabAlpha,[path2,'_round_',num2str(roundNo),'/tabAlpha'],'WriteVariableNames',true)
                            %---------- fin du calcul des alphas ----------%
                            
                            
                            %------------- Calcul des cnn -------------%
                            for numKit=1:length(list_kits_to_consider)
                                tic;
                                kit = list_kits_to_consider{numKit}; k = indexKit(kit); disp(['kit=',k])
                                alpha = alphas.(k); msgAlpha_k = msgAlphas.(k);
                                if ~isequal(msgAlpha_k(end-1:end),'-0')    %alpha n'a pas été trouvé
                                    elimCostbyModel(numKit,:).kit = [kit{:}];
                                else
                                    disp([' Calcul des c_i avec i dans le kit {',k,'}, alpha_k=',num2str(alpha)])
                                    [cnn,~,msgCnn] = find_Cnn_kit(paramTab,mu,b,paramRho,f,kit,aprioriBndsC,alpha,afficherOutput,log_path_id);
                                    %cnn.HIV=3;cnn.syphilis=4; cnn.Ct=1;cnn.Ng=2;
                                    %msgCnn='0';
                                    elimCostbyModel(numKit,:).kit = [kit{:}];
                                    for dis=kit
                                        elimCostbyModel(numKit,:).(dis{:}) = cnn.(dis{:});
                                    end
                                end
                                elimCostbyModel(numKit,:).msgCnnK    = msgCnn;
                                elimCostbyModel(numKit,:).timeCompil = toc;
                                elimCostbyModel(numKit,:)
                            end
                        end
                        %------------- fin du calcul des cnn -------------%
                        
                        %------------- Calcul des srats -------------%
                        elimCostsByStrat = table('Size',[length(vecStratNos),5],...
                            'VariableNames',{'HIV','syphilis','Ct','Ng','timeCompil'},...
                            'VariableTypes',{'double', 'double', 'double', 'double', 'double'});
                        for noStrat=vecStratNos %strategies
                            strat = strategies{noStrat};
                            if ~byModels
                                tic;
                                [elimCost_strat] = findCnn_strategy(paramTab,paramRho,mu,b,f,strat,aprioriBndsC,vecAlphas,afficherOutput);
                                timeStrat = toc;
                            else
                                timeStrat=0;
                                elimCost_strat = elimCost_stratTemplate;
                                
%                                 for dis=["HIV","syphilis","Ct","Ng"]
%                                     disp(msgSol.(dis))
%                                     [strat{:}]
%                                 end
                                
                                
                                if 0%noResultForThisStrat
                                    %elimCost_strat(1,["HIV","syphilis","Ct","Ng","timeCompil"]) = array2table([NaN,NaN,NaN,NaN,NaN]);
                                else
                                    for kit=strat
                                        kit=kit{:};
                                        temp = elimCostbyModel(elimCostbyModel.kit==[kit{:}],:);
                                        timeStrat = timeStrat+temp.timeCompil;
                                        for dis=kit
                                            disp(dis)
                                            elimCost_strat(:,:).(dis{:}) = temp.(dis{:});
                                        end
                                    end
                                end
                            end
                            elimCostsByStrat(noStrat,["HIV","syphilis","Ct","Ng"]) = elimCost_strat(:,["HIV","syphilis","Ct","Ng"]);
                            elimCostsByStrat(noStrat,:).timeCompil = timeStrat;
                            %In the big table per round
                            tabElim((simRndNo-1)*nbStrats+noStrat,:)
                            
                            tabElim(((simRndNo-1)*nbStrats+noStrat),:).IDech = IDech_id;
                            tabElim(((simRndNo-1)*nbStrats+noStrat),:).strat = noStrat;
                            tabElim(((simRndNo-1)*nbStrats+noStrat),["HIV","syphilis","Ct","Ng"]) = elimCostsByStrat(noStrat,["HIV","syphilis","Ct","Ng"]);
                            tabElim((simRndNo-1)*nbStrats+noStrat,:).p = pHIV;
                            tabElim((simRndNo-1)*nbStrats+noStrat,:).f = f;
                            tabElim((simRndNo-1)*nbStrats+noStrat,:).timeCompil = round(elimCostsByStrat(noStrat,:).timeCompil,1);
                            
                        end
                    end
                    
                    writetable(array2table(recapSimID,'VariableNames', {'IDech','roundNo','f','p','simNo'}),[path2,'_round_',num2str(roundNo),...
                        '/recapSimID'],'WriteVariableNames',true)
                    fileName = ['_round_',num2str(roundNo),'/elim_cost.txt'];
                    writetable(tabElim,[path2,fileName],'WriteVariableNames',true)
                end
            end
        end
    end
else
    error([backUpFolder, ' already exists. Change the name of the back-up folder'])
end


