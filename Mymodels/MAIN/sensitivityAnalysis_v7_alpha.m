%difference with sensitivity_analysis_v6 : we only compute alpha of
%submodels of a given model/kit, when it is necessary

%clear all;
close all;
%parametrizationNo = 1;

%Reading parameterSets
%setsPath = ['./ParameterAnalysis/paramSets_',num2str(parametrizationNo)];
setsPath = pathParams;
%
paramSetsCt_tot=[];paramSetsNg_tot=[];paramSetsHIV_tot=[];paramSetsS_tot=[];
for roundNo = vecRounds
    pathFile = [setsPath,'round_',num2str(roundNo)];
    paramSetsCt  = readtable([pathFile,'/allParametersSets_Ct.txt']);
    paramSetsNg  = readtable([pathFile,'/allParametersSets_Ng']);
    paramSetsHIV = readtable([pathFile,'/allParametersSets_HIV']);
    paramSetsS   = readtable([pathFile,'/allParametersSets_syphilis']);
    
    paramSetsCt_tot = [paramSetsCt_tot;paramSetsCt];
    paramSetsNg_tot = [paramSetsNg_tot;paramSetsNg];
    paramSetsHIV_tot = [paramSetsHIV_tot;paramSetsHIV];
    paramSetsS_tot = [paramSetsS_tot;paramSetsS];
end

%% Calcul alpha avec le modèle à 4 infections
% ampl_models_dir = [pwd,'/AMPL_models/'];
% log_path = ['./Resultats_knitro/parametrization_',num2str(parametrizationNo),'/'];
% backupPath = ['./ParameterAnalysis/results_',num2str(parametrizationNo),'/'];
% paramSolverAlpha.tolP0  = 0.5e-4;
% paramSolverAlpha.maxBndAlpha = 50;
% paramSolverAlpha.nbRelanceMax = 2;
% paramSolverAlpha.timeLimit = 30;
% paramSolverAlpha.method_alpha = 'dicho';
% paramSolverAlpha.tolAlpha = 1e-3;
% paramSolverAlpha.iterMaxDicho = 20;
%
% f=1;
%
% afficherOutput=0;
%
% %addpath './';
% addpath './Initialisation';
% addpath './AMPL_models'
% setupOnce;

%attention a ne pas ecraser le fichier tabAlpha

createParamRho;
paramRho.rho_h = 0;
paramRho.rho_s = 0;
paramRho.rho_c = 0;
paramRho.rho_g = 0;
paramRho.VTunderART = VTunderART;

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
    %elimCostbyModel = createelimCostbyModel(list_kits_to_consider);
    tabAlpha        = createTabAlpha2(list_kits_to_consider,0);
    
    vecID_ech = unique(paramCt_round.IDech_id);
    for id_ech=vecID_ech'
        paramCt_ech  = paramCt_round(paramCt_round.IDech_id==id_ech,:);
        paramNg_ech  = paramNg_round(paramNg_round.IDech_id==id_ech,:);
        paramS_ech   = paramS_round(paramS_round.IDech_id==id_ech,:);
        paramHIV_ech = paramHIV_round(paramHIV_round.IDech_id==id_ech,:);
        vecP = unique([paramHIV_ech.p]);  % j avais ajouté 0  mais je l'enleve
        b  = paramS_ech.pi;
        mu = paramS_ech.mu;
        
        for p=vecP'
            noEch = paramCt_ech.nbEch;
            pHIV  = p;

            if eval(cond_part_A)
                disp(['ID_ech=',num2str(id_ech),', p=',num2str(p),', nEch=',num2str(noEch)]);
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
                tabAlpha = [tabAlpha;array2table(nan(1,size(tabAlpha,2)),'VariableNames',tabAlpha.Properties.VariableNames)];
                %alphas = array2table(nan(1,size(tabAlpha,2)),'VariableNames',tabAlpha.Properties.VariableNames);
                %alphas.IDech                = id_ech;
                tabAlpha(numSim,:).IDech    = id_ech;
                %alphas.roundNo              = roundNo;
                tabAlpha(numSim,:).roundNo  = roundNo;
                %alphas.nbEch                = nbEch;
                tabAlpha(numSim,:).nbEch    = nbEch;
                %alphas.p = p;
                tabAlpha(numSim,:).p        = p;
                msg = [];
                
                %check Feasability
                [infElim,msgSolver,msgSol] = checkFeasability_v2(paramTab,paramRho,b,mu,f,'hscg',afficherOutput,paramSolverAlpha,log_path_round,ampl_models_dir);
                
                tabAlpha(numSim,["elim_h","elim_s","elim_c","elim_g"]) = struct2table(infElim);
                
                tic;
                                
                for numKit=1:length(list_kits_to_consider)
                    kit = list_kits_to_consider{numKit}; k = indexKit(kit);
                    if ~contains(k,'h')
                        mod = ['h',k]; k_mod = [k,'_',mod];
                    else
                        mod = k; k_mod = [k,'_',k];
                    end
                    doNotRunModel=0;
                    for inf='hscg'
                        if contains(mod,inf) && infElim.(inf)
                            doNotRunModel = 1;
                        end
                    end
                    if doNotRunModel
                        alpha        =  0;
                        msgAlpha_mod = '1';
                        elim_i       = '0';
                        tabAlpha(numSim,:).([k_mod,'_elim'])= elim_i;%table({elim_i}, 'VariableNames',{[k_mod,'_elim']});
                        tabAlpha(numSim,:).(k_mod)      = alpha;
                        msg                             = [msg, ' ',msgAlpha_mod];  
                    else            %on va calculer pour tous les sous-modeles
                        k_suivant       = k;
                        k_mod_suivant   = k_mod;
                        kit_suivant     = kToKit(k_suivant);
                        disp(k_mod_suivant);
                        
                        %Est-ce que les alphas du sous-kit ont deja ete calcules ?
                        isTableCol = @(thisCol) ismember(thisCol, tabAlpha.Properties.VariableNames);
                        while ~isempty(k_suivant)
                            if ~isTableCol(k_mod_suivant)  %la colonne n'exist pas 
                                k_suivant_non_calcule = 1;
                            else
                                if isnan(tabAlpha(numSim,:).(k_mod_suivant))
                                    k_suivant_non_calcule = 1;
                                else
                                    k_suivant_non_calcule = 0; 
                                    k_suivant = ''; %on ne continue pas la boucle
                                end
                            end
                            
                            if k_suivant_non_calcule
                                [alpha,P,~,msgAlpha_mod,elim_i] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit_suivant,k_mod_suivant,afficherOutput,...
                                    log_path_round,paramSolverAlpha,ampl_models_dir);
                                tabAlpha(numSim,:).([k_mod_suivant,'_elim'])= {elim_i};%table({elim_i}, 'VariableNames',{[k_mod_suivant,'_elim']});
                                tabAlpha(numSim,:).(k_mod_suivant)          = alpha;
                                msg                                         = [msg, ' ',msgAlpha_mod];
                                for inf=elim_i
                                    k_suivant       = erase(k_suivant,inf);
                                    k_mod_suivant   = erase(k_mod_suivant,inf);
                                    kit_suivant     = kToKit(k_suivant);
                                end
                            end
                        end
                    end
                    
%                     if ~contains(k,'h')
%                         k_mod2 = [k,'_',k];
%                         lineAlpha0 = tabAlpha(tabAlpha.IDech==id_ech & tabAlpha.p==0 & tabAlpha.roundNo==roundNo,:);
%                         if ~infElim.h && ~isempty(lineAlpha0) %if HIV is not eliminated, the model is similar to the one with hiv but p=0 (i.e. sg_sg is similar to sg_hsg(p=0))
%                             tabAlpha(numSim,[k_mod2,'_elim'])  = lineAlpha0.([k_mod,'_elim']);
%                             tabAlpha(numSim,:).(k_mod2)        = lineAlpha0.(k_mod);
%                         else           %if HIV is eliminated, idem si ~infElim.h normalement
%                             [alpha,P,~,msgAlpha_mod,elim_i] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit,k_mod2,afficherOutput,...
%                                 log_path_round,paramSolverAlpha,ampl_models_dir);
%                             tabAlpha(numSim,[k_mod2,'_elim'])  = table({elim_i}, 'VariableNames',{[k_mod2,'_elim']});
%                             tabAlpha(numSim,:).(k_mod2)        = alpha;
%                         end
%                         msg = [msg, ' ',msgAlpha_mod];
%                     end
                end
                
                tabAlpha(numSim,:).timeCompil = round(toc,2);
                tabAlpha(numSim,:).msg = msg;
                
                %alphas
                %tabAlpha = join(tabAlpha,alphas);%,{'roundNo','p','IDech','nbEch'});
            end
            writetable(tabAlpha,[backupPath,'_round_',num2str(roundNo),'/tabAlpha'],'WriteVariableNames',true)
        end
    end
end