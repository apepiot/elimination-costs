% Issues identification and solving
clear all;
%------------------------------
vecRounds    =  1:320;      %[-[1,3:15,40:44,50:61],1:45,46,50,51,52];
paramFolder  = '610';
nameFile     = arrayAsACompactedString(vecRounds);
mainPath     = '..\ParameterAnalysis\';
pathMatlab  = 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
%pathMatlab   = 'C:\Users\pepiot\Documents\Artelys\AMPL-13.1.20220703-Win-64\amplapi\examples\matlab';
ampl_models_dir = 'C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\MAIN\AMPL_models\';
%ampl_models_dir = 'C:\Users\pepiot\Documents\PhD\codes\multi-voluntary-testing\Mymodels\MAIN\AMPL_models\';
addpath('..\INITIALISATION\')
addpath('..\RESULTS_ANALYSIS\')
addpath('..\MAIN\')
addpath('..')

lower_bound = -2;
upper_bound = 1;
%-------------------------------
pathAnalysis = [mainPath,'analysis_',paramFolder,'\rds_',nameFile,'\'];
pathResults  = [mainPath,'results_',paramFolder,'\'];
pathParameters= [mainPath,'paramSets_',paramFolder,'\'];
fileCosts    = [pathAnalysis,'elimCosts_concatenated_rds_',nameFile,'.txt'];
fileAlphas   = [pathAnalysis,'tabAlpha_concatenated_rds_',nameFile,'.txt'];
log_path     = '~/logs';
opts_c       = detectImportOptions(fileCosts);
cost_table   = readtable(fileCosts,opts_c);
opts_a       = detectImportOptions(fileAlphas);
opts.VariableTypes(contains(opts_a.VariableNames,'_elim')) = {'char'};
tabAlpha     = readtable(fileAlphas,opts_a);

paramSolverAlpha.tolP0        = 0.5e-4;
paramSolverAlpha.maxBndAlpha  = 50;
paramSolverAlpha.nbRelanceMax = 2;
paramSolverAlpha.timeLimit    = 30;
paramSolverAlpha.method_alpha = 'dicho';
paramSolverAlpha.tolAlpha     = 1e-4;
paramSolverAlpha.iterMaxDicho = 20;

%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ALPHAS %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%histogram(tabAlpha.g_g)
%tabAlpha(tabAlpha.g_g>0.7,:) %seems ok

%onlyAlphas = table2array(tabAlpha(:,6:27));
%round(max(onlyAlphas),3)>0.7
%find(onlyAlphas==0.5);


%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%% COSTS %%%%%%%
%%%%%%%%%%%%%%%%%%%%%

% If already exist, use the costs that have already been computed
pathfileCorr = [pathAnalysis,'_recalcul_costs_outsiders.txt'];
if isfile(pathfileCorr)
    selected_last_time = readtable(pathfileCorr,opts_c);
    %selected_last_time.Var18 =  [];
    %selected_last_time.roundNo_old = selected_last_time.roundNo;
    already_corrected  = deleteCopies(selected_last_time);%(selected_last_time.HIV<upper_bound & selected_last_time.HIV>lower_bound,:);
    updated_cost_table = replaceLineInTable(cost_table,already_corrected);
else
    updated_cost_table = cost_table;
    selected_last_time = [];
    already_corrected  = [];
end

%% c_i==Inf or = 1 sur HIV or nan (and infections of the kit not already eliminated)
ech0 = cost_table(cost_table.HIV>=upper_bound | cost_table.HIV<=lower_bound,:);  %cost that should have been corrected at the beginning
disp(size(ech0))
ech  = updated_cost_table(updated_cost_table.HIV>=upper_bound | updated_cost_table.HIV<=lower_bound | isnan(updated_cost_table.HIV),:); %[ech,already_corrected] = ech0
disp(size(ech))
n_nan = sum(isnan(ech.HIV));
disp(size(ech,1)-n_nan)
ech(~isnan(ech.HIV),:)

addpath(pathMatlab); 
setupOnce;

% To compute P
verbose = 0;
paramSolver.tolP0 = 0.5e-4;
paramSolver.nbRelanceMax = 10;
paramSolver.timeLimit = 30; %seconds
mod = 'hscg';

%To compute cnn
paramsC         = paramSolver;
paramsC.sup     = 1;
paramsC.inf     = -2;
paramsC.tolC    = 1e-3;
paramsC.iterMax = 10;
afficherOutput  = 0;
paramSolver.verbose       = 0;
paramSolver.sup_bnd_alpha = 20;
paramSolver.timeSolver    = 20;
iter_not_nan=1;
for i=1:length(ech.IDech)
    disp([num2str(i),'/',num2str(length(ech.IDech)), '(tot), ',num2str(iter_not_nan),'/',num2str(length(ech.IDech)-n_nan)])

    currSet = ech(i,:);
    ID_ech  = currSet.IDech; 
    pHIV    = currSet.p;
    roundNo = currSet.roundNo;
    f       = currSet.f;
    kit     = currSet.kit{:};
   
    problemFound = 0;

    for inf=kit
        if ech(i,:).(['elim_',inf]) == 1
            %disp('une des infections est deja eliminee')
            problemFound=1;
        end
    end
    
    if problemFound
        continue;
    else
       iter_not_nan=iter_not_nan+1;
    end
    
    %Reading parameters    
    pathParam = [pathParameters,'\round_',num2str(roundNo),'\'];
    paramH = readtable([pathParam, 'allParametersSets_HIV.txt']);
    paramS = readtable([pathParam, 'allParametersSets_syphilis.txt']);
    paramC = readtable([pathParam, 'allParametersSets_Ct.txt']);
    paramG = readtable([pathParam, 'allParametersSets_Ng.txt']);
    paramTab{1} = table2struct(paramC(paramC.IDech_id==ID_ech,:));
    paramTab{2} = table2struct(paramG(paramG.IDech_id==ID_ech,:));
    paramTab{3} = table2struct(paramH(paramH.IDech_id==ID_ech & paramH.p==pHIV,:));
    paramTab{4} = table2struct(paramS(paramS.IDech_id==ID_ech,:));
    mu = paramTab{3}.mu; b = paramTab{3}.pi;
    
    %Prevalence computation at rho_k =0
    pathRes = [pathResults,'_round_',num2str(roundNo),'\'];
    paramRho = table2struct(readtable([pathRes ,'paramRho.txt']));
    ampl = 0;
    [P,ES,msg,ampl] = P_mod_v8(paramTab,paramRho,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl);
    ampl.close();   
    disp(P)
    for inf=kit
        if P.(inf)<10*paramSolver.tolP0
            disp('une des infections est bien eliminee')
            problemFound=1;
            ech(i,:).(['elim_',inf]) = 1;
            ech(i,:).HIV = NaN;
            ech(i,:).syphilis = NaN;
            ech(i,:).Ct = NaN;
            ech(i,:).Ng = NaN;
        end
    end
    if ~problemFound
        %recalcul de c_i
        
        fileAlphas = [pathResults,'\_round_',num2str(roundNo),'\tabAlpha.txt'];
        opts       = detectImportOptions(fileCosts);
        cost_table = readtable(fileCosts,opts);
        opts       = detectImportOptions(fileAlphas);
        opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
        tabAlpha   = readtable(fileAlphas,opts);
        
        alphas = tabAlpha(tabAlpha.IDech==ID_ech & tabAlpha.p==pHIV,:); disp(alphas)
        
        if isnan(alphas.h_h) || isnan(alphas.([kit,'_',kit]))
            list_k_mod = opts.VariableNames(contains(opts.VariableNames,'_elim'));
            for j = 1:length(list_k_mod)
                k_mod = erase(list_k_mod{j},'_elim');
                [alpha_k_mod,~,~,~,~] = findAlpha_v4(paramTab,paramRho,mu,b,f,kToKit(kit),k_mod,verbose,log_path,paramSolverAlpha,ampl_models_dir);
                alphas.(k_mod) = alpha_k_mod;
            end
            disp(alphas)
        end
        
        [cnn,~,msgCnn] = find_Cnn_kit_2(paramTab,mu,b,paramRho,f,kToKit(kit),paramsC,alphas,afficherOutput,log_path,paramSolver,[],ampl_models_dir);
        
        for inf=kit
            bigInf = cell2char(kToKit(inf));
            ech(i,:).(bigInf) = cnn.(bigInf);
        end
    end
    writetable([ech;already_corrected],[pathAnalysis,'_recalcul_costs_outsiders.txt'])
end

updated_cost_table = replaceLineInTable(updated_cost_table,ech);
%updated_cost_table(updated_cost_table.HIV>1000,:).HIV = NaN;

writetable(updated_cost_table,[pathAnalysis,'elimCosts_concatenated_rds_',nameFile,'_correction.txt'])

histogram(ech0.HIV,'BinWidth',0.4)
hold on
newEch = [ech;already_corrected];
histogram(newEch.HIV,'BinWidth',0.4)
