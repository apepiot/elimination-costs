% Issues identification and solving
clear all;
vecRounds   =  1:237; %[-[1,3:15,40:44,50:61],1:45,46,50,51,52];
paramFolder='3';
nameFile   = arrayAsACompactedString(vecRounds);
pathBackup = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\analysis_',paramFolder,'\rds_',nameFile,'\'];
fileCosts   = [pathBackup,'elimCosts_concatenated_rds_',nameFile,'.txt'];
fileAlphas   = [pathBackup,'tabAlpha_concatenated_rds_',nameFile,'.txt'];
opts       = detectImportOptions(fileCosts);
tableCosts_rdns = readtable(fileCosts,opts);
opts       = detectImportOptions(fileAlphas);
opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
tabAlpha      = readtable(fileAlphas,opts);

%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ALPHAS %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
histogram(tabAlpha.g_g)
tabAlpha(tabAlpha.g_g>0.7,:) %seems ok

onlyAlphas = table2array(tabAlpha(:,6:27));
round(max(onlyAlphas),3)>0.7
find(onlyAlphas==0.5);


%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% COSTS %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%



%% c_i==Inf or = 1 sur HIV
ech0 = tableCosts_rdns(tableCosts_rdns.HIV>=1,:);
ech = ech0;
paramFolder_set = paramFolder;
%usually it is because one infection has already been eliminated
addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
setupOnce;

% To compute P
verbose = 0;
paramSolver.tolP0 = 0.5e-4;
paramSolver.nbRelanceMax = 10;
paramSolver.timeLimit = 30; %seconds
log_path    = 'C:/Users/Moi/Desktop/Temporaire/tests';
ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
mod = 'hscg';

%To compute cnn
paramsC         = paramSolver;
paramsC.sup     = 1;
paramsC.inf     = -5;
paramsC.tolC    = 1e-3;
paramsC.iterMax = 60;
afficherOutput  = 0;
paramSolver.verbose       = 0;
paramSolver.sup_bnd_alpha = 20;
paramSolver.timeSolver    = 20;

for i=1:length(ech.IDech)
    currSet = ech(i,:);
    ID_ech  = currSet.IDech; %tester cette simu avec hscg, round 9
    pHIV    = currSet.p;
    roundNo = currSet.roundNo;
    f       = currSet.f;
    kit     = currSet.kit{:};
    
    %Reading parameters
    if roundNo<0
        roundNo = abs(roundNo);
        paramFolder_set = '1';
    else
        paramFolder_set = '2';
    end
    pathParam = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\paramSets_',paramFolder_set,'\round_',num2str(roundNo),'\'];
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
    pathRes = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\results_',paramFolder,'\_round_',num2str(roundNo),'\'];
    paramRho = table2struct(readtable([pathRes ,'paramRho.txt']));
    ampl = 0;
    [P,ES,msg,ampl] = P_mod_v8(paramTab,paramRho,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl);
    ampl.close();   
    
    disp(i)
    disp(P)
    problemFound=0;

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
        
        fileAlphas = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\results_',paramFolder,'\_round_',num2str(roundNo),'\tabAlpha.txt'];
        opts       = detectImportOptions(fileCosts);
        tableCosts_rdns = readtable(fileCosts,opts);
        opts       = detectImportOptions(fileAlphas);
        opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
        tabAlpha      = readtable(fileAlphas,opts);
        
        alphas = tabAlpha(tabAlpha.IDech==ID_ech & tabAlpha.p==pHIV,:);
        [cnn,~,msgCnn] = find_Cnn_kit_2(paramTab,mu,b,paramRho,f,kToKit(kit),paramsC,alphas,afficherOutput,log_path,paramSolver,[]);
        
        for inf=kit
            bigInf = cell2char(kToKit(inf));
            ech(i,:).(bigInf) = cnn.(bigInf);
        end
    end
end

% ajouter ces informations dans un fichier txt 
% tout ech
writetable(ech,[pathBackup,'_recalcul_costs_outsiders_high.txt'])



%% c_i<... (too low)
%histogram(tableCosts_rdns.HIV)
ech0 = tableCosts_rdns(tableCosts_rdns.HIV<=-0.5,:)

ech = ech0;
%usually it is because one infection has already been eliminated
setupOnce;

% To compute P
verbose = 0;
paramSolver.tolP0 = 0.5e-4;
paramSolver.nbRelanceMax = 10;
paramSolver.timeLimit = 30; %seconds
log_path    = 'C:/Users/Moi/Desktop/Temporaire/tests';
ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
mod = 'hscg';

%To compute cnn
paramsC         = paramSolver;
paramsC.sup     = 1;
paramsC.inf     = -5;
paramsC.tolC    = 1e-3;
paramsC.iterMax = 60;
afficherOutput  = 0;
paramSolver.verbose       = 0;
paramSolver.sup_bnd_alpha = 20;
paramSolver.timeSolver    = 20;

for i=1:length(ech.IDech)
    currSet = ech(i,:);
    ID_ech  = currSet.IDech; %tester cette simu avec hscg, round 9
    pHIV    = currSet.p;
    roundNo = currSet.roundNo;
    f       = currSet.f;
    kit     = currSet.kit{:};
    
    %Reading parameters
    if roundNo<0
        %roundNo = abs(roundNo);
        paramFolder_set = '1';
    else
        paramFolder_set='2';
    end
    pathParam = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\paramSets_',paramFolder_set,'\round_',num2str(abs(roundNo)),'\'];
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
    
    pathRes = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\results_',paramFolder,'\_round_',num2str(roundNo),'\'];
    paramRho = table2struct(readtable([pathRes ,'paramRho.txt']));
    ampl = 0;
    [P,ES,msg,ampl] = P_mod_v8(paramTab,paramRho,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl);
    ampl.close();   
    
    disp(i)
    disp(P)
    problemFound=0;

    for inf=kit
        if P.(inf)<10*paramSolver.tolP0
            disp('une des infections est bien eliminee')
            problemFound=1;
            ech(i,:).(['elim_',inf]) = 1;
        end
    end
    if ~problemFound
        fileAlphas = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\results_',paramFolder,'\_round_',num2str(roundNo),'\tabAlpha.txt'];
        opts       = detectImportOptions(fileCosts);
        tableCosts_rdns = readtable(fileCosts,opts);
        opts       = detectImportOptions(fileAlphas);
        opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
        tabAlpha      = readtable(fileAlphas,opts);
        %recalcul de c_i
        alphas = tabAlpha(tabAlpha.IDech==ID_ech & tabAlpha.p==pHIV,:);
        tic 
        [cnn,~,msgCnn] = find_Cnn_kit_2(paramTab,mu,b,paramRho,f,kToKit(kit),paramsC,alphas,afficherOutput,log_path,paramSolver,[]);
        toc;
        for inf=kit
            bigInf = cell2char(kToKit(inf));
            ech(i,:).(bigInf) = cnn.(bigInf);
        end
    end
end

% ajouter ces informations dans un fichier txt 
% tout ech
writetable(ech,[pathBackup,'_recalcul_costs_outsiders_low.txt'])

histogram(ech0.HIV,'BinWidth',0.4)
hold 
histogram(ech.HIV,'BinWidth',0.4)

