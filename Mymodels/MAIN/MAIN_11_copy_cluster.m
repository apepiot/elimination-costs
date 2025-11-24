%% MAIN_11 : alterne les echantillons, un dossier par k_mod
%clear all;
%-------------------------------------------------------
parametrizationNo = 11;
f                 = 1;
VTunderART        = 1;
afficherOutput    = 0;
vecRounds_global  = 1:10; %%.
kk                = 1;          %premier kit à lancer sur les 1:15
nosKitsAlpha      = 1;%[1,5,6,7,11,12,13,15]; %only HIV
nosKitsCnn        = nosKitsAlpha;
calculer_alpha    = 1;   %if f>1 : alphas already computed
cond_part_A       = '1'; %'pHIV==0.1 && noEch==2'; %'1' means all.
calculer_costs    = 0;
cond_part_Cnn     = '1';
ecraser_elimCosts = 1;
%--------------------------------------------------------

addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
%%.addpath '/home/pepiota/utils/ampl-z-13.1.20220703-Linux-64/amplapi/examples/matlab';
addpath './INITIALISATION/parametersSets/';
addpath './MAIN/AMPL_models'

list_kits = {{'HIV'},{'syphilis'},{'Ct'},{'Ng'},...
    {'HIV','syphilis'},{'HIV','Ct'},{'HIV','Ng'}, {'syphilis','Ct'},...
    {'syphilis','Ng'},{'Ct','Ng'}, {'HIV','syphilis','Ct'},{'HIV','syphilis','Ng'},...
    {'HIV','Ct','Ng'},{'syphilis','Ct','Ng'}, {'HIV','syphilis','Ct','Ng'}};


%Parameters for the calculation of alphas
ampl_models_dir = [pwd,'/MAIN/AMPL_models/']; 
log_path   = ['./Resultats_knitro/parametrization_',num2str(parametrizationNo),'/'];
backupPath = ['./ParameterAnalysis/results_',num2str(parametrizationNo),'/'];
paramSolverAlpha.tolP0        = 0.5e-4;
paramSolverAlpha.maxBndAlpha  = 50;
paramSolverAlpha.nbRelanceMax = 2;
paramSolverAlpha.timeLimit    = 30;
paramSolverAlpha.method_alpha = 'dicho';
paramSolverAlpha.tolAlpha     = 1e-4;
paramSolverAlpha.iterMaxDicho = 20;

%Parameters for the calculation of cnns
backUpFolder        = 'Runtest_1';
byModels            = 1;
paramSolver.tolP0   = paramSolverAlpha.tolP0;
paramsC.sup         = 1;
paramsC.inf         = -5;
paramsC.tolC        = 1e-3;
paramsC.iterMax     = 30;
paramSolver.verbose = afficherOutput;
paramSolver.sup_bnd_alpha = 20;
paramSolver.timeSolver    = 20;
pathBackup    = ['./ParameterAnalysis/results_',num2str(parametrizationNo),'/'];
pathParams    = ['./ParameterAnalysis/paramSets_',num2str(parametrizationNo),'/'];
pathKnitroLog = ['./Resultats_knitro/parametrization_',num2str(parametrizationNo),'/knitro_out_rounds/',backUpFolder,'/'];

setupOnce;
for vecRounds=vecRounds_global
    disp(kk)
    disp(list_kits{kk})
    %nkit = nosKitsCnn(kk);
    if calculer_alpha
        disp(['###########== round ', num2str(vecRounds),' - calcul alpha ==###########'])
        %list_kits_to_consider = combInCell(list_kits{kk}); %le kit en question et tous ses sous-kits
        list_kits_to_consider = {list_kits{kk}};
        sensitivityAnalysis_v7_alpha;
    else
      warning('not computing alphas')
    end
    
    if calculer_costs
        disp(['###########== round ', num2str(vecRounds),' - calcul costs ==###########'])
        list_kits_to_consider = {list_kits{kk}};
        findCosts;
    end
    clear paramRho;
    %kk=kk+1; %on passe au kit suivant au round suivant
    ikk = find(nosKitsCnn==kk);
    
    if ikk+1>length(nosKitsCnn)
        kk = nosKitsCnn(1);
    else
        kk = nosKitsCnn(ikk+1);
    end
    
end




