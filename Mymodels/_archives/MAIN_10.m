%% MAIN_10
clear all;
vecRounds_global = 100;
%-------------------------------------------------------
f               = 1;
VTunderART      = 1;
afficherOutput  = 0;
nosKitsAlpha    = 1:4;%15;%1:15;
nosKitsCnn      = 1:15;
calculer_alpha  = 0;
cond_part_A     = '1';%'pHIV==0.1 && noEch==2'; %'1' means all.
calculer_costs  = 1;
cond_part_Cnn   = '(Id_ech>=11251846140001 && pHIV>0.8)|| Id_ech>11251846140001';
ecraser_elimCosts = 0;
%--------------------------------------------------------

addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
parametrizationNo = 1;
list_kits = {{'HIV'},{'syphilis'},{'Ct'},{'Ng'},...
    {'HIV','syphilis'},{'HIV','Ct'},{'HIV','Ng'}, {'syphilis','Ct'},...
    {'syphilis','Ng'},{'Ct','Ng'}, {'HIV','syphilis','Ct'},{'HIV','syphilis','Ng'},...
    {'HIV','Ct','Ng'},{'syphilis','Ct','Ng'}, {'HIV','syphilis','Ct','Ng'}};


%Parameters for the calculation of alphas
ampl_models_dir = [pwd,'/AMPL_models/'];
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


addpath './Initialisation';
addpath './AMPL_models'
setupOnce;

for vecRounds=vecRounds_global
    if calculer_alpha
        disp(['###########== round ', num2str(vecRounds),' - calcul alpha ==###########'])
        list_kits_to_consider = list_kits(nosKitsAlpha);
        sensitivityAnalysis_v6_alpha;
    else
      warning('alphas not computed')
    end
    
    if calculer_costs
        disp(['###########== round ', num2str(vecRounds),' - calcul costs ==###########'])
        list_kits_to_consider = list_kits(nosKitsCnn);
        findCosts;
    end
    clear paramRho;
    
end