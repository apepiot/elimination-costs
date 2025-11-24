clear all;
%------------------------------
vecRounds    =  1:237;      %[-[1,3:15,40:44,50:61],1:45,46,50,51,52];
paramFolder  = '6';
nameFile     = arrayAsACompactedString(vecRounds);
mainPath     = './ParameterAnalysis/';
pathMatlab   = '/home/pepiota/utils/ampl-z-13.1.20220703-Linux-64/amplapi/examples/matlab/';
%addpath '/home/pepiota/utils/ampl-z-13.1.20220703-Linux-64/amplapi/examples/matlab/';
ampl_models_dir = [pwd,'/AMPL_models/'];
%addpath '/home/pepiota/utils/ampl-z-13.1.20220703-Linux-64/amplapi/matlab/';
lower_bound = -2;
upper_bound = 1;
%addpath '/home/pepiota/utils/ampl-z-13.1.20220703-Linux-64:/home/pepiota/utils/knitro-14.0.0-Linux64/knitroampl';
%-------------------------------
pathAnalysis = [mainPath,'analysis_',paramFolder,'/rds_',nameFile,'/'];
pathResults  = [mainPath,'results_',paramFolder,'/'];
pathParameters= [mainPath,'paramSets_',paramFolder,'/'];
fileCosts    = [pathAnalysis,'elimCosts_concatenated_rds_',nameFile,'.txt'];
fileAlphas   = [pathAnalysis,'tabAlpha_concatenated_rds_',nameFile,'.txt'];
log_path     = './logs';
opts         = detectImportOptions(fileCosts);
cost_table   = readtable(fileCosts,opts);
opts         = detectImportOptions(fileAlphas);
opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
tabAlpha     = readtable(fileAlphas,opts);
%addpath './AMPL_models/';
