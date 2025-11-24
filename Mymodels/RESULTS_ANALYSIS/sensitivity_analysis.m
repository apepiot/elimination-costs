
close all; clear all;
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\MAIN')
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\INITIALISATION');

%% MERGE parameter Sets by infection
if(0) %already done
    clear all; close all;
    %-----------------------%
    paramFile = '48';
    nameFile  = '1-320';
    vecRounds = 1:320;
    %-----------------------%

    pathParameters = ['..\ParameterAnalysis\paramSets_',paramFile];

    tabHIV = [];
    tabS   = [];
    tabCt  = [];
    tabNg  = [];

    for roundNo = vecRounds
        tabHIV_rnd = readtable([pathParameters,'\round_',num2str(roundNo),'\allParametersSets_HIV']);
        tabHIV = [tabHIV;tabHIV_rnd];

        tabS_rnd = readtable([pathParameters,'\round_',num2str(roundNo),'\allParametersSets_syphilis']);
        tabS = [tabS;tabS_rnd];

        tabCt_rnd = readtable([pathParameters,'\round_',num2str(roundNo),'\allParametersSets_Ct']);
        tabCt = [tabCt;tabCt_rnd];

        tabNg_rnd = readtable([pathParameters,'\round_',num2str(roundNo),'\allParametersSets_Ng']);
        tabNg = [tabNg;tabNg_rnd];
    end

    pathAnalysis  = ['..\ParameterAnalysis\analysis_',paramFile,'\rds_',nameFile];
    writetable(tabHIV,[pathAnalysis,'\tabHIV_param'])
    writetable(tabS,  [pathAnalysis,'\tabS_param'])
    writetable(tabCt, [pathAnalysis,'\tabCt_param'])
    writetable(tabNg, [pathAnalysis,'\tabNg_param'])
end

%%
clear all; close all;
%-----------------------%
paramFile = '48';
nameFile  = '1-320';
%-----------------------%
pathAnalysis  = ['..\ParameterAnalysis\analysis_',paramFile,'\rds_',nameFile];
TABALPHA = readtable([pathAnalysis,'\tabAlpha_concatenated_rds_',nameFile]);
TABPARAM_HIV = readtable([pathAnalysis,'\tabHIV_param']);
TABPARAM_S   = readtable([pathAnalysis,'\tabS_param']);
TABPARAM_Ct  = readtable([pathAnalysis,'\tabCt_param']);
TABPARAM_Ng  = readtable([pathAnalysis,'\tabNg_param']);


%% alpha'
close all
pHIV = 20;
kit  = 'h';

alpha_k_p = readtable([pathAnalysis,'\tabAlpha_kit_',kit,'_p_',num2str(pHIV),'_elim']);

alpha_k_p.IDech_id = alpha_k_p.IDech;

alpha_k_p_params = join(alpha_k_p,TABPARAM_HIV);

% variables = {'Ptot_prep_0','Ptot_prep_base','Ptot_prev_0', 'Ptot_prev_base', 'Pun_prep_0','Pun_prep_base', ...
%     'Pund_prev_0', 'Pund_prev_base',
variables ={'Ptot_prev_base','R_prep_0','R_prep_base', 'R_prev_0', 'R_prev_base', 'alpha_p0', ...
    'betaC', 'betaI', 'rhob', 'sigma', 'theta0','zeta','mu'};

fig = figure('Renderer', 'painters', 'Position', [10 10 900 650]);
tiledlayout(4,4);
for param = variables
    nexttile

    plot(alpha_k_p_params.(param{:}), alpha_k_p_params.h,'.')
%     ax=gca;
%     ax.XAxis.FontSize  = 0.7*fontSize;
%     ax.YAxis.FontSize  = 0.7*fontSize;
%     ax.XLabel.FontSize = 0.7*fontSize;
%     ax.YLabel.FontSize = 0.7*fontSize;
%     ax.Title.FontSize  = 0.7*fontSize;
    title(param{:})
end



%% c'_hiv