%Analyse des résultats %voir aussi result_analysis_v5_bis
clear all
f=1; %beep off
for strat=1:15
    path  = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\_sauvegarde du 23-03-22\_b_1\';
    path2 = [path,'_round_1_3','\Strat_',num2str(strat),'_b_',num2str(f),'.txt'];
    A = readtable(path2,'ReadVariableNames', true);
    tot = A;
    for round=2:4
        path2 = [path,'_round_',num2str(round),'_100','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_100','.txt'];
        A = readtable(path2,'ReadVariableNames', true);
        tot = [tot;A];
    end
    writetable(tot,[path,'Strat_',num2str(strat),'_b_',num2str(f),'_concatenated.txt'])
end

%% Creation des fichiers summary pour chaque strategie
clear all
f=1;
for strat=1:15
    path  = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\_sauvegarde du 23-03-22\_b_1\';
    tot = readtable([path,'Strat_',num2str(strat),'_b_',num2str(f),'_concatenated.txt']);
    
    %histogram(tot.x1) %Ct
    %histogram(tot.x2) %Ng
    %histogram(tot.x3) %HIV
    %histogram(tot.x4) %syph
    
    Infection =       ["Ct";               "Ng";               "HIV";              "Syph"];
    Moyenne   = round([mean(tot.x1);       mean(tot.x2);       mean(tot.x3);       mean(tot.x4)],3);
    Mediane   = round([median(tot.x1);     median(tot.x2);     median(tot.x3);     median(tot.x4)],3);
    Std       = round([std(tot.x1);        std(tot.x2);        std(tot.x3);        std(tot.x4)],3); %standard deviation (t-student)
    N         = round([length(tot.x1);     length(tot.x2);     length(tot.x3);     length(tot.x4)],3);
    Min       = round([min(tot.x1);        min(tot.x2);        min(tot.x3);        min(tot.x4)],3);
    Max       = round([max(tot.x1);        max(tot.x2);        max(tot.x3);        max(tot.x4)],3);
    low_95    = round([prctile(tot.x1,05); prctile(tot.x2,05); prctile(tot.x3,05); prctile(tot.x4,05)],3);
    high_95   = round([prctile(tot.x1,95); prctile(tot.x2,95); prctile(tot.x3,95); prctile(tot.x4,95)],3);
    low_90    = round([prctile(tot.x1,10); prctile(tot.x2,10); prctile(tot.x3,10); prctile(tot.x4,10)],3);
    high_90   = round([prctile(tot.x1,90); prctile(tot.x2,90); prctile(tot.x3,90); prctile(tot.x4,90)],3);
    
    sumup_tab = table(Infection,Moyenne,Mediane,Std,N,Min,Max,low_95,high_95,low_90,high_90);
    writetable(sumup_tab,[path,'Strat_',num2str(strat),'_b_',num2str(f),'_summary.txt'])
end


%% Comparaison des moyennes des couts (Student test)
% Distribution de reference : strategie [1], [2], [3], [4]
clear all; 
f=1;
path  = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\_sauvegarde du 23-03-22\_b_1\';
strat1 = readtable([path,'Strat_1_b_',num2str(f),'_concatenated.txt']);

recapTstudent = table('Size',[14,10],'VariableTypes',["string",repmat("double",1,9)],...
    'VariableNames',["tStudent","z_Ct","p_Ct","z_Ng","p_Ng","z_HIV","p_HIV","z_s","p_s","N_i"]);

for strat=1:15
    tot = readtable([path,'Strat_',num2str(strat),'_b_',num2str(f),'_concatenated.txt']);
    recapTstudent(strat,"tStudent") = {['strat1_vs_strat_',num2str(strat)]};
    %[H1,P1,CI1,STATS1] = ttest2(strat1.x1,tot.x1,'Vartype','unequal');
    [z1,p1] = myStudentTest(strat1.x1,tot.x1);
    recapTstudent(strat,"z_Ct") = {z1};
    recapTstudent(strat,"p_Ct") = {p1};
    
    %[H2,P2,CI2,STATS2] = ttest2(strat1.x2,tot.x2,'Vartype','unequal');
    [z2,p2] = myStudentTest(strat1.x2,tot.x2);
    recapTstudent(strat,"z_Ng") = {z2};
    recapTstudent(strat,"p_Ng") = {p2};
    
    %[H3,P3,CI3,STATS3] = ttest2(strat1.x3,tot.x3,'Vartype','unequal');
    [z3,p3] = myStudentTest(strat1.x3,tot.x3);
    recapTstudent(strat,"z_HIV") = {z3};
    recapTstudent(strat,"p_HIV") = {p3};
    
    %[H4,P4,CI4,STATS4] = ttest2(strat1.x4,tot.x4,'Vartype','unequal');
    [z4,p4] = myStudentTest(strat1.x4,tot.x4);
    recapTstudent(strat,"z_s") = {z4};
    recapTstudent(strat,"p_s") = {p4};
    recapTstudent(strat,"N_i") = {length(tot.x4)};
    
    recapTstudent{:,2:10} = round(recapTstudent{:,2:10},3);
end
% diff significative si p<0.05
writetable(recapTstudent,[path,'_tStudent_versus_strategy_1.txt'])

%% Comparaison des distributions des couts  (Mann-Whitney-Wilcoxon)
% Distribution de reference : strategie [1], [2], [3], [4]
clear all; 
f=1;
path  = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\_sauvegarde du 23-03-22\_b_1\';
strat1 = readtable([path,'Strat_1_b_',num2str(f),'_concatenated.txt']);

recapMWW = table('Size',[14,10],'VariableTypes',["string",repmat("double",1,9)],...
    'VariableNames',["MWW","z_Ct","p_Ct","z_Ng","p_Ng","z_HIV","p_HIV","z_s","p_s","N_i"]);

for strat=1:15
    tot = readtable([path,'Strat_',num2str(strat),'_b_',num2str(f),'_concatenated.txt']);
    recapMWW(strat,"MWW") = {['strat1_vs_strat_',num2str(strat)]};
    %[H1,P1,CI1,STATS1] = ttest2(strat1.x1,tot.x1,'Vartype','unequal');
    [z1,p1] = myMWWTest(strat1.x1,tot.x1);
    recapMWW(strat,"z_Ct") = {z1};
    recapMWW(strat,"p_Ct") = {p1};
    
    %[H2,P2,CI2,STATS2] = ttest2(strat1.x2,tot.x2,'Vartype','unequal');
    [z2,p2] = myMWWTest(strat1.x2,tot.x2);
    recapMWW(strat,"z_Ng") = {z2};
    recapMWW(strat,"p_Ng") = {p2};
    
    %[H3,P3,CI3,STATS3] = ttest2(strat1.x3,tot.x3,'Vartype','unequal');
    [z3,p3] = myMWWTest(strat1.x3,tot.x3);
    recapMWW(strat,"z_HIV") = {z3};
    recapMWW(strat,"p_HIV") = {p3};
    
    %[H4,P4,CI4,STATS4] = ttest2(strat1.x4,tot.x4,'Vartype','unequal');
    [z4,p4] = myMWWTest(strat1.x4,tot.x4);
    recapMWW(strat,"z_s") = {z4};
    recapMWW(strat,"p_s") = {p4};
    recapMWW(strat,"N_i") = {length(tot.x4)};
    
    recapMWW{:,2:10} = round(recapMWW{:,2:10},3);
end
% diff significative si p<0.05
writetable(recapMWW,[path,'_MWW_versus_strategy_1.txt'])


