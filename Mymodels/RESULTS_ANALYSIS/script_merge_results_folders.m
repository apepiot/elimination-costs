clear all
%------------%
paramNo_1 = 6; n_1 = 320;
paramNo_2 = 10; n_2 = 320;
paramNo_1_2 = 610;
f           = 10;
%------------%
foldersPath_1 = ['..\ParameterAnalysis\results_',num2str(paramNo_1),'\'];
foldersPath_2 = ['..\ParameterAnalysis\results_',num2str(paramNo_2),'\'];

foldersPath_param1 = ['..\ParameterAnalysis\paramSets_',num2str(paramNo_1),'\'];
foldersPath_param2 = ['..\ParameterAnalysis\paramSets_',num2str(paramNo_2),'\'];

foldersPath_new = ['..\ParameterAnalysis\results_',num2str(paramNo_1_2),'\'];
foldersPath_new_param = ['..\ParameterAnalysis\paramSets_',num2str(paramNo_1_2),'\'];


%parcourir les fichiers de paramNo_1 et copier dans foldersPath_new, les
%dossiers pour lesquels le kit ne contient pas s

%list_kits = {{'HIV'},{'syphilis'},{'Ct'},{'Ng'},...
%    {'HIV','syphilis'},{'HIV','Ct'},{'HIV','Ng'}, {'syphilis','Ct'},...
%    {'syphilis','Ng'},{'Ct','Ng'}, {'HIV','syphilis','Ct'},{'HIV','syphilis','Ng'},...
%    {'HIV','Ct','Ng'},{'syphilis','Ct','Ng'}, {'HIV','syphilis','Ct','Ng'}};
%[1,5,6,7,11,12,13,15] only HIV
%[-,5,-,-,11,12,-,15] HIV and syphilis (1->2,2->5,3->6,4->8)

for k=1:n_1
    if sum((mod(k,8)==[1,3,4,7]))>=1
        currName = ['_round_',num2str(k),'\'];
        copyfile([foldersPath_1,currName],[foldersPath_new,currName]);       
        
        currName_param = ['round_',num2str(k),'\'];
        copyfile([foldersPath_param1,currName_param],[foldersPath_new_param,currName_param]);       

    end 
end

warning('a changer ici')
% n_2=84;
% corresp_fold_1 = [reshape([1:4]'+([0:n_2/4]*7),1,n_2+4); reshape([2;5;6;8]+[0:n_2/4]*8,1,n_2+4)];
% corresp_fold_2 = [165:240; reshape([2;5;6;8]+[0:18]*8,1,76)+176];
n_2=88;
corresp_fold_1 = [reshape([1:4]'+([0:n_2/4]*7),1,n_2+4); reshape([2;5;6;8]+[0:n_2/4]*8,1,n_2+4)];
corresp_fold_2 = [161:240; reshape([6;8]+[0:39]*8 + 184,1,80)];
corresp_fold_3 = [241:320; reshape([2;5]+[0:39]*8 + 184,1,80)];

corresp_fold = [corresp_fold_1,corresp_fold_2,corresp_fold_3]';

%changer roundNo dans les fichiers qui viennent de paramNo_1_2.
for k=1:size(corresp_fold,1)
    
    currName_2 = ['_round_',num2str(corresp_fold(k,1)),'\'];
    
    if isfolder([foldersPath_2,currName_2])    
        %copie des dossiers results_i\_round_j
        nameFold_12 = ['_round_',num2str(corresp_fold(k,2)),'\'];
        copyfile([foldersPath_2,currName_2],[foldersPath_new,nameFold_12]);     
        
        %copie des dossiers params_i\round_j
        currName_param_2 = ['round_',num2str(corresp_fold(k,1)),'\'];
        currName_param_12 = ['round_',num2str(corresp_fold(k,2)),'\'];
        copyfile([foldersPath_param2,currName_param_2],[foldersPath_new_param,currName_param_12]);       
    end
end


%% supprimer les jeux de parametres non utilises dans paramSets_8, 9, 10
paramNo = 9;
%------------%
foldersPath = ['..\ParameterAnalysis\results_',num2str(paramNo),'\'];

S = dir(foldersPath);
noms_fichiers = {S.name};
n = nnz(~ismember({S.name},{'.','..'})&[S.isdir]);
nom_dossiers = noms_fichiers([S.isdir]);
list_num =[];
for dossier=1:n
    currDossier = nom_dossiers(dossier);
    newNum = str2num(erase(currDossier{:},'_round_'));
    list_num = [list_num,newNum];
end


num_dossiers_to_delete = setdiff(1:155,sort(list_num));
foldersPath = ['..\ParameterAnalysis\paramSets_',num2str(paramNo),'\'];

for k=num_dossiers_to_delete
    doss_to_deleted = ['round_',num2str(k),''];
    rmdir([foldersPath,doss_to_deleted], 's')
end




