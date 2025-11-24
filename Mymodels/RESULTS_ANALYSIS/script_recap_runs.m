%% Recap des fichiers qui ont tourné
%------------%
paramNo = 11;
f       = 1;
%------------%
foldersPath = ['..\ParameterAnalysis\results_',num2str(paramNo),'\'];

S = dir(foldersPath);
%n = nnz(~ismember({S.name},{'.','..'})&[S.isdir]);
n = 5;

tabRecap = table('Size',[n,5],'VariableTypes',{'double','char','double','double','char'},...
    'VariableNames',{'roundNo','kit','lastId_a','lastId_c','date'});
list_p=0:0.1:0.5;
n_p = length(list_p);
tabRecap_p = table('Size',[320*n_p,7],'VariableTypes',{'double','char','double','double','double','double','char'},...
    'VariableNames',{'roundNo','kit','p','n_by_p','lastId_a','lastId_c','date'});
j=1;
for k=1:n
    %currName = S(i).name;
    %no = str2double(erase(currName,'_round_'));
    currName = ['_round_',num2str(k)];
    %[a,b] = ;
    if isfolder([foldersPath,currName,'\'])
        A = readtable([foldersPath,currName,'\tabAlpha']);
        lastId_a = A.nbEch(end);
        
        if isfile([foldersPath,currName,'\elimCosts_f_',num2str(f),'.txt'])
            B = readtable([foldersPath,currName,'\elimCosts_f_',num2str(f),'.txt']);
            lastId_c = B.noEch(end);
            kit = B.kit(1);
        else
            kit={''};
            lastId_c = 0;
        end
        tabRecap(k,:).kit=kit;
        tabRecap(k,:).lastId_a = lastId_a;
        tabRecap(k,:).lastId_c = lastId_c;
        %tabRecap(k,:).date={S(i).date};
        
        for p=list_p
            tabRecap_p(j,:).p = p;
            tabRecap_p(j,:).kit      = kit;
            tabRecap_p(j,:).lastId_a = lastId_a;
            tabRecap_p(j,:).lastId_c = lastId_c;
            tabRecap_p(j,:).n_by_p   = sum(A.p==p);
            tabRecap_p(j,:).roundNo=k;
            j=j+1;
        end
    end
    tabRecap(k,:).roundNo=k;
end

sortrows(tabRecap,{'roundNo'},{'ascend'})


%% Somme par modele
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\RESULTS_ANALYSIS')
list_kits = unique(tabRecap.kit)';
%list_p = unique(
for kit=list_kits
    for p=list_p
        if ~isempty(kit{:})
        thiskit=isCharEqCell(tabRecap_p.kit,kit{:})' & tabRecap_p.p==p;
        disp(['kit=',kit{:},', p=',num2str(p),' : ',num2str(sum(tabRecap_p(thiskit,:).n_by_p))])
        disp(' ')
        end
    end
end

%
writematrix(tabRecap(tabRecap.lastId_c>0,:).roundNo',[foldersPath,'roundNos',num2str(paramNo),'.txt'])







