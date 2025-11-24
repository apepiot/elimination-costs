% Results analysis of MAIN_11
close all; clear all;
%addpath('C:\Users\pepiot\Documents\PhD\codes\multi-voluntary-testing\Mymodels\MAIN')
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\MAIN')
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\INITIALISATION');
%------------------------
parametrizationFile = '11';
%------------------------
pathSims   = ['..\ParameterAnalysis\results_',parametrizationFile,'\'];
pathBackup = ['..\ParameterAnalysis\analysis_',parametrizationFile,'\'];
[status, msg, msgID] = mkdir(pathBackup);

rougeHIV = [215, 0, 0]/255;     cols.h = rougeHIV;
jauneS   = [250, 215, 0]/255;   cols.s = jauneS;
bleuCt   = [56, 57, 186]/255;   cols.c = bleuCt;
vertNg   = [43, 152, 38]/255;   cols.g = vertNg;
vecRounds = readmatrix([pathSims,'roundNos',num2str(parametrizationFile),'.txt']);

%% Copies des dossiers de results_1 dans results_1_2_mix
% on ajoute un "-" devant
%
if (0)
    vecRounds = [1,3:15,40:44,50:61];
    for roundNo=vecRounds
        patFile = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\results_1\_round_',num2str(roundNo),'\'];
        mkdir(patFile);
        pathCopy = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\results_1_2_mix\_round_-',num2str(roundNo),'\'];
        mkdir(pathCopy)
        if isfile([patFile,'elimCosts_f_1.txt'])
            elimCosts = readtable([patFile,'elimCosts_f_1.txt']);
            copyfile([patFile,'elimCosts_f_1.txt'],[pathCopy,'elimCosts_f_1_all.txt']);
        else
            elimCosts = readtable([patFile,'elimCosts.txt']);
            copyfile([patFile,'elimCosts.txt'],[pathCopy,'elimCosts_f_1_all.txt']);
        end
        elimCosts = elimCosts(contains(elimCosts.kit,'hscg'),:);
        elimCosts.roundNo = - elimCosts.roundNo;
        writetable(elimCosts,[pathCopy,'elimCosts_f_1.txt']);
        
        copyfile([patFile,'tabAlpha.txt'],[pathCopy,'tabAlpha.txt']);
        
        pathTabAlpha = [patFile,'tabAlpha.txt'];
        opts       = detectImportOptions(pathTabAlpha);
        opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
        opts.VariableTypes(contains(opts.VariableNames,'msg')) = {'char'};
        tabAlpha   = readtable(pathTabAlpha,opts);
        
        if 0 %on ne prend pas en compte les modeles reduits
            %on retire les colonnes k_mod des modeles reduits
            Abis = tabAlpha;
            %Identification du modele le plus gros calcule
            mods_elim = Abis.Properties.VariableNames(contains(Abis.Properties.VariableNames,'_elim'));
            cellsz = cellfun(@length,mods_elim,'uni',false);
            mod_elim = mods_elim(max(cell2mat(cellsz)) <= cell2mat(cellsz)+1); %s_hs, s_s
            k_mod = erase(mod_elim,'_elim'); disp(k_mod);
            %on retire les colonnes k_mod et k_mod_elim des modeles reduits
            varUnderscore = Abis.Properties.VariableNames(contains(Abis.Properties.VariableNames,'_') &...
                ~contains(Abis.Properties.VariableNames,'elim_'));
            colToReplace = varUnderscore(~contains(varUnderscore,k_mod));
            
            for i=1:length(colToReplace)
                Abisbis = removevars(Abis,{colToReplace{i}});
                %Abis(:,:).(colToReplace{i}) = NaN(size(Abis,1),1);
                Abis = Abisbis;
            end
            %         varnotUnderscore = Abis.Properties.VariableNames(~(contains(Abis.Properties.VariableNames,'_') &...
            %             ~contains(Abis.Properties.VariableNames,'elim_')));
            %         Abisbis = Abis(:,varnotUnderscore);
            
            Abis.roundNo = -Abis.roundNo;
            writetable(Abis,[pathCopy,'tabAlpha_hscg.txt']);
        end
        copyfile([patFile,'paramRho.txt'],[pathCopy,'paramRho.txt']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANALYSE DES ALPHAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Global table of alphas
close all;
%--------------------------------%
vecRounds = 1:8 %readmatrix([pathSims,'roundNos',num2str(parametrizationFile),'.txt']);
%reshape(8*[0:29] + [1:5]',1,30*5);     %[-[1,3:15,40:44,50:61],1:45,46,50,51,52,56,57];
f         = 1;
%--------------------------------%

for roundNo=vecRounds
    fileName = ['_round_',num2str(roundNo),'\tabAlpha.txt'];
    pathTabAlpha = [pathSims,fileName];
    
    %     if roundNo>0
    fileName = ['_round_',num2str(roundNo),'\tabAlpha.txt'];
    pathTabAlpha = [pathSims,fileName];
    opts       = detectImportOptions(pathTabAlpha);
    opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
    opts.VariableTypes(contains(opts.VariableNames,'msg')) = {'char'};
    opts.Delimiter = ",";
    A   = readtable(pathTabAlpha,opts);
    A.roundNo_old = A.roundNo;
    A.roundNo     = roundNo*ones(size(A,1),1);
    
    if roundNo<0
        A.roundNo = -abs(A.roundNo);
    end
    
    if roundNo~=vecRounds(1)
        B = outerjoin(B,A,'MergeKeys', true);
    else
        B=A;
    end
    %     else
    %         fileName = ['_round_',num2str(roundNo),'\tabAlpha_hscg.txt'];
    %         pathTabAlpha = [pathSims,fileName];
    %         opts       = detectImportOptions(pathTabAlpha);
    %         opts.VariableNamesLine = 1;
    %         opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
    %         opts.VariableTypes(contains(opts.VariableNames,'msg')) = {'char'};
    %         opts.VariableTypes(opts.VariableNames=="Var1") = {'double'};
    %         opts.VariableTypes(opts.VariableNames=="Var8") = {'char'};
    %         opts.VariableTypes(opts.VariableNames=="Var14") = {'char'};
    %         opts.VariableTypes(opts.VariableNames=="Var13") = {'char'};
    %         opts.Delimiter = ",";
    %         A   = readtable(pathTabAlpha,opts);
    %         writetable(A,[pathSims,'_round_',num2str(roundNo),'\tabAlpha.txt']);
    %         if roundNo~=vecRounds(1)
    %             B = outerjoin(B,A,'MergeKeys', true);
    %         else
    %             B=A;
    %         end
    %     end
end

B(:,contains(B.Properties.VariableNames,"Var"))=[];

res = arrayAsACompactedString(vecRounds);
mkdir([pathBackup,'rds_',res])
writetable(B,[pathBackup,'rds_',res,'\tabAlpha_concatenated_rds_',res,'.txt'])

%% Table de la proportion d'infections eliminees sans dépistage volontaire

warning('ici utiliser uniquement les jeux de parametres de 8')
clear all; close all;
%-----------------------%
paramFile = '11';
f         = 1;
nameFolder  = '1-8';% '1-320';
%-----------------------%
pathBackup = ['..\ParameterAnalysis\analysis_',paramFile,'\'];
TABALPHA_0 = readtable([pathBackup,'rds_',nameFolder,'\tabAlpha_concatenated_rds_',nameFolder,'.txt']);

TABALPHA = TABALPHA_0(TABALPHA_0.roundNo~=TABALPHA_0.roundNo_old,:);

vecP = 0.6:0.1:1; %0:0.1:0.5;

tabElim = zeros(4,length(vecP));
for i=1:length(vecP)
    p=vecP(i);
    disp(p)
    N = size(TABALPHA(TABALPHA.p==p,:),1);
    tabElim(1,i) = round(sum(TABALPHA(TABALPHA.p==p,:).elim_h)/N,2)*100;
    tabElim(2,i) = round(sum(TABALPHA(TABALPHA.p==p,:).elim_s)/N,2)*100;
    tabElim(3,i) = round(sum(TABALPHA(TABALPHA.p==p,:).elim_c)/N,2)*100;
    tabElim(4,i) = round(sum(TABALPHA(TABALPHA.p==p,:).elim_g)/N,2)*100;
end
disp(tabElim)



%%
clear all; close all;
%-----------------------%
paramFile = '610';
f         = 10;
nameFile  = '1-320';
%-----------------------%
pathCopy  = ['..\ParameterAnalysis\results_',paramFile,'\'];
vecRounds = readmatrix([pathCopy,'roundNos',num2str(paramFile),'.txt']);
pathWrite = ['..\ParameterAnalysis\analysis_',paramFile,'\rds_',nameFile,'\'];
%pathW_recap = ['C:\Users\pepiot\Documents\RecapPhD\tables\rds_',nameFile,'_',num2str(f),'\'];%,[fileName,'_',infi],'.tex'];
pathW_recap = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\tables\rds_',nameFile,'_',num2str(f),'\'];%,[fileName,'_',infi],'.tex'];

% tabModels = [vecRounds;...
%     [15*ones(1,length([1,3:15,40:44,50:61])),[1:15,1:15,1:15]],1,5,6,7,11,12]';

makeHist  = 0;
makeTable = 1;
vecP = 0:0.1:1;
%-----------------------------------%
list_mod_h = [1,5,6,7,11,12,13,15];
tabModels = [vecRounds;list_mod_h(mod(vecRounds-1,8)+1)]';

%TO DO : pour chaque modele, recuperer le taux d'elimination, de chaque
%infection du modele
%e.g. kit = {'h','c'}, on regarde quelle infection est eliminee en premier
%et on recupere h_h ou c_c pour l'autre c'est la valeur de rho_hc
%------------------------%
list_kits = {{'HIV'},{'syphilis'},{'Ct'},{'Ng'},...
    {'HIV','syphilis'},{'HIV','Ct'},{'HIV','Ng'}, {'syphilis','Ct'},...
    {'syphilis','Ng'},{'Ct','Ng'}, {'HIV','syphilis','Ct'},{'HIV','syphilis','Ng'},...
    {'HIV','Ct','Ng'},{'syphilis','Ct','Ng'}, {'HIV','syphilis','Ct','Ng'}};
% list_kits = {{'HIV'},{'HIV','syphilis'},{'HIV','Ct'},{'HIV','Ng'},...
%     {'HIV','syphilis','Ct'}};
for nmod=1:15
    vecRds = tabModels(tabModels(:,2)==nmod,1);
    k = indexKit(list_kits{nmod});
    for p=vecP
        B=[];
        for roundNo=vecRds'
            fileName  = [pathCopy,'_round_',num2str(roundNo),'\tabAlpha.txt'];
            opts         = detectImportOptions(fileName);
            opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
            opts.VariableTypes(contains(opts.VariableNames,'msg')) = {'char'};
            opts.Delimiter = ",";
            A    = readtable(fileName,opts);
            Abis = A(abs(A.p-p)<1e-5,:);
            B    = [B;Abis];
        end
        if ~isempty(B)
            warning('interet de k ?')
            writetable(B, [pathWrite,'tabAlpha_kit_',k,'_p_',num2str(round(p*100,0)),'.txt'])
        end
    end
end
%
for nmod=1:15
    k = indexKit(list_kits{nmod});
    disp(k)
    for p=vecP
        fileName  = [pathWrite,'tabAlpha_kit_',k,'_p_',num2str(round(p*100,0)),'.txt'];
        if exist(fileName)
            opts      = detectImportOptions(fileName);
            opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
            opts.VariableTypes(contains(opts.VariableNames,'msg')) = {'char'};
            opts.Delimiter = ",";
            tabAlpha = readtable(fileName,opts);
            
            if ~isempty(tabAlpha)
                for j=1:size(tabAlpha,1)
                    thisSet = tabAlpha(j,:);
                    nk = length(k);
                    if contains(k,'h')
                        k_mod = [k,'_',k];
                    else
                        k_mod = [k,'_h',k];
                    end
                    sub_k = k;
                    while nk>0
                        k_mod_elim = thisSet.([k_mod,'_elim']);
                        if ~isempty(k_mod_elim) && ~strcmp(k_mod_elim{:},'0') && ~isempty(k_mod_elim{:})
                            nk = nk-length(k_mod_elim{:});
                            inf_eliminated = k_mod_elim{:};
                            alpha_k_mod = thisSet.(k_mod);
                            for inf=inf_eliminated
                                k_mod = erase(k_mod,inf);
                                thisSet.(inf) = alpha_k_mod;
                                sub_k = erase(sub_k,inf);
                            end
                        elseif k_mod_elim{:}=='0'
                            for inf=k
                                if thisSet.(['elim_',inf]) %infection deja eliminee
                                    k_mod = erase(k_mod,inf);
                                    thisSet.(inf) = 0;
                                    nk = nk-1;
                                    sub_k = erase(sub_k,inf);
                                else
                                    thisSet.(inf)=double(NaN);
                                    %error('stop here')
                                end
                            end
                            nk=0;
                        else
                            error('stop here')
                        end
                    end
                    if j==1
                        recapAlpha = thisSet;
                    else
                        recapAlpha = [recapAlpha;thisSet];
                    end
                end
            else
                recapAlpha=cell2table(cell(0,length(recapAlpha.Properties.VariableNames)), 'VariableNames', recapAlpha.Properties.VariableNames);%table('Size',[0,length(recapAlpha.Properties.VariableNames)],...
                %'VariableNames',recapAlpha.Properties.VariableNames,...
                %'VariableTypes',class(recapAlpha));
            end
            writetable(recapAlpha, [pathWrite,'tabAlpha_kit_',k,'_p_',num2str(round(p*100,0)),'_elim.txt'])
        end
        disp('')
    end
end
%
rougeHIV = [215, 0, 0]/255;     cols.h = rougeHIV;
jauneS   = [250, 215, 0]/255;   cols.s = jauneS;
bleuCt   = [56, 57, 186]/255;   cols.c = bleuCt;
vertNg   = [43, 152, 38]/255;   cols.g = vertNg;
%vecP=0:0.1:0.5
%%
mkdir(pathW_recap)
vecP=0:0.1:0.5;
for nmod=1:15
    k = indexKit(list_kits{nmod});
    disp(k)
    for inf=k
        M = [];
        i=1;
        for p=vecP
            pathAlpha=[pathWrite,'tabAlpha_kit_',k,'_p_',num2str(round(p*100,0)),'_elim.txt'];
            if exist(pathAlpha)
                tabAlpha_p_k = readtable(pathAlpha);
                alphas = tabAlpha_p_k.(inf);
                if isempty(alphas)
                    alphas=NaN;
                end
                if makeTable
                    M(i,1) = p;
                    n      = length(alphas);            M(i,2) = n;
                    avg    = mean(alphas,'omitnan');    M(i,3) = avg;
                    med    = median(alphas,'omitnan');  M(i,4) = med;
                    q025   = quantile(alphas,0.025);    M(i,5) = q025;
                    q975   = quantile(alphas,0.975);    M(i,6) = q975;
                    min_a  = min(alphas);               M(i,7) = min_a;
                    max_a  = max(alphas);               M(i,8) = max_a;
                    %n_0    = sum(alphas<1e-5);          M(i,10)= round(n_0/n*100,0);
                    %tabAlpha_p_k(table2array(sum(tabAlpha_p_k(:,strcat('elim_',num2cell(k))),2))>0,:)
                    subTab = table2array(tabAlpha_p_k(:,strcat('elim_',num2cell(k))));
                    n_0    = sum(sum(subTab,2)>0);      M(i,9)= round(n_0/n*100,0);
                    n_nan = sum(isnan(alphas) &  sum(subTab,2)==0);         
                                                        M(i,10) = round(n_nan/n*100,0);

                    i      = i+1;
                    disp([k,', p=',num2str(p)])
                    disp(round([M(:,3),M(:,5),M(:,6)],2))
                    %
                    %disp([round(tabAlpha_k_inf.avg,2),round(tabAlpha_k_inf.q025,2),round(tabAlpha_k_inf.q975,2)])

                end
                if makeHist
                    fig = histogram(alphas,'FaceColor',cols.(inf),'EdgeColor','none','Normalization','probability');
                    title(['$\rho^\prime_',inf,'$ ',...
                        num2str(round(mean(alphas,'omitna'),3)),' [',...
                        num2str(round(quantile(alphas,0.025),3)),',',...
                        num2str(round(quantile(alphas,0.975),3)),']',', n=',num2str(n)], 'Interpreter','latex')
                    mkdir([pathWrite,'_histograms\_alphas\'])
                    saveas(fig,[pathWrite,'_histograms\_alphas\hist_alphas_p_',num2str(round(p*100)),'_',k,'_',inf,'_rdns_',nameFile,'.png'])
                end
                close all;
                clear alphas;
            end
        end
        
        if makeTable
            if ~isempty(M)
                % save M
                tabAlpha_k_inf = array2table(M,'VariableNames',{'p','n','avg','med','q025','q975','min','max','pct_0','pct_nan'});
                %disp(tabAlpha_k_inf)
                caption = ['$\rho_{',k,',',inf,'}^\prime $'];
                label = ['tab:alpha',k,inf];
                nameFile_a = ['recapAlpha_',k,'_',inf,'.tex'];
                mkdir([pathWrite,'_summaries\']);
                tableToTex_col_row(tabAlpha_k_inf,[pathWrite,'_summaries\',nameFile_a],caption,label,0)
                tableToTex_col_row(tabAlpha_k_inf,[pathW_recap,nameFile_a],caption,label,0)
           end
        end
        disp('')
    end
end

%% Recapitulative table of rho_k,i p in column, kit in line
pathRead = ['..\ParameterAnalysis\analysis_610\rds_1-320\'];
vecP=0:0.1:0.5;
list_k = {'h','s','c','g','hs','hc','hg','sc','sg','cg','hsc','hsg','hcg','scg','hscg'};
for i={'h','s','c','g'}
    inf = i{:};
    list_k_i = list_k(isCharInCell(list_k,'h') & isCharInCell(list_k,inf));
    recap = table('Size',[8,7],'VariableTypes',repmat({'string'},1,7));
    j=1;
    for k=list_k_i
        kit = k{:};
        recap(j,:).Var1 = kit;
        m = 2;
        for p=vecP
            recapAlpha =  readtable([pathRead,'tabAlpha_kit_',kit,'_p_',num2str(round(p*100,0)),'_elim.txt']);
            alphas = recapAlpha.(inf);
            %disp(length(alphas));
            moy = round(mean(alphas,'omitnan'),2);
            ic_low = round(quantile(alphas,0.025),2);
            ic_up  = round(quantile(alphas,0.975),2);
            compact_str = [num2str(moy),' [',num2str(ic_low),',',num2str(ic_up),']'];
            recap(j,:).(['Var',num2str(m)]) = compact_str;
            m = m+1;
        end         
        j=j+1; 
    end
    disp(recap)
end

%% Summary table of alphas
% close all;
% clearvars -except pathBackup
% %--------------------------------%
% vecRounds = [-[1,3:15,40:44,50:61],1:45];
% f=1;
% recalculerRES=0;
%
% %faire en sorte que les roundNo<1 prennent juste hscg
% %--------------------------------%
%
% list_kits = {{'HIV'},{'syphilis'},{'Ct'},{'Ng'},...
%     {'HIV','syphilis'},{'HIV','Ct'},{'HIV','Ng'}, {'syphilis','Ct'},...
%     {'syphilis','Ng'},{'Ct','Ng'}, {'HIV','syphilis','Ct'},{'HIV','syphilis','Ng'},...
%     {'HIV','Ct','Ng'},{'syphilis','Ct','Ng'}, {'HIV','syphilis','Ct','Ng'}};
%
% list_kits_to_consider = list_kits;
%
% fileRds = arrayAsACompactedString(vecRounds);
%
% filePath = [pathBackup,'rds_',fileRds,'\tabAlpha_concatenated_rds_',fileRds,'.txt'];
% opts       = detectImportOptions(filePath);
% opts.VariableTypes(contains(opts.VariableNames,'_elim')) = {'char'};
% tableAlpha_rdns   = readtable(filePath,opts);
%
% if recalculerRES
%     RES = table('Size',[0,13],'VariableNames', {'IDech','roundNo','nbEch','p','f','h','s', 'c', 'g','k_mod','kit','model','infElim'},...
%         'VariableTypes',{'double','double','double','double','double','double','double','double','double','char','char','char','char'});
%
%     for j=1:size(tableAlpha_rdns,1)
%         if ceil(j/10)==j/10
%             disp([num2str(j),'/',num2str(size(tableAlpha_rdns,1))]);
%         end
%         temp  = tableAlpha_rdns(j,:);
%         temp2 = temp(:,1:4);
%         temp2.f = f;
%         for i=1:length(list_kits_to_consider)
%             k   = indexKit(list_kits_to_consider{i});
%             mod = k;
%             res = findCnn_inf_mod(temp,k,mod);
%             res.k_mod = [k,'_',mod];
%             res.kit   = k;
%             res.model = mod;
%             res.infElim = temp.([res.k_mod,'_elim']);
%             RES = [RES;[temp2,res]];
%
%             if ~contains(k,'h')
%                 mod = ['h',k];
%                 res = findCnn_inf_mod(temp,k,mod);
%                 res.k_mod = [k,'_',mod];
%                 res.kit   = k;
%                 res.model = mod;
%                 res.infElim = temp.([res.k_mod,'_elim']);
%                 RES = [RES;[temp2,res]];
%             end
%         end
%     end
% end
%
% writetable(RES, [pathBackup,'rds_',fileRds,'\tabAlpha_concatenated_rds_',fileRds,'_all_k_mod.txt'] );
%
%
% %% Par k_mod
% %RES2 = RES(contains(RES.k_mod,'h'),:);
% %% ici il y a des problemes lies au fait qu'on prend des echantillons independants
% close all; clear all;
%
% %--------------------------------%
% vecRounds = [-[1,3:15,40:44,50:61],1:45];
% makeHists  = 0;
% makeTables = 1;
% f          = 1;
% parametrizationFile ='1_2_mix';
% %--------------------------------%
%
% nameFile = arrayAsACompactedString(vecRounds);
% pathBackup = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFile,'\'];
% RES = readtable([pathBackup,'\tabAlpha_concatenated_rds_',nameFile,'_all_k_mod.txt'] );
% namesVar = {'k_mod','p','b','n','dfe','med_h','avg_h','q025_h','q975_h','nan_h','n0_h',...
%     'med_s','avg_s','q025_s','q975_s','nan_s','n0_s',...
%     'med_c','avg_c','q025_c','q975_c','nan_c','n0_c',...
%     'med_g','avg_g','q025_g','q975_g','nan_g','n0_g',...
%     'h_elim','s_elim','c_elim','g_elim'};
% typesVar = cell(1,length(namesVar)); typesVar(:) = {'double'};
% typesVar(1) = {'char'};
% summaryTable = table('Size',[11*22,length(namesVar)],'VariableTypes',typesVar,'VariableNames',namesVar);
% i=1;
% list_k_mods = unique(RES.k_mod);
% BinWidth    = 0.01;
% for k_mod = list_k_mods'
%     for pHIV=0:0.1:1
%         disp(k_mod)
%         RES_k_mod = RES(ismember(RES.k_mod,k_mod) & abs(RES.p-pHIV)<1e-10 & RES.f==1,:);
%         n = size(RES_k_mod.h,1);
%         if makeHists
%             mkdir([pathBackup,'_histograms\_alphas\']);
%             if sum(sum(isnan(table2array(RES_k_mod(:,{'h','s','c','g'}))),1)) < 4*n
%                 fig = figure();
%                 subplot(2,2,1)
%                 histogram(RES_k_mod.h,'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
%                 title(['$\rho\prime_h$ ',num2str(round(mean(RES_k_mod.h,'omitna'),3)),' [',num2str(round(quantile(RES_k_mod.h,0.025),3)),',',num2str(round(quantile(RES_k_mod.h,0.975),3)),']'], 'Interpreter','latex')
%                 subplot(2,2,2)
%                 histogram(RES_k_mod.s,'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
%                 title(['$\rho\prime_s$ ',num2str(round(mean(RES_k_mod.s,'omitna'),3)),' [',num2str(round(quantile(RES_k_mod.s,0.025),3)),',',num2str(round(quantile(RES_k_mod.s,0.975),3)),']'], 'Interpreter','latex')
%                 subplot(2,2,3)
%                 histogram(RES_k_mod.c,'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
%                 title(['$\rho\prime_c$ ',num2str(round(mean(RES_k_mod.c,'omitna'),3)),' [',num2str(round(quantile(RES_k_mod.c,0.025),3)),',',num2str(round(quantile(RES_k_mod.c,0.975),3)),']'], 'Interpreter','latex')
%                 subplot(2,2,4)
%                 histogram(RES_k_mod.g,'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
%                 title(['$\rho\prime_g$ ',num2str(round(mean(RES_k_mod.g,'omitna'),3)),' [',num2str(round(quantile(RES_k_mod.g,0.025),3)),',',num2str(round(quantile(RES_k_mod.g,0.975),3)),']'], 'Interpreter','latex')
%
%                 sgtitle(['Kit \{',RES_k_mod.kit{1}, '\}, Model \{',RES_k_mod.model{1},'\}, p=',num2str(pHIV),...
%                     ', NaN h:',num2str(round(sum(isnan(RES_k_mod.h))/n*100)),'% ',...
%                     's:',num2str(round(sum(isnan(RES_k_mod.s))/n*100)),'%, ',...
%                     'c:',num2str(round(sum(isnan(RES_k_mod.c))/n*100)),'%, ',...
%                     'g:',num2str(round(sum(isnan(RES_k_mod.g))/n*100)),'%, '])
%
%                 saveas(fig,[pathBackup,'_histograms\_alphas\hist_alphas_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',k_mod{:},...
%                     '_rdns_',nameFile,'.png'])
%             end
%         end
%         if makeTables
%             % Big table
%             summaryTable(i,:).k_mod    = k_mod;
%             summaryTable(i,:).p        = pHIV;
%             summaryTable(i,:).b        = f;
%             summaryTable(i,:).n        = 0;
%             for inf='hscg'
%                 alphas = RES_k_mod.(inf);
%                 summaryTable(i,:).n        = max(sum(~isnan(alphas)),summaryTable(i,:).n);
%                 summaryTable(i,:).(['med_',inf])    = median(alphas,'omitna');
%                 summaryTable(i,:).(['avg_',inf])    = mean(alphas,'omitna');
%                 summaryTable(i,:).(['q025_',inf])   = quantile(alphas,0.025);
%                 summaryTable(i,:).(['q975_',inf])   = quantile(alphas,0.975);
%                 summaryTable(i,:).(['nan_',inf])    = round(sum(isnan(alphas))./summaryTable(i,:).n,2)*100;
%                 summaryTable(i,:).(['n0_',inf])     = round(sum(abs(alphas)<1e-6)./summaryTable(i,:).n,2)*100;
%                 summaryTable(i,:).([inf,'_elim'])= round(sum(ismember(RES_k_mod.infElim,inf))./summaryTable(i,:).n,2)*100;
%             end
%             summaryTable(i,:).dfe      = round(sum(isCharInCell(RES_k_mod.infElim,'0'))./summaryTable(i,:).n*100,0);
%             i = i+1;
%         end
%         close all
%     end
%    %disp('coucou')
% end
% close all
%
% if makeTables
%     mkdir([pathBackup,'_summaries\']);
%     writetable(summaryTable,[pathBackup,'_summaries\alphas_by_p_and_k_mod','_rdns_',nameFile,'.txt'])
% end
%
% % Table by k_mod
% if makeTables
%     for k_mod = list_k_mods'
%         varToKeep  = summaryTable.Properties.VariableNames;
%         [k,mod] = retrieve_K_MOD(k_mod{:});
%         for inf='hscg'
%             if ~contains(k,inf)
%                 varToKeep(contains(varToKeep,['_',inf]))=[];
%             end
%         end
%         k_modTable = summaryTable(ismember(summaryTable.k_mod,k_mod),varToKeep);
%         writetable(k_modTable,[pathBackup,'_summaries\alphas_',k_mod{:},'_rdns_',nameFile,'.txt'])
%     end
% end


%% alphas_k_mod_rdns_vecRounds.txt to .tex
%     clear all;
%     close all;
%     %--------------------------------%
%     parametrizationFile = '1_2_mix';
%     vecRounds = [-[1,3:15,40:44,50:61],1:45];
%     f = 1;
%     %--------------------------------%
%     nameRds = arrayAsACompactedString(vecRounds);
%     pathBackup = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameRds,'\_summaries\'];
%     mkdir(['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\tables\rds_',nameRds,'\']);
%
%     list_k_mod = {'h_h','s_s','s_hs','c_c','c_hc','g_g','g_hg','hs_hs','hc_hc','hg_hg','sc_sc','sc_hsc','sg_sg','sg_hsg','cg_cg','cg_hcg','hsc_hsc','hsg_hsg','scg_scg','scg_hscg','hscg_hscg'};
%     for i=1:length(list_k_mod)
%         k_mod = list_k_mod{i};
%
%         fileName = ['alphas_',k_mod,'_rdns_',nameRds];
%
%         disp(['\input{tables/rds_1-45/',fileName,'.tex}'])
%         %disp(k_mod)
%         filePath = [pathBackup,fileName,'.txt'];
%         tab = readtable(filePath);
%         tab_copy = tab(:,2:(end-4));
%
%         %Splitting by infection
%         varNames = tab_copy.Properties.VariableNames;
%         m = length(varNames);
%         ninf = (m-4)/6;
%         for j=1:ninf
%             tab_copy2 = tab_copy(:,[1:4,4+6*(j-1)+ (1:6)]);
%
%             med_inf = varNames(6*(j-1)+4+1);
%             infi     = erase(med_inf{:},'med_');
%
%             %Table of elimination
%             pathW = [pathBackup,fileName,'_',infi,'.tex'];
%             pathW_recap = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\tables\rds_',nameRds,'\',[fileName,'_',infi],'.tex'];
%             caption = ['$\rho_{',infi,',',transformCharTex(k_mod),'}^\prime $'];
%             label = ['alpha',infi,k_mod];
%             tableToTex(tab_copy2,pathW,caption,label)
%             tableToTex(tab_copy2,pathW_recap,caption,label)
%         end
%     end
%


%% Temps de calcul des alphas
% par modele
if (0)
close all;
clear all;
%--------------------------------%
vecRounds = 1:240;
f         = 3;
parametrizationFile = '5';
%--------------------------------%
nameFile = arrayAsACompactedString(vecRounds);
pathBackup = ['..\ParameterAnalysis\analysis_',parametrizationFile,'\'];
tabAlpha = readtable([pathBackup,'rds_',nameFile,'\tabAlpha_concatenated_rds_',nameFile,'.txt']);

vecP        = unique(tabAlpha.p)';
% tableTimes = table('Size',[length(vecP)*4,6],...
%     'VariableNames',{'p','n_k','t_avg','t_min','t_max','n'},...
%     'VariableTypes',{'double','double','double','double','double','double'});
tableTimes = table('Size',[length(vecP)*4,4],...
    'VariableNames',{'p','n_k','t_avg','t_tot'},...
    'VariableTypes',{'double','double','double','double'});
i=1;
list_k_mod = {'h_h','s_s','c_c','g_g',...
    'hs_hs','hc_hc','hg_hg','sc_sc','sg_sg','cg_cg',...
    'hsc_hsc','hsg_hsg','scg_scg','hscg_hscg'};
for pHIV=vecP
    for j=1:4
        ech2=[];
        for k_mod = list_k_mod
            ech1 = tabAlpha(abs(tabAlpha.p-pHIV)<1e-15 & ~isnan(tabAlpha.(k_mod{:})),:);
            if length(k_mod{:})==j*2+1
                ech2 = [ech1;ech2];
                
            end
        end
        tableTimes(i,:).p = pHIV;
        tableTimes(i,:).t_avg = round(mean(ech2.timeCompil),0);
        %tableTimes(i,:).med = round(median(ech2.timeCompil),0);
        %tableTimes(i,:).t_min = round(min(ech2.timeCompil),0);
        %tableTimes(i,:).t_max = round(max(ech2.timeCompil),0);
        %tableTimes(i,:).n   = sum(~isnan(ech2.timeCompil));
        tableTimes(i,:).n_k = j;
        
        i = i+1;
    end
    tableTimes(i-1,:).t_tot = sum(tableTimes([i-4:i-1],:).t_avg.*[4;6;4;1]);
end

for p=0:0.1:0.5
    res = tableTimes(abs(tableTimes.p-p)<1e-10,:);
    disp([num2str(round(sum(res.t_avg.*[4;6;4;1])*100/60/60,0)),' hours for p=',num2str(p), ' for 100 simus']) %100 simus
end

pathW_recap = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\tables\rds_',nameFile,'_',num2str(f),'\','times_alpha.tex'];
capt        =  'Compilation times of $\rho_{\texttt{k}}^\prime$';
if 1
    tableToTex(tableTimes,pathW_recap,capt,'tab:timesalpha',0)
end



% Estimation temps juste HIV par p
tableTimes = table('Size',[length(vecP),3],...
    'VariableNames',{'p','n_k','t_rho'},...
    'VariableTypes',{'double','double','double'});
list_k_mod = {'h_h','hs_hs','hc_hc','hg_hg','hsc_hsc','hsg_hsg','hscg_hscg'};
i=1;
for pHIV=vecP
    ech1 = tabAlpha(abs(tabAlpha.p-pHIV)<1e-15,:);
    tps_p = sum(ech1.timeCompil);
    n_p   = length(ech1.timeCompil);
    
    tableTimes(i,:).p     = pHIV;
    tableTimes(i,:).('t_rho') = round(tps_p/60/60,0);
    tableTimes(i,:).n_k   = n_p/8;
    i = i+1;
end

disp(tableTimes)
pathW_recap = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\tables\rds_',nameFile,'_',num2str(f),'\','times_alpha_hiv.tex'];
capt        =  'Computation times of $\rho_{\texttt{k}}^\prime$';
tableToTex(tableTimes,pathW_recap,capt,'tab:timesalpha_hiv',0)
writetable(tableTimes,[pathBackup,'rds_',nameFile,'/','times_alpha_hiv.txt'])
end %if (0)









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANALYSE DES COUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Summary table of cnn
close all; clear all;
%--------------------------------%
parametrizationFile = '610';
f         = 10;
%--------------------------------%
pathSims   = ['..\ParameterAnalysis\results_',parametrizationFile,'\'];
vecRounds = readmatrix([pathSims,'roundNos',num2str(parametrizationFile),'.txt']);
nameFile = arrayAsACompactedString(vecRounds);
pathBackup = ['..\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFile,'\'];

tabAlpha = readtable([pathBackup,'tabAlpha_concatenated_rds_',nameFile,'.txt']);
B=[];

for roundNo=vecRounds 
    fileName = ['_round_',num2str(roundNo),'\elimCosts_f_',num2str(f),'.txt'];
    pathTabCosts = [pathSims,fileName];
    if exist(pathTabCosts)
        opts       = detectImportOptions(pathTabCosts);
        opts.VariableTypes(contains(opts.VariableNames,'msgCnnK')) = {'char'};
        tabCosts = readtable(pathTabCosts,opts);
        
        tabCosts.roundNo_old = tabCosts.roundNo;
        tabCosts.roundNo     = ones(size(tabCosts,1),1)*roundNo;
        
        tabCosts.nbEch = tabCosts.noEch;
        %tabCosts.roundNo = abs(tabCosts.roundNo);
        tabCosts = join(tabCosts, tabAlpha(:,{'IDech','p','roundNo','elim_h','elim_s','elim_c','elim_g','nbEch'}));
        B = [B;tabCosts];
    end
end

size(B(B.HIV>=1,:))
size(B(B.HIV<=-1,:))

writetable(B,[pathBackup,'elimCosts_concatenated_rds_',nameFile,'.txt'])

%% TO DO : lancer le recalcul dans issuesIdentification2.m
warning('lancer issuesIdentification2.m')

%% Histograms of the costs by k_mod and by p
close all; clear all;

%--------------------------------%
parametrizationFile = '48';
f                   = 1;
%--------------------------------%
pathSims   = ['..\ParameterAnalysis\results_',parametrizationFile,'\'];
vecRounds = readmatrix([pathSims,'roundNos',num2str(parametrizationFile),'.txt']);

rougeHIV=[215, 0, 0]/255;
jauneS  =[250, 215, 0]/255;
bleuCt  =[56, 57, 186]/255;
vertNg  =[43, 152, 38]/255;

nameFile = arrayAsACompactedString(vecRounds);
pathBackup = ['..\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFile,'\'];

filePath = [pathBackup,'elimCosts_concatenated_rds_',nameFile,'_correction.txt'];
opts       = detectImportOptions(filePath);
tableCosts_rdns   = readtable(filePath,opts);

mkdir([pathBackup,'_histograms/_costs/']);

list_k_mods = unique(tableCosts_rdns.kit);
minC    = -0.5;
binWidth = (1-(-1))/20;
%BinWidth=0.01;
for k_mod = list_k_mods'
    if ~isempty(k_mod{:}) %&& contains(k_mod{:},'h') && length(k_mod{:})<5
        for pHIV=0:0.1:0.5
            disp(k_mod)
            RES_k_mod = tableCosts_rdns(ismember(tableCosts_rdns.kit,k_mod) & abs(tableCosts_rdns.p-pHIV)<1e-15 & tableCosts_rdns.f==f,:);
            n = size(RES_k_mod.HIV,1);
            
            costs_h = RES_k_mod.HIV;
            costs_s = RES_k_mod.syphilis;
            costs_c = RES_k_mod.Ct;
            costs_g = RES_k_mod.Ng;
            %c_h
            if sum(isnan(costs_h))~=size(RES_k_mod,1)
                fig = figure();
                histogram(costs_h,'FaceColor',rougeHIV,'EdgeColor','none','Normalization','probability')%,'BinWidth',binWidth); %HIV
                sgtitle(['Kit \{',RES_k_mod.kit{1}, '\}, Model \{',RES_k_mod.kit{1},'\}, p=',num2str(pHIV),...
                    ', NaN : ',num2str(round(sum(isnan(costs_h))/n*100)),'%, ',...
                    num2str(round(mean(costs_h,'omitnan'),3)),' [',num2str(round(prctile(costs_h,2.5),3)),',',num2str(round(prctile(costs_h,97.5),3)),']'])
                xlabel('$c_h^\prime$','Interpreter','latex')
                %xlim([minC,1])
                saveas(fig,[pathBackup,'_histograms\_costs\hist_costs_h_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',k_mod{:},...
                    '_rdns_',nameFile,'.png'])
            end
            
            %c_s
            if sum(isnan(costs_s))~=size(RES_k_mod,1)
                fig = figure();
                histogram(costs_s,'FaceColor',jauneS,'EdgeColor','none','Normalization','probability')%,'BinWidth',binWidth);
                sgtitle(['Kit \{',RES_k_mod.kit{1}, '\}, Model \{',RES_k_mod.kit{1},'\}, p=',num2str(pHIV),...
                    ', NaN : ',num2str(round(sum(isnan(costs_s))/n*100)),'%, ',...
                    num2str(round(mean(costs_s,'omitnan'),3)),' [',num2str(round(prctile(costs_s,2.5),3)),',',num2str(round(prctile(costs_s,97.5),3)),']'])
                xlabel('$c_s^\prime$','Interpreter','latex')
                %xlim([minC,1])
                saveas(fig,[pathBackup,'_histograms\_costs\hist_costs_s_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',k_mod{:},...
                    '_rdns_',nameFile,'.png'])
            end
            
            %c_c
            if sum(isnan(costs_c))~=size(RES_k_mod,1)
                fig = figure();
                histogram(costs_c,'FaceColor',bleuCt,'EdgeColor','none','Normalization','probability')%,'BinWidth',binWidth);
                sgtitle(['Kit \{',RES_k_mod.kit{1}, '\}, Model \{',RES_k_mod.kit{1},'\}, p=',num2str(pHIV),...
                    ', NaN : ',num2str(round(sum(isnan(costs_c))/n*100)),'%, ',...
                    num2str(round(mean(costs_c,'omitnan'),3)),' [',num2str(round(prctile(costs_c,2.5),3)),',',num2str(round(prctile(costs_c,97.5),3)),']'])
                xlabel('$c_c^\prime$','Interpreter','latex')
                %xlim([minC,1])
                saveas(fig,[pathBackup,'_histograms\_costs\hist_costs_c_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',k_mod{:},...
                    '_rdns_',nameFile,'.png'])
            end
            
            %c_g
            if sum(isnan(costs_g))~=size(RES_k_mod,1)
                fig = figure();
                histogram(costs_g,'FaceColor',vertNg,'EdgeColor','none','Normalization','probability')%,'BinWidth',binWidth);
                sgtitle(['Kit \{',RES_k_mod.kit{1}, '\}, Model \{',RES_k_mod.kit{1},'\}, p=',num2str(pHIV),...
                    ', NaN : ',num2str(round(sum(isnan(costs_g))/n*100)),'%, ',...
                    num2str(round(mean(costs_g,'omitnan'),3)),' [',num2str(round(prctile(costs_g,2.5),3)),',',num2str(round(prctile(costs_g,97.5),3)),']'])
                xlabel('$c_g^\prime$','Interpreter','latex')
                %xlim([minC,1])
                saveas(fig,[pathBackup,'_histograms\_costs\hist_costs_g_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',k_mod{:},...
                    '_rdns_',nameFile,'.png'])
            end
        end
        close all;
    end
    close all;
end

% Quick analysis
% Check for simulations that didn't work
%- the first infection eliminated has not been found
%- error msgCnnk diff de 0 à la fin ?
median(tableCosts_rdns.HIV,'omitnan')
outliders=tableCosts_rdns(tableCosts_rdns.HIV<-2,:)

head(tableCosts_rdns(isnan(tableCosts_rdns.HIV),:))

%head(tableCosts_rdns)

%1.put NAN in the lines when one disease has been eliminated and computing
%the model does not make sense
%-> already done in issueseidentification2
% 
% for i=1:size(tableCosts_rdns,1)
%     kit = tableCosts_rdns(i,:).kit{:};
%     for inf=kit
%         if tableCosts_rdns(i,:).(['elim_',inf])
%             for inf2='hscg'
%                 bigInf = kToKit(inf2);
%                 tableCosts_rdns(i,:).(bigInf{:}) = NaN;
%             end
%         end
%     end
% end

size(tableCosts_rdns(isnan(tableCosts_rdns.HIV) & tableCosts_rdns.p<0.3,:))



%% Summary table of costs by k_mod and by p
close all; clear all;

%--------------------------------%
parametrizationFile = '59';
f   = 3;
inf = 'HIV';
%pathRecap = 'C:\Users\pepiot\Documents\RecapPhD\';
pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\';
%--------------------------------%
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\MAIN\');

pathSims   = ['..\ParameterAnalysis\results_',parametrizationFile,'\'];
vecRounds = readmatrix([pathSims,'roundNos',num2str(parametrizationFile),'.txt']);

nameFile = arrayAsACompactedString(vecRounds);
pathBackup = ['..\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFile,'\'];

filePath = [pathBackup,'elimCosts_concatenated_rds_',nameFile,'_correction.txt'];
opts       = detectImportOptions(filePath);
tableCosts_rdns   = readtable(filePath,opts);

list_strat = {{'h','s','c','g'},{'h','s','cg'},{'h','sc','g'},{'hs','c','g'},...
    {'h','sg','c'},{'hc','s','g'},{'hg','s','c'},{'hs','cg'},{'hc','sg'},{'hg','sc'},...
    {'hsc','g'},{'hsg','c'},{'hcg','s'},{'h','scg'},{'hscg'}};
%
all_kits = unique(tableCosts_rdns.kit);
res = isCharInCell(all_kits,indexKit(inf));
list_k_mods = all_kits(res);
summaryTable_costs = table('Size',[length(list_k_mods)*6,12],...
    'VariableTypes',{'char','char','double','double','double','double','double','double','double','double','double','double'},...
    'VariableNames',{'kit','strats','p','b','n','med','avg','q25','q75','q025','q975','pct_nan'});
i=1;
for k=list_k_mods'
    disp(k)
    
    %Strategies concerned by the kit k_mod
    strats = [''];
    for j=1:length(list_strat)
        s = list_strat{j};
        if ismember(k,s)
            strats = [strats,' ',num2str(j)];
        end
    end
    
    if ~isempty(k{:})
        for pHIV=0:0.1:0.5
            RES_k_mod = tableCosts_rdns(ismember(tableCosts_rdns.kit,k) & abs(tableCosts_rdns.p-pHIV)<1e-15 & tableCosts_rdns.f==f,:);
            n = size(RES_k_mod,1);
            costs_i = RES_k_mod.(inf);
            summaryTable_costs(i,:).med    = median(costs_i,'omitna');
            summaryTable_costs(i,:).avg    = mean(costs_i,'omitna');
            summaryTable_costs(i,:).q25    = quantile(costs_i,0.25);
            summaryTable_costs(i,:).q75    = quantile(costs_i,0.75);
            summaryTable_costs(i,:).q025   = quantile(costs_i,0.025);
            summaryTable_costs(i,:).q975   = quantile(costs_i,0.975);
            summaryTable_costs(i,:).kit    = k;
            summaryTable_costs(i,:).pct_nan = round(sum(isnan(costs_i))/n*100,0);
            summaryTable_costs(i,:).p      = pHIV;
            summaryTable_costs(i,:).b      = f;
            summaryTable_costs(i,:).n      = n;
            summaryTable_costs(i,:).strats = {strats};
            i=i+1;
        end
    end
end
writetable(summaryTable_costs,[pathBackup,'_summaries\costs_',inf,'_b_',num2str(f),...
    '_rdns_',nameFile,'.txt'])


inf   = 'HIV';
fileName = ['costs_',inf,'_b_',num2str(f),'_rdns_',nameFile];
% path0 = [pathBackup,'_summaries\',fileName,'.txt'];
pathW = [pathBackup,'_summaries\',fileName,'.tex'];
% 
% tab   = readtable(path0,ReadVariableNames=true);

tab = summaryTable_costs;
tableToTex_col_row(tab,pathW,['Summary table of costs $c_',indexKit(inf),'$'],'tab:costs_mod_h',1)
pathW_recap = [pathRecap,'tables\rds_',nameFile,'_',num2str(f),'\',fileName,'.tex'];
tableToTex_col_row(tab,pathW_recap,['Summary table of costs $c_',indexKit(inf),'$ for b=',num2str(f)],['tab:costs_mod_h_b_',num2str(f)],1)


%% Recapitulative table of c_k,i p in column, kit in line
close all; clear all;
%--------------------------------%
parametrizationFile = '610';
f   = 10;
inf = 'HIV';
%pathRecap = 'C:\Users\pepiot\Documents\RecapPhD\';
pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\';
%--------------------------------%
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\MAIN\');

pathSims   = ['..\ParameterAnalysis\results_',parametrizationFile,'\'];
vecRounds = readmatrix([pathSims,'roundNos',num2str(parametrizationFile),'.txt']);
nameFile = arrayAsACompactedString(vecRounds);
pathBackup = ['..\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFile,'\'];
filePath = [pathBackup,'elimCosts_concatenated_rds_',nameFile,'_correction.txt'];
opts       = detectImportOptions(filePath);
tableCosts_rdns   = readtable(filePath,opts);

tableCosts_rdns(tableCosts_rdns.HIV==Inf,:).HIV = [NaN;NaN];

vecP=0:0.1:0.5;
list_k = {'h','s','c','g','hs','hc','hg','sc','sg','cg','hsc','hsg','hcg','scg','hscg'};
for i={'h'}%,'s','c','g'}
    inf = i{:};
    list_k_i = list_k(isCharInCell(list_k,'h') & isCharInCell(list_k,inf));
    recap = table('Size',[8,7],'VariableTypes',repmat({'string'},1,7));
    j=1;
    for k=list_k_i
        kit = k{:};
        recap(j,:).Var1 = kit;
        m = 2;
        costs_k = tableCosts_rdns(isCharEqCell(tableCosts_rdns.kit,kit),:);
        for p=vecP
            costs_k_p = costs_k(costs_k.p==p ,:).HIV;
            %disp(length(costs_k_p))
            moy = round(mean(costs_k_p,'omitnan'),2);
            ic_low = round(quantile(costs_k_p,0.025),2);
            ic_up  = round(quantile(costs_k_p,0.975),2);
            compact_str = [num2str(moy),' [',num2str(ic_low),',',num2str(ic_up),']'];
            recap(j,:).(['Var',num2str(m)]) = compact_str;
            m = m+1;
        end        
        j=j+1;
    end
    disp(recap)
end


%% Histograms par kit
%comparaison des modeles contenant HIV 2 à 2, une figure par valeur de p
close all
%-------------------------------------------------------------
%myFilter = @ (data) (data(data>=-1 & ~isnan(data) & data<1));
%myFilter = @(data) (data(data<1000));
myFilter = @(data) (data(data>=quantile(data,0.02) & data<=quantile(data,0.99)));
parametrizationFile = '48';
f = 1;
myTest = @kstest2;
vecP = 0:0.1:0.6;
%-------------------------------------------------------------
%pathRecap = 'C:\Users\pepiot\Documents\RecapPhD\';
pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\';
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels');
path3      = [pathRecap,'graphes\rds_',nameFile,'_',num2str(f),'\'];

list_k_mods = {'h';'s';'c';'g';'hs';'hc';'hg';'sc';'sg';'cg';'hsc';'hsg';'hcg';'scg';'hscg'};
list_k_mods = list_k_mods(isCharInCell(list_k_mods,'h'));
nMod = size(list_k_mods,1);

minBinWidth=1;
path2 = ['..\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFile,'\_histograms\_costs\'];

%to define x-axis histograms :
% minHIV=[];maxHIV=[];
% minS=[];maxS=[];
% nStrat=15;

fex=figure();
Lf = 1000; %fex.Position(3); %largeur figure
Ls = Lf/(nMod+6); %largeur subplot
a  = (Lf - nMod*Ls)/(nMod+1+1); %espace horizontal entre les subplots;

hf = 750;%fex.Position(4);  %hauteur figure
hb = 20; % espace adidtionnel entre le haut de la figure et la permière ligne
hs = (hf-hb)/(nMod+10); %hauteur subplot (10)
b  = (hf - nMod*hs-hb)/(nMod+2); %(10) %espace vertical entre les subplots;
close all

%if ~isempty(k{:})
for pHIV=vecP
    for kitA=1:nMod
        KITA = list_k_mods(kitA);
        RES_k_mod_A = tableCosts_rdns(ismember(tableCosts_rdns.kit,KITA) & abs(tableCosts_rdns.p-pHIV)<1e-15 & tableCosts_rdns.f==f,:);
        
        echA = RES_k_mod_A.HIV;
        echA_f = myFilter(echA);
        for kitB=kitA:nMod
            plotPos = [round(a+(kitB-1)*(a+Ls),3),round(-b + hf-b-hs-(kitA-1)*(hs+b)- hb,3), Ls hs];
            KITB = list_k_mods(kitB);
            if kitA==kitB
                fig3 = figure(3);
                fig3.Position = [10 10 Lf hf];
                h3 = axes('Position', plotPos);
                binWidth = max([(max(echA_f)-min(echA_f))/15,1e-3]);
                
                if ~(isempty(echA) || sum((isnan(echA)))>=length((echA))-1)
                    h=histogram(echA_f,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',binWidth);
                    set(h3,'Unit','pixels','Position',plotPos);
                    if ~isempty(echA_f)
                        n=2;
                        dec_min=0; dec_max=dec_min;
                        while dec_max<=dec_min
                            dec_min = ceil(min(echA_f)*10^n)/10^n;
                            dec_max = floor(max(echA_f)*10^n)/10^n;
                            n=n+1;
                        end
                        if abs(dec_max-dec_min)<5e-4
                            dec_max = ceil(min(echA)*10^2)/10^2;
                            dec_min = floor(max(echA)*10^2)/10^2;
                            xlim(sort([dec_min,dec_max]))
                        end
                        %xlim(sort([dec_min,dec_max]))
                        %xticks([dec_min,dec_max])
                        %xticks([])
                        ylim([0,max(h.Values)])

                    end
                    %axis off;
                    set(gca,'Fontsize',8)
                    set(gca,'xtick',[])

                    %juste mettre l'axe des abscisses
                    set(h3, 'box', 'off')
                    set(gca, 'box', 'off', 'ytick', [], 'xcolor', [0.3,0.3,0.3], 'ycolor', [1,1,1], 'color', 'none');
                end
            else
                RES_k_mod_B = tableCosts_rdns(ismember(tableCosts_rdns.kit,KITB) & abs(tableCosts_rdns.p-pHIV)<1e-15 & tableCosts_rdns.f==f,:);
                echB = RES_k_mod_B.HIV;
                echB_f = myFilter(echB);
                echAB_f = [echA_f;echB_f];

                if ~(isempty(echAB_f) || sum(isnan(echAB_f))==length(echAB_f))
                    figure(3);
                    h7 = axes('Position',plotPos);

                    %p-value
                    if sum(isnan(echA))~=length(echA) && sum(isnan(echB))~=length(echB)
                        [z,pv] = myTest(echA,RES_k_mod_B.HIV);
                        str_p = sprintf('%.2f',pv);
                    else
                        str_p = '-';
                    end
                    binWidth = max((max(echAB_f)-min(echAB_f))/15,1e-3);
                    
                    histogram(echA_f,'EdgeColor','none','Normalization','probability','BinWidth',binWidth);
                    hold on
                    histogram(echB_f,'EdgeColor','none','Normalization','probability','BinWidth',binWidth);
                    set(h7,'Unit','pixels','Position',plotPos);
                    
                    dec_min=0; dec_max=dec_min;
                    while dec_max<=dec_min
                        dec_min = ceil(min(echAB_f)*10^n)/10^n;
                        dec_max = floor(max(echAB_f)*10^n)/10^n;                            
                        n       = n+1;
                    end
                    if abs(dec_max-dec_min)<5e-4
                        dec_max = ceil(min(echAB_f)*10^2)/10^2;
                        dec_min = floor(max(echAB_f)*10^2)/10^2;
                        xlim(sort([dec_min,dec_max]))
                    end
                     
                    set(gca,'xtick',[])
                    set(gca,'ytick',[])

                    %juste mettre l'axe des abscisses
                    set(h7, 'box', 'off')
                    set(gca, 'box', 'off', 'ytick', [], 'xcolor', [0.3,0.3,0.3], 'ycolor', [1,1,1], 'color', 'none');

                    %ajout de la p-value
                    title(str_p,'FontSize',13,'Interpreter','latex')
                    %axis off;
                    %legend(['Kit ',KITA{:}],['Kit ',KITB{:}])
                    %saveas(h7,[path2,'hist_HIV_kits_',KITA{:},'_vs_',...
                    %    KITB{:},'_b_',num2str(f),'_p_',num2str(100*pHIV),'.png'])
                end
            end
            
        end
    end
    clear k
    % adding names of the kits
    for numFig=3
        for nk=1:nMod
            kit = list_k_mods(nk);
            figure(numFig)
            x1 = (nk*a+(nk-1)*Ls)/Lf;
            y1 = 1-(hf-1.8*b-hb-hs-(nMod-1)*(hs+b))/hf;
%             annotation('textbox', [Ct y1 Ls/Lf hs/hf], 'string', ['$c_{',kit{:},',h}^\prime$'],...
%                 'VerticalAlignment','middle','HorizontalAlignment','center',...
%                 'LineStyle','none','FontSize',15,'Interpreter','latex')
             annotation('textbox', [x1 y1 Ls/Lf hs/hf], 'string', ['$',kit{:},'$'],...
                'VerticalAlignment','middle','HorizontalAlignment','center',...
                'LineStyle','none','FontSize',18,'Interpreter','latex')
            x2 = 1-1.5*a/Lf;
            y2 = (-b + hf-b-hs-(nk-1)*(hs+b)- hb)/hf;%
            %y2 = 1-(4.5*b+b/3+nk*b+(nk-1)*hs)/hf;
%             annotation('textbox', [Ng y2 Ls/Lf hs/hf], 'string', ['$c_{',kit{:},',h}^\prime$'],...
%                 'VerticalAlignment','middle','HorizontalAlignment','center',...
%                 'LineStyle','none','FontSize',15,'Interpreter','latex')
            annotation('textbox', [x2 y2 Ls/Lf hs/hf], 'string', ['$',kit{:},'$'],...
                'VerticalAlignment','middle','HorizontalAlignment','center',...
                'LineStyle','none','FontSize',18,'Interpreter','latex')
        end
    end
    
    annotation('line', [0. 1], 1-(hf-1.5*b-hb-hs-(nMod-1)*(hs+b))/hf*[1 1], 'Color', 'black', 'LineWidth', 1);
    annotation('line', 1-1.5*a/Lf*[1 1], [0. 1], 'Color', 'black', 'LineWidth', 1);
    
    exportgraphics(gcf,[path2,'hist_HIV_kits_b_',num2str(f),'_p_',num2str(100*pHIV),'.png'],'Resolution',300)
    exportgraphics(gcf,[path2,'hist_HIV_kits_b_',num2str(f),'_p_',num2str(100*pHIV),'.pdf'],'Resolution',300)
    %exportgraphics(gcf,[path3,'hist_HIV_kits_b_',num2str(f),'_p_',num2str(100*pHIV),'.png'],'Resolution',300)
    close all;
end
%end

%%
close all;
clearvars -except pathBackup

%--------------------------------%
parametrizationFile = '610';
f         = 10;
%myFilter = @ (data) (data(data>=-0.8 & ~isnan(data) & data<1));
myFilter = @ (data) (data(data<1000));
%--------------------------------%
pathSims   = ['..\ParameterAnalysis\results_',parametrizationFile,'\'];
vecRounds = readmatrix([pathSims,'roundNos',num2str(parametrizationFile),'.txt']);
nameFile   = arrayAsACompactedString(vecRounds);
pathBackup = ['..\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFile,'\'];
filePath   = [pathBackup,'elimCosts_concatenated_rds_',nameFile,'_correction.txt'];
opts       = detectImportOptions(filePath);
tableCosts_rdns   = readtable(filePath,opts);

list_k_mods = unique(tableCosts_rdns.kit);
list_k_mods = list_k_mods(isCharInCell(list_k_mods,'h'));

%BinWidth=0.01;
for k_mod = list_k_mods'
    if ~isempty(k_mod{:}) %&& contains(k_mod{:},'h') && length(k_mod{:})<5
        vecP = 0:0.1:0.6;
        RES_k_mod = tableCosts_rdns(ismember(tableCosts_rdns.kit,k_mod) & tableCosts_rdns.f==f,:);
        c_min = min(myFilter(RES_k_mod.HIV));
        c_max = max(myFilter(RES_k_mod.HIV));
        binW = abs((c_max+c_min))./15;
        for pHIV=vecP
            RES_k_mod_p = RES_k_mod((RES_k_mod.p-pHIV)<1e-15,:);
            disp(k_mod)
            
            fig = figure(1);
            histogram(myFilter(RES_k_mod_p.HIV),'FaceAlpha',0.2,'EdgeColor','auto','Normalization','probability','BinWidth',binW); %HIV
            %title('$c^\prime_h$', 'Interpreter','latex')
            hold on
        end
        n = size(RES_k_mod.HIV,1);
        legend("p="+vecP)
        sgtitle(['Kit \{',RES_k_mod.kit{1}, '\}, Model \{',RES_k_mod.kit{1},'\}',...
            ', NaN: h:',num2str(round(sum(isnan(RES_k_mod.HIV))/n*100)),'% '])
        
        saveas(fig,[pathBackup,'_histograms\_costs\hist_costs_h_b_',num2str(f),'_',k_mod{:},...
            '_rdns_',num2str(vecRounds(1)),'-',num2str(vecRounds(end)),'.png'])
        close all
    end
end


%% Comparaison entre les modèles/kits
%------------------------------------
inf  = 'HIV';
vecP = 0:0.1:0.5;
f    = 1;
parametrizationFile = '48';
myTest = @kstest2;     %@myStudentTest, @ttest2, @myMWWTest, @kstest2;
%------------------------------------
pathSims   = ['..\ParameterAnalysis\results_',parametrizationFile,'\'];
vecRounds = readmatrix([pathSims,'roundNos',num2str(parametrizationFile),'.txt']);
nameFolder = arrayAsACompactedString(vecRounds);
pathBackup = ['..\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFolder,'\'];
%fileRecap  = 'C:\Users\pepiot\Documents\RecapPhD\';
fileRecap  = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\';
tableCosts_rdns = readtable([pathBackup,'elimCosts_concatenated_rds_',nameFolder,'_correction.txt']);

list_k_mods = unique(tableCosts_rdns.kit);
list_k_mods = list_k_mods(isCharInCell(list_k_mods,indexKit(inf)));

n = length(list_k_mods);
tableTest = nan(n,n);

for pHIV=vecP
    disp(pHIV)
    for i=1:n
        k1 = list_k_mods{i};
        ech1 = tableCosts_rdns(ismember(tableCosts_rdns.kit,k1) & abs(tableCosts_rdns.p-pHIV)<1e-15 & tableCosts_rdns.f==f,:).(inf);
        for j=i:n
            k2 = list_k_mods{j};
            ech2 = tableCosts_rdns(ismember(tableCosts_rdns.kit,k2) & abs(tableCosts_rdns.p-pHIV)<1e-15 & tableCosts_rdns.f==f,:).(inf);
            
            % t-test
            if sum(isnan(ech1))~=length(ech1) && sum(isnan(ech2))~=length(ech2)
                [z,pv] = myTest(ech1,ech2);
            else
                pv = NaN;
            end
            tableTest(i,j) = pv;
        end
    end  
    %m=8;
    %disp(tableTest>10^(-m-1) & tableTest<10^(-m))
    %disp(tableTest)
    fileName = ['costs_',indexKit(inf),'_comparison_kits_',func2str(myTest),'_b_',num2str(f),'_p',num2str(pHIV*100)];
    writematrix(tableTest,[pathBackup,fileName,'mat.txt']);   
    
    % To LaTeX
    simplifedTab = round(tableTest,3);
    newTable     = table(char(list_k_mods{:}));
    newTable(:,2:(n+1)) = array2table(simplifedTab);
    newTable.Properties.VariableNames = {'k',list_k_mods{:}};
    
    % Table with strings, to add p-value < 10^-m
    vartypes    = cell(1, size(newTable,2));
    vartypes(:) = {'string'};
    table_str = table('Size',size(newTable),'VariableTypes',vartypes,'VariableNames',newTable.Properties.VariableNames);
    table_str{:,1} = string(newTable{:,1});
    for j=2:size(newTable,2)
        table_str{:,j} = string(simplifedTab(:,j-1));
    end
    
    for i=1:n
        for j=1:n
            if tableTest(i,j)<10^(-15)
                table_str{i,j+1} = "{\footnotesize $<$}1e-15";
            end
            if tableTest(i,j)<10^(-4) & tableTest(i,j)>=10^(-15)
                table_str{i,j+1} = string(sprintf('%.e ',tableTest(i,j)));
            end
%             for m = 15:4
%                 if newTable{i,j+1}<10^(-m+1) && newTable{i,j+1}>=10^(-m)
%                     table_str{i,j+1} = join(["$<10^{-",m,"}$"]);
%                 end
%             end
        end
    end
    disp(table_str)
    pathW        = [pathBackup,fileName,'.txt'];
    pathW_recap  = [fileRecap,'tables\rds_',nameFolder,'_',num2str(f),'\',fileName,'.tex'];
    capt         = ['p-values of the two-sample KS test performed on $c_',indexKit(inf),...
        '$ values by two-by-two kits for $p=',num2str(pHIV),'$ and $b=',num2str(f),'$. '];
    tableToTex_damier(newTable,pathW,capt,['tab:costs_kit_p',num2str(pHIV),'_b',num2str(f)])
    tableToTex_damier(newTable,pathW_recap,capt,['tab:costs_kit_p',num2str(pHIV),'_b',num2str(f)])
    
    pathW_recap  = [fileRecap,'tables\rds_',nameFolder,'_',num2str(f),'\',fileName,'str.tex'];
    tableToTex_damier(table_str,pathW_recap,capt,['tab:costs_kit_p',num2str(pHIV),'_b',num2str(f)])

end



%% KS test by strategy (read from models)

clear all;
%----------------------------
inf  = 'HIV';
vecP = 0:0.1:0.5;
f    = 10;
parametrizationFile = '610';
myTest = @kstest2; %@myStudentTest, @ttest2, @myMWWTest, @kstest2;
%----------------------------
%pathRecap = 'C:\Users\pepiot\Documents\RecapPhD\';
pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\';
pathSims  = ['..\ParameterAnalysis\results_',parametrizationFile,'\'];
vecRounds = readmatrix([pathSims,'roundNos',num2str(parametrizationFile),'.txt']);

list_strat = {{'h','s','c','g'},{'h','s','cg'},{'h','sc','g'},{'hs','c','g'},...
    {'h','sg','c'},{'hc','s','g'},{'hg','s','c'},{'hs','cg'},{'hc','sg'},{'hg','sc'},...
    {'hsg','g'},{'hsg','c'},{'hcg','s'},{'h','scg'},{'hscg'}};

nameFolder = arrayAsACompactedString(vecRounds);
pathBackup = ['..\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFolder,'\'];
tableCosts_rdns = readtable([pathBackup,'elimCosts_concatenated_rds_',nameFolder,'_correction.txt']);

n = length(list_strat);
tableTest = nan(n,n);

for pHIV=vecP
    for i=1:n
        s1 = list_strat{i};
        k1 = s1{isCharInCell(s1,indexKit(inf))};
        ech1 = tableCosts_rdns(ismember(tableCosts_rdns.kit,k1) & abs(tableCosts_rdns.p-pHIV)<1e-15 & tableCosts_rdns.f==f,:).(inf);
        for j=i:n
            s2 = list_strat{j};
            k2 = s2{isCharInCell(s2,indexKit(inf))};
            ech2 = tableCosts_rdns(ismember(tableCosts_rdns.kit,k2) & abs(tableCosts_rdns.p-pHIV)<1e-15 & tableCosts_rdns.f==f,:).(inf);
            
            % t-test
            if sum(isnan(ech1))~=length(ech1) && sum(isnan(ech2))~=length(ech2)
                [z,pv] = myTest(ech1,ech2);
            else
                pv = NaN;
            end
            tableTest(i,j) = pv;
        end
    end
    
    fileName = ['costs_',indexKit(inf),'_comparison_strats_',func2str(myTest),'_b_',num2str(f),'_p',num2str(pHIV*100)];
    writematrix(tableTest,[pathBackup,fileName,'.txt']);
%     disp(['b=',num2str(f),' p=',num2str(pHIV)])
%     disp(tableTest)

    tableTest = round(tableTest,2);

    % To LaTeX
    simplifedTab = round(tableTest,3);
    newTable     = table(['\texttt{s}_{'+string([1:n]')+'}']);
    newTable(:,2:(n+1)) = array2table(simplifedTab);
    newTable.Properties.VariableNames = ['s','\texttt{s}_{'+string(1:n)+'}'];

    % With strings, to add p-value < 10^-m
    vartypes    = cell(1, size(newTable,2));
    vartypes(:) = {'string'};
    table_str = table('Size',size(newTable),'VariableTypes',vartypes,'VariableNames',newTable.Properties.VariableNames);
    table_str{:,1} = string(newTable{:,1});
    for j=2:size(newTable,2)
        table_str{:,j} = string(simplifedTab(:,j-1));
    end   
    for i=1:n
        for j=1:n
            if tableTest(i,j)<10^(-15)
                %table_str{i,j+1} = "{\footnotesize $<$}1e-15";
                table_str{i,j+1} = "$0^*$";
            end
            if tableTest(i,j)<10^(-4) & tableTest(i,j)>=10^(-15)
                table_str{i,j+1} = string(sprintf('%.e ',tableTest(i,j)));
            end
        end
    end  
    
    pathW        = [pathBackup,fileName,'.txt'];
    pathW_recap  = [pathRecap,'\tables\rds_',nameFolder,'_',num2str(f),'\',fileName,'.tex'];
    capt         = ['p-values of the two-sample K-S test performed on $c_',indexKit(inf),...
        '$ values by two-by-two strategies $\texttt{s}_i$ for $p=',num2str(pHIV),'$ and $b=',num2str(f),'$. $0^*$ means p-value $<$ 1e-15'];
    tableToTex_damier(newTable,pathW,capt,['tab:costs_s_p',num2str(pHIV),'_b',num2str(f)])
    tableToTex_damier(newTable,pathW_recap,capt,['tab:costs_s_p',num2str(pHIV),'_b',num2str(f)])
    pathW_recap  = [pathRecap,'\tables\rds_',nameFolder,'_',num2str(f),'\',fileName,'str.tex'];
    tableToTex_damier(table_str,pathW_recap,capt,['tab:costs_s_p',num2str(pHIV),'_b',num2str(f)])

end



%% Histograms of costs strat A vs strat B

clear all; close all;
%---------%
inf='HIV';
vecP = 0:0.1:0.5;
vecF = 10;
parametrizationFile = '610';
%myFilter = @ (data) (data(data>=-0.5 & ~isnan(data) & data<1));
%myFilter = @ (data) (data(data<1000));
myFilter = @(data) (data(data>quantile(data,0.01) & data<quantile(data,0.99)));
myFilterBnds = @(data) (data(data>quantile(data,0.005) & data<quantile(data,0.995)));

%---------%

pathSims   = ['..\ParameterAnalysis\results_',parametrizationFile,'\'];
vecRounds = readmatrix([pathSims,'roundNos',num2str(parametrizationFile),'.txt']);
%pathRecap = 'C:\Users\pepiot\Documents\RecapPhD\';
pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\';
list_strat = {{'h','s','c','g'},{'h','s','cg'},{'h','sc','g'},{'hs','c','g'},...
    {'h','sg','c'},{'hc','s','g'},{'hg','s','c'},{'hs','cg'},{'hc','sg'},{'hg','sc'},...
    {'hsg','g'},{'hsg','c'},{'hcg','s'},{'h','scg'},{'hscg'}};
nameFolder = arrayAsACompactedString(vecRounds);
pathBackup = ['..\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFolder,'\'];
tableCosts_rdns = readtable([pathBackup,'elimCosts_concatenated_rds_',nameFolder,'_correction.txt']);
for p=vecP
    for f=vecF
        mkdir([pathRecap,'\graphes\rds_',nameFolder,'_',num2str(f)]);
        close all
        %maxNBins=15;
        minBinWidth=0.05;
        figure(1);
        
        %to define x-axis histograms :
        minCt=[];maxCt=[];
        minNg=[];maxNg=[];
        minHIV=[];maxHIV=[];
        minS=[];maxS=[];
        nStrat=length(list_strat);
        
        fex=figure();
        Lf = 1000;%fex.Position(3); %largeur figure
        Ls = Lf/(nStrat+7); %largeur subplot
        a  = (Lf - nStrat*Ls)/(nStrat+1+1); %espace horizontal entre les subplots;
        
        hf = 700;%fex.Position(4);  %hauteur figure
        hs = hf/(nStrat+8); %hauteur subplot (10)
        b  = (hf - nStrat*hs)/(nStrat+3); %(10) %espace vertical entre les subplots;
        close all
        
        for stratA=1:nStrat
            sA = list_strat{stratA};
            kA = sA{isCharInCell(sA,indexKit(inf))};
            echA = tableCosts_rdns(ismember(tableCosts_rdns.kit,kA) & abs(tableCosts_rdns.p-p)<1e-15 & tableCosts_rdns.f==f,:);
            for stratB=stratA:nStrat
                plotPos = [round(a+(stratB-1)*(a+Ls),3),round(-2*b + hf-b-hs-(stratA-1)*(hs+b),3), Ls hs];
                if stratA==stratB
                    %if(0)
                    %                 fig1 = figure(1);
                    %                 fig1.Position = [10 10 Lf hf];
                    %                 %h1 = subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    %                 h1 = axes('Position', plotPos);
                    %                 %nBins = min(length([min(stratAfile.Ct):0.002:max(stratAfile.Ct)]),maxNBins);
                    %                 BinWidth=min([minBinWidth,(max(stratAfile.Ct)-min(stratAfile.Ct))./10]);
                    %                 histogram(stratAfile.Ct,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
                    %                 %minCt = min(minCt,min(stratAfile.Ct));
                    %                 %maxCt = max(maxCt,max(stratAfile.Ct));
                    %                 %set(gca,'xtick',[])
                    %                 set(gca,'ytick',[])
                    %                 set(gca,'fontsize',4)
                    %                 set(h1,'Unit','pixels','Position',plotPos);
                    %                 %hold on
                    %                 %plot([0,0],[0;100],'k')
                    %
                    %                 fig2 = figure(2);
                    %                 fig2.Position = [10 10 Lf hf];
                    %                 %h2=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    %                 h2 = axes('Position', plotPos);
                    %                 %nBins = min(length([min(stratAfile.Ng):0.002:max(stratAfile.Ng)]),maxNBins);
                    %                 BinWidth=min([minBinWidth,(max(stratAfile.Ng)-min(stratAfile.Ng))./10]);
                    %                 histogram(stratAfile.Ng,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
                    %                 minNg = min(minNg,min(stratAfile.Ng));
                    %                 maxNg = max(maxNg,max(stratAfile.Ng));
                    %                 %set(gca,'xtick',[])
                    %                 set(h2,'Unit','pixels','Position',plotPos);
                    %                 set(gca,'ytick',[])
                    %                 set(gca,'fontsize',4)
                    %                 %end
                    
                    echA_h = myFilter(echA.HIV);
                    if ~isempty(echA_h)
                        fig3 = figure(3);
                        fig3.Position = [10 10 Lf hf];
                        %h3 = subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                        h3 = axes('Position', plotPos);
                        %nBins = min(length([min(stratAfile.HIV):0.002:max(stratAfile.HIV)]),maxNBins);
                        BinWidth=max(min([minBinWidth,(max(myFilter(echA_h))-min(myFilter(echA_h)))./12]),1e-3);
                        histogram(echA_h,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
%                         minHIV = min(minHIV,min(echA_h));
%                         maxHIV = max(maxHIV,max(echA_h));
                        dec_min=0;dec_max=dec_min;
                        n=2;
                        while dec_max<=dec_min
                            dec_min = ceil(min(myFilter(echA_h))*10^n)/10^n;
                            dec_max = floor(max(myFilter(echA_h))*10^n)/10^n;
                            n=n+1;
                        end
                        xticks([dec_min,dec_max])
                        xlim([min(myFilterBnds(echA_h)),max(myFilterBnds(echA_h))])
                        %set(gca,'xtick',[])
                        set(h3,'Unit','pixels','Position',plotPos);
                        set(gca,'ytick',[])
                        set(gca,'fontsize',8)
                    end
                    
                    %                 fig4 = figure(4);
                    %                 fig4.Position = [10 10 Lf hf];
                    %                 %h4=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    %                 h4=axes('Position', plotPos);
                    %                 %nBins = min(length([min(stratAfile.syph):0.002:max(stratAfile.syph)]),maxNBins);
                    %                 BinWidth=min([minBinWidth,(max(stratAfile.syph)-min(stratAfile.syph))./10]);
                    %                 histogram(stratAfile.syph,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
                    %                 minS = min(minS,min(stratAfile.syph));
                    %                 maxS = max(maxS,max(stratAfile.syph));
                    %                 %set(gca,'xtick',[])
                    %                 set(h4,'Unit','pixels','Position',plotPos);
                    %                 set(gca,'ytick',[])
                    %                 set(gca,'fontsize',4)
                    
                else
                    sB = list_strat{stratB};
                    kB = sB{isCharInCell(sB,indexKit(inf))};
                    echB = tableCosts_rdns(ismember(tableCosts_rdns.kit,kB) & abs(tableCosts_rdns.p-p)<1e-15 & tableCosts_rdns.f==f,:);
                    
                    %Ct
                    %                 figure(1);
                    %                 %h5=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    %                 h5 = axes('Position', plotPos);
                    %                 %nBins = min(length([min(stratAfile.Ct):0.002:max(stratAfile.Ct)]),maxNBins);
                    %                 %BinWidth=max([minBinWidth,min([(max(stratBfile.Ct)-min(stratBfile.Ct))./10,(max(stratAfile.Ct)-min(stratAfile.Ct))./10])]);
                    %                 BinWidth= (max([stratBfile.Ct;stratAfile.Ct])-min([stratAfile.Ct;stratBfile.Ct]))/12;
                    %                 histogram(h5,stratAfile.Ct,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
                    %                 hold on
                    %                 %nBins = min(length([min(stratBfile.Ct):0.002:max(stratBfile.Ct)]),maxNBins,'Normalization','probability');
                    %                 histogram(h5,stratBfile.Ct,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth)
                    %                 %axis off
                    %                 set(gca,'xtick',[])
                    %                 set(gca,'ytick',[])
                    %                 set(h5,'Unit','pixels','Position',plotPos);
                    %                 fprintf(['x=',num2str(a+(stratB-1)*(a+Ls)), ' y=',num2str(hf-b-hs-(stratA-1)*(hs+b)),'\n'])
                    %                 %h=gca; h.XAxis.TickLength = [0 0];
                    %                 %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                    %                 %saveas(hCt,[path2,'hist_Ct_Strat_',num2str(stratA),'(',num2str(length(stratAfile.Ct)),')_vs_',...
                    %                 %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.Ct)),')',...
                    %                 %    '_b_',num2str(f),'.png'])
                    
                    %                 figure(2);
                    %                 %h6=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    %                 h6 = axes('Position',plotPos);
                    %                 %nBins = min(length([min(stratAfile.Ng):0.002:max(stratAfile.Ng)]),maxNBins);
                    %                 %BinWidth=min([minBinWidth,(max(stratBfile.Ng)-min(stratBfile.Ng))./10,(max(stratAfile.Ng)-min(stratAfile.Ng))./10]);
                    %                 BinWidth= (max([stratBfile.Ng;stratAfile.Ng])-min([stratAfile.Ng;stratBfile.Ng]))/12;
                    %                 histogram(stratAfile.Ng,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
                    %                 hold on
                    %                 %nBins = min(length([min(stratBfile.Ng):0.002:max(stratBfile.Ng)]),maxNBins);
                    %                 histogram(stratBfile.Ng,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth);
                    %                 %axis off
                    %                 set(h6,'Unit','pixels','Position',plotPos);
                    %                 set(gca,'xtick',[])
                    %                 set(gca,'ytick',[])
                    %                 %set(gca,'xtick',[])
                    %                 %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                    %                 %saveas(hNg,[path2,'hist_Ng_Strat_',num2str(stratA),'(',num2str(length(stratAfile.Ng)),')_vs_',...
                    %                 %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.Ng)),')',...
                    %                 %    '_b_',num2str(f),'.png'])
                    %
                    
                    echB_h = myFilter(echB.HIV);
                    if ~isempty(echA_h) && ~isempty(echB_h)
                        figure(3);
                        %h7=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                        h7=axes('Position', plotPos);
                        %BinWidth=min([minBinWidth,(max(stratBfile.HIV)-min(stratBfile.HIV))./10,(max(stratAfile.HIV)-min(stratAfile.HIV))./10]);
                        BinWidth= max((max([myFilter(echB_h);myFilter(echA_h)])-min([myFilter(echA_h);myFilter(echB_h)]))/12,1e-3);
                        %nBins = min(length([min(stratAfile.HIV):0.002:max(stratAfile.HIV)]),maxNBins);
                        histogram(h7,echA_h,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
                        hold on
                        %nBins = min(length([min(stratBfile.HIV):0.002:max(stratBfile.HIV)]),maxNBins);
                        histogram(h7,echB_h,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth);
                        %axis off
                        set(gca,'XTick',[])
                        set(gca,'ytick',[])
                        set(h7,'Unit','pixels','Position',plotPos);
                        %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                        %saveas(hHIV,[path2,'hist_HIV_Strat_',num2str(stratA),'(',num2str(length(stratAfile.HIV)),')_vs_',...
                        %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.HIV)),')',...
                        %    '_b_',num2str(f),'.png'])
                    end
                    %                 figure(4);
                    %                 %h8=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    %                 h8 = axes('Position', plotPos);
                    %                 %BinWidth=min([minBinWidth,(max(stratBfile.syph)-min(stratBfile.syph))./10,(max(stratAfile.syph)-min(stratAfile.syph))./10]);
                    %                 BinWidth=(max([stratBfile.syph;stratAfile.syph])-min([stratAfile.syph;stratBfile.syph]))/12;
                    %                 %nBins = min(length([min(stratAfile.syph):0.002:max(stratAfile.syph)]),maxNBins);
                    %                 histogram(stratAfile.syph,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
                    %                 hold on
                    %                 %nBins = min(length([min(stratBfile.syph):0.002:max(stratBfile.syph)]),maxNBins);
                    %                 histogram(stratBfile.syph,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth);
                    %                 %axis off
                    %                 set(gca,'xtick',[])
                    %                 set(gca,'ytick',[])
                    %                 set(h8,'Unit','pixels','Position',plotPos);
                    %                 %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                    %                 %saveas(hS,[path2,'hist_S_Strat_',num2str(stratA),'(',num2str(length(stratAfile.syph)),')_vs_',...
                    %                 %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.syph)),')',...
                    %                 %    '_b_',num2str(f),'.png'])
                    
                    %close all;
                    %pause(1)
                    
                end
            end
        end
        
        %add numbers of strategies
        for k=1:length(list_strat)
            for numFig=3 %1:4
                figure(numFig)
                Ct = (k*a+(k-1)*Ls )/Lf;
                y1 = 1-(hf-b-hs-14*(hs+b))/hf;
                annotation('textbox', [Ct y1 Ls/Lf hs/hf], 'string', ['$\texttt{s}_{',num2str(k),'}$'],...
                    'VerticalAlignment','middle','HorizontalAlignment','center',...
                    'LineStyle','none','Interpreter','latex','FontSize',16)
                Ng = 1-2*a/Lf;
                y2 = 1-(4*b+b/3+k*b+(k-1)*hs)/hf;
                annotation('textbox', [Ng y2 Ls/Lf hs/hf], 'string', ['$\texttt{s}_{',num2str(k),'}$'],...
                    'VerticalAlignment','middle','HorizontalAlignment','center',...
                    'LineStyle','none','Interpreter','latex','FontSize',16)
            end
        end
        
        %         %Ct
        %         saveas(fig1,[path2,'hist_Ct_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.png'])
        %         set(fig1,'PaperOrientation','landscape');
        %         saveas(fig1,[pathRecap,'hist_Ct_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.pdf'])
        %
        %         %Ng
        %         saveas(fig2,[path2,'hist_Ng_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.png'])
        %         set(fig2,'PaperOrientation','landscape');
        %         saveas(fig2,[pathRecap,'hist_Ng_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.pdf'])
        %
        %HIV
%         saveas(fig3,[pathBackup,'_histograms\_costs\hist_costs_strats_h_b_',num2str(f),'_p_',num2str(round(p*100)),...
%             '_rdns_',num2str(vecRounds(1)),'-',num2str(vecRounds(end)),'.png'])
        set(fig3,'PaperOrientation','landscape');
%         saveas(fig3,[pathBackup,'_histograms\_costs\hist_costs_strats_h_b_',num2str(f),'_p_',num2str(round(p*100)),...
%             '_rdns_',num2str(vecRounds(1)),'-',num2str(vecRounds(end)),'.pdf'])
%         
%         saveas(fig3,[pathRecap,'graphes\rds_',nameFolder,'_',num2str(f),'\hist_costs_strats_h_b_',num2str(f),'_p_',num2str(round(p*100)),...
%             '_rdns_',num2str(vecRounds(1)),'-',num2str(vecRounds(end)),'.pdf'])
%         
       exportgraphics(fig3,[pathRecap,'graphes\rds_',nameFolder,'_',num2str(f),'\hist_costs_strats_h_b_',num2str(f),'_p_',num2str(round(p*100)),...
            '_rdns_',num2str(vecRounds(1)),'-',num2str(vecRounds(end)),'.pdf'],'Resolution',300)

        %         %Syphilis
        %         saveas(fig4,[path2,'hist_S_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.png'])
        %         set(fig4,'PaperOrientation','landscape');
        %         saveas(fig4,[pathRecap,'hist_S_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.pdf'])
    end
end

close all;


%% Temps de calcul des couts
% par modele
close all;
clear all;
%--------------------------------%
vecRounds = 1:240;
f         = 1;
parametrizationFile = '4';
%myFilter = @ (data) (data(data>=-0.8 & ~isnan(data) & data<1));
%--------------------------------%
nameFile = arrayAsACompactedString(vecRounds);
pathBackup = ['C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFile,'\'];
filePath = [pathBackup,'elimCosts_concatenated_rds_',nameFile,'_correction.txt'];
tableCosts_rdns = readtable([pathBackup,'elimCosts_concatenated_rds_',nameFile,'_correction.txt']);


list_k_mods = unique(tableCosts_rdns.kit);
vecP        = unique(tableCosts_rdns.p)';
n = length(list_k_mods);
tableTimes = table('Size',[n*length(vecP),7],...
    'VariableNames',{'k_mod','p','n','avg','med','min','max'},...
    'VariableTypes',{'char','double','double','double','double','double','double'});
j=1;
for pHIV=vecP
    for i=1:n
        k1 = list_k_mods{i};
        ech1 = tableCosts_rdns(ismember(tableCosts_rdns.kit,k1) & abs(tableCosts_rdns.p-pHIV)<1e-15 & tableCosts_rdns.f==f,:);
        tableTimes(j,1) = list_k_mods(i);
        tableTimes(j,:).avg = round(mean(ech1.timeCompil),2);
        tableTimes(j,:).med = round(median(ech1.timeCompil),2);
        tableTimes(j,:).min = round(min(ech1.timeCompil),2);
        tableTimes(j,:).max = round(max(ech1.timeCompil),2);
        tableTimes(j,:).p = pHIV;
        tableTimes(j,:).n = length(ech1.IDech);
        j = j+1;
    end
end

%sum(tableTimes.avg)

%Regroupement par taille du modèle
func = @(A) (length(A{:}));
tableTimes.nmod = arrayfun(func,tableTimes.k_mod);
M=[]; i=1;
for p=0:0.1:0.5
    for n=1:4
        disp(['p=',num2str(p), ', n=',num2str(n)])
        res = tableTimes(tableTimes.nmod==n & (tableTimes.p-p)<1e-10,:);
        
        disp(sum(res.n.*res.avg)/(sum(res.n)))
        disp(' ')
        M(i,1) = n;
        M(i,2) = p;
        M(i,3) = sum(res.n.*res.avg)/(sum(res.n));
        M(i,4) = 0;
        i=i+1;
    end
    M(i-1,4) = round(M(i-4,3)*4+M(i-3,3)*6+M(i-2,3)*4+M(i-1,3),0);
end

timesTable = array2table(M, 'VariableNames',{'n_k','p','t_avg_c','t_tot_c'});
pathW_recap = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\','tables\rds_',nameFile,'_',num2str(f),'\','times_costs.tex'];
%pathW_recap = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\tables\rds_',nameFile,'_',num2str(f),'\','times_costs.tex'];
capt        =  'Compilation times of $c_{\texttt{k}}^\prime$';
if 1
    tableToTex(timesTable,pathW_recap,capt,'tab:timescosts',0)
end


% Only HIV models
tableTimes = table('Size',[length(vecP),3],...
    'VariableNames',{'p','n_k','t_cost'},...
    'VariableTypes',{'double','double','double'});
list_k_mod = {'h_h','hs_hs','hc_hc','hg_hg','hsc_hsc','hsg_hsg','hscg_hscg'};
i=1;
for pHIV=vecP
    ech1 = tableCosts_rdns(abs(tableCosts_rdns.p-pHIV)<1e-15,:);
    tps_p = sum(min(ech1.timeCompil,500));
    n_p   = length(ech1.timeCompil);
    
    tableTimes(i,:).p     = pHIV;
    tableTimes(i,:).('t_cost') = round(tps_p/60/60,0);
    tableTimes(i,:).n_k   = n_p/8;
    i = i+1;
end

disp(tableTimes)
pathW_recap = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\tables\rds_',...
    nameFile,'_',num2str(f),'\','times_cost_hiv.tex'];
capt        =  'Compilation times of $c_{\texttt{k},h}^\prime$';
tableToTex(tableTimes,pathW_recap,capt,'tab:timescost_hiv',0)


%Merge with time alpha
timesAlpha = readtable([pathBackup,'times_alpha_hiv.txt']);
mergedTabletimes = join(timesAlpha,tableTimes);
mergedTabletimes.('t_all') = mergedTabletimes.('t_rho') + mergedTabletimes.('t_cost');
disp(mergedTabletimes)
pathW_recap = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\tables\rds_',...
    nameFile,'_',num2str(f),'\','times_runs_hiv.tex'];
capt        =  'Compilation times of $\rho_{\texttt{k}}^\prime$ and $c_{\texttt{k},h}^\prime$';
tableToTex(mergedTabletimes,pathW_recap,capt,'tab:times_hiv',0)




%% Regroupement de strategies
clear all;
close all;
%--------------------------------%
parametrizationFile = '48';
f        = 1;
%myFilter = @ (data) (data(data>=-0.8 & ~isnan(data) & data<=1));
myFilter = @ (data) (data(data<1000));
myTest   = @kstest2;
vecP = 0:0.1:0.5;
vecP=0.1;
inf  = 'HIV';
thresh = 0.05;
%--------------------------------%
pathSims   = ['.\ParameterAnalysis\results_',parametrizationFile,'\'];
vecRounds  = readmatrix([pathSims,'roundNos',num2str(parametrizationFile),'.txt']);
nameFile   = arrayAsACompactedString(vecRounds);
pathBackup = ['.\ParameterAnalysis\analysis_',parametrizationFile,'\rds_',nameFile,'\'];
filePath   = [pathBackup,'elimCosts_concatenated_rds_',nameFile];%,'_correction.txt'];
tableCosts_rdns = readtable([pathBackup,'elimCosts_concatenated_rds_',nameFile,'_correction.txt']);
%pathRecap  = 'C:\Users\pepiot\Documents\RecapPhD\';
pathRecap  = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\';
pathW_recap = [pathRecap,'tables\rds_',nameFile,'_',num2str(f),'\','groups_mod_b_',num2str(f),'.tex'];
list_strat = {{'h','s','c','g'},{'h','s','cg'},{'h','sc','g'},{'hs','c','g'},...
    {'h','sg','c'},{'hc','s','g'},{'hg','s','c'},{'hs','cg'},{'hc','sg'},{'hg','sc'},...
    {'hsc','g'},{'hsg','c'},{'hcg','s'},{'h','scg'},{'hscg'}};

list_k_mods = unique(tableCosts_rdns.kit);
list_k_mods = list_k_mods(isCharInCell(list_k_mods,indexKit(inf)));

tableRecapStrat = table('Size',[6*8,8],'VariableTypes',{'char','double','double','double','double','double','double','cell'},...
    'VariableNames',{'group','p','avg','ic_low','ic_high','n','nan','strat'});
k=0;
for p=vecP
    disp(['p=',num2str(p)])
    fileName = ['costs_',indexKit(inf),'_comparison_kits_',func2str(myTest),'_b_',num2str(f),'_p',num2str(p*100)];
    res = readmatrix([pathBackup,fileName,'mat.txt']);
    
    %tabRes = {};
    
    %1. On regarde les correspondances par modèle
    for i=1:size(res,1)
        modSameAsi = [''];
        ech = [];
        mods_alike = [''];
        list_strat_i = [];
        for j=1:size(res,1)
            res(j,i) = res(i,j);
            if res(i,j) > thresh
                %modSameAsi = [modSameAsi,num2str(j),' '];
                newEch = tableCosts_rdns(ismember(tableCosts_rdns.kit,list_k_mods{j}) &...
                    abs(tableCosts_rdns.p-p)<1e-15 & tableCosts_rdns.f==f,:).('HIV');
                ech = [ech;newEch];
                disp(list_k_mods{j});
                mods_alike = [mods_alike,list_k_mods{j},','];
                
                list_strat_i = unique([list_strat_i,findModInStrats(list_k_mods{j},list_strat)]);
                
            end
        end
        mods_alike = mods_alike(1:end-1);
        k=k+1;
        ech(ech==Inf,:) = NaN;
        disp([num2str(sum(~isnan(ech))),' ',num2str(round(mean(ech,'omitna'),3)),' [', num2str(round(prctile(ech,2.5),3)),',',num2str(round(prctile(ech,97.5),3)),']'])
        disp(['  '])
        tableRecapStrat(k,:).p      = p;
        tableRecapStrat(k,:).avg    = round(mean(ech,'omitna'),3);
        tableRecapStrat(k,:).ic_low = round(prctile(ech,2.5),3);
        tableRecapStrat(k,:).ic_high= round(prctile(ech,97.5),3);
        tableRecapStrat(k,:).nan    = sum(isnan(ech));
        tableRecapStrat(k,:).group  = {mods_alike};
        tableRecapStrat(k,:).n      = sum(~isnan(ech));
        tableRecapStrat(k,:).strat  = {num2str(list_strat_i)};
        %tabRes{i} = modSameAsi;
    end  
    disp(' ')
end

% 2. On enlève les resultats redondants
[a,b] = unique(tableRecapStrat(:,2:7),'rows');
T = tableRecapStrat(b,:);
pp=cellfun(@isempty,T.group);
T = T(~pp,:);
caption = ['Groups of model for b=',num2str(f)];
tableToTex(T,pathW_recap,caption,['tab:groups_b_',num2str(f)],0)

%3. Faire un tableau plus propre
t = table('Size',[size(T,1),5], 'VariableTypes',{'char','cell','double','double','char'},...
    'VariableNames',{'group of kits','group of strategies','p','n','mean [95\% CI]'});

for i=1:size(T,1)
    t(i,:).('group of kits')        = {['(',T(i,:).group{:},')']};
    t(i,:).('group of strategies')  = T(i,:).strat;
    t(i,:).p                        = T(i,:).p;
    t(i,:).n                        = T(i,:).n;
    inter = [num2str(T(i,:).avg),' [',num2str(T(i,:).ic_low),',',num2str(T(i,:).ic_high),']'];
    t(i,:).('mean [95\% CI]') = {inter};
end
caption = ['Groups of models for b=',num2str(f)];
tableToTex_col_row(t,pathW_recap,caption,['tab:groups_b_',num2str(f)],0)


%4. un tableau par p
for p=vecP
    tp = t(t.p==p,[1,2,4,5]);
    pathW_recap = [pathRecap,'\tables\rds_',nameFile,'_',num2str(f),'\','groups_mod_b_',num2str(f),'_',num2str(100*p),'.tex'];
    caption = ['Groups of models for b=',num2str(f),', p=',num2str(p)];
    tableToTex_col_row(tp,pathW_recap,caption,['tab:groups_',num2str(100*p),'_b_',num2str(f)],0)
end










%% Comparaison des couts par b
clear all; close all;
inf ='HIV';
vecP = 0:0.1:0.5;
filter = @(data) data(data>-0.6 & data<1000);
list_k = {'h','hs','hc','hg','hsc','hsg','hcg','hscg'};
fontSize = 24;
set(groot,'defaultAxesTickLabelInterpreter','latex');  


path1   = '..\ParameterAnalysis\analysis_48\rds_1-320\';
costs1  = readtable([path1,'elimCosts_concatenated_rds_1-320_correction.txt']);
path3   = '..\ParameterAnalysis\analysis_59\rds_1-320\';
costs3  = readtable([path3,'elimCosts_concatenated_rds_1-320_correction.txt']);
path10  = '..\ParameterAnalysis\analysis_610\rds_1-320\';
costs10 = readtable([path10,'elimCosts_concatenated_rds_1-320_correction.txt']);


%pathRecap  = 'C:\Users\pepiot\Documents\RecapPhD\graphes\comparaison_b\';
pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\comparaison_b\';
pathComp = ['..\ParameterAnalysis\comparaison_b\'];

%1. histograms par kit 
n = length(list_k);

for p=vecP
    j=1;
    %figure()
    fig = figure('Renderer', 'painters', 'Position', [10 10 900 650])
    tiledlayout(3,3);
    for i=1:n
        k = list_k{i};
        ech1  = costs1(ismember(costs1.kit,k)   & abs(costs1.p-p)<1e-15 & costs1.f==1,:).(inf);
        ech3  = costs3(ismember(costs3.kit,k)   & abs(costs3.p-p)<1e-15 & costs3.f==3,:).(inf);
        ech10 = costs10(ismember(costs10.kit,k) & abs(costs10.p-p)<1e-15 & costs10.f==10,:).(inf);
        
        ech1_bis = filter(ech1(~isnan(ech1)));
        ech3_bis = filter(ech3(~isnan(ech3)));
        ech10_bis = filter(ech10(~isnan(ech10)));
        e1 = (quantile(ech1_bis,0.99)-quantile(ech1_bis,0.05))/10;
        e3 = (quantile(ech3_bis,0.99)-quantile(ech3_bis,0.05))/10;
        e10= (quantile(ech10_bis,0.99)-quantile(ech10_bis,0.05))/10;
                
        binW = max([e1,e3,e10]);

        if ~isnan(binW) & ~isempty(binW)
            if binW>0
            %subplot(3,3,i)
            nexttile
            histogram(ech1_bis,'BinWidth',binW,'EdgeColor','none','Normalization','probability','FaceAlpha',0.5)
            hold on 
            histogram(ech3_bis,'BinWidth',binW,'EdgeColor','none','Normalization','probability','FaceAlpha',0.5)
            histogram(ech10_bis,'BinWidth',binW,'EdgeColor','none','Normalization','probability','FaceAlpha',0.5)
            
            echtot = [ech1_bis;ech3_bis;ech10_bis];%
            
            low_bnd = quantile(echtot,0.01);
            up_bnd  = quantile(echtot,0.995);
            xlim([low_bnd,up_bnd])
            xticks(round([(0.9*low_bnd+0.1*up_bnd),(low_bnd+up_bnd)/2,(0.1*low_bnd+0.9*up_bnd)],2))
            %h.Annotation.LegendInformation.IconDisplayStyle = 'on'; % Ensure the icon is displayed
            end
            ax=gca;
            ax.XAxis.FontSize  = 0.7*fontSize;
            ax.YAxis.FontSize  = 0.7*fontSize;
            ax.XLabel.FontSize = 0.7*fontSize;
            ax.YLabel.FontSize = 0.7*fontSize;
            ax.Title.FontSize  = 0.7*fontSize;
            ax.YAxis.TickLength = [0,0];
            %set(gca,'xtick',[])
            set(gca,'ytick',[])
            ylabel('Frequency','Interpreter','latex')
            xlabel('$c_h^\prime$','Interpreter','latex')
            %juste mettre l'axe des abscisses
            %set(h7, 'box', 'off')
            set(gca, 'box', 'off', 'ytick', [], 'xcolor', [0.3,0.3,0.3], 'ycolor', [0.3,0.3,0.3], 'color', 'none');
           	j=j+1;

        end
        title(['$',k,'$'],'Interpreter','latex')  
    end
    lgd=legend(' $~b=1$',' $~b=3$',' $~b=10$','Location','NorthWest','FontSize',17,'Interpreter','latex','Box','off');
    lgd.Layout.Tile = j;
    x1=25 ; x2=25;
    lgd.ItemTokenSize = [x1,x2];
    %lgd.Layout.Tile = 'east';
    sgtitle(['$c_h^\prime$, $p=',num2str(p),'$'],'Interpreter','latex','FontSize',20)

   warning('pas de sauvegarde')
    if(0)
        exportgraphics(fig,[pathComp,'hist_HIV_b_1-10_p_',num2str(round(p*100)),'.pdf'],'Resolution',300)
        exportgraphics(fig,[pathRecap,'hist_HIV_b_1-10_p_',num2str(round(p*100)),'.pdf'],'Resolution',300)
    end
    close all;
end



%% Comparaison des couts par p

clear all; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');  

inf = 'HIV';
vecP = 0:0.1:0.5;
vecB = 1%[1,3,10];
fontSize = 24;
filter = @(data) data(data>-0.5 & data<1000);
list_k = {'h','hs','hc','hg','hsc','hsg','hcg','hscg'};

colors = [[0.6350 0.0780 0.1840];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];...
    [0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.4940 0.1840 0.5560]];

path1   = '..\ParameterAnalysis\analysis_48\rds_1-320\';
costs1  = readtable([path1,'elimCosts_concatenated_rds_1-320_correction.txt']);
path3   = '..\ParameterAnalysis\analysis_59\rds_1-320\';
costs3  = readtable([path3,'elimCosts_concatenated_rds_1-320_correction.txt']);
path10  = '..\ParameterAnalysis\analysis_610\rds_1-320\';
costs10 = readtable([path10,'elimCosts_concatenated_rds_1-320_correction.txt']);

cost_tot = [costs1;costs3;costs10];

%pathRecap  = 'C:\Users\pepiot\Documents\RecapPhD\graphes\comparaison_b\';
pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\comparaison_b\';
pathComp  = ['..\ParameterAnalysis\comparaison_b\'];

%1. histograms par kit
n = length(list_k);
for b=vecB
    fig = figure('Renderer', 'painters', 'Position', [10 10 900 650]);
    
    tiledlayout(3,3);
    ech2  = cost_tot(cost_tot.f==b & cost_tot.(inf)<1000,:);
    for i=1:n
        k = list_k{i};
        ech3 = ech2(ismember(ech2.kit,k),:);
        %binW = (max(filter(ech3.(inf)))-min(filter(ech3.(inf))))/50;
        binW = (quantile(filter(ech3.(inf)),0.995) - quantile(filter(ech3.(inf)),0.005))/15;
        ax = nexttile;
        %hA = axes(fig);
        for j=1:length(vecP)
            p=vecP(j);
            ech = ech3(abs(ech3.p-p)<1e-15 ,:).(inf);
            ech_bis = filter(ech(~isnan(ech)));
            if ~isempty(binW) && ~isnan(binW)
                if binW>0
                    histogram(ech_bis,'BinWidth',binW,'EdgeColor','none','FaceAlpha',0.75,'FaceColor',colors(j,:))
                    hold on
                end
            end
        end
        %box(ax,'off')

        low_bnd = quantile(filter(ech3.(inf)),0.005);
        up_bnd = quantile(filter(ech3.(inf)),0.995);
        xlim([low_bnd,up_bnd])
        xticks(unique(round([0.95*low_bnd+0.05*up_bnd, (low_bnd+up_bnd)/2, 0.05*low_bnd+0.95*up_bnd],3)))
        title(['$',k,'$'],'Interpreter','latex','FontSize',fontSize)
        
        ax=gca;
        ax.XAxis.FontSize  = 0.5*fontSize;
        ax.YAxis.FontSize  = 0.5*fontSize;
        ax.XLabel.FontSize = 0.5*fontSize;
        ax.YLabel.FontSize = 0.5*fontSize;
        ax.Title.FontSize  = 0.7*fontSize;
        ax.YAxis.TickLength = [0,0];
        set(gca,'ytick',[])
        set(gca, 'box', 'off', 'ytick', [], 'xcolor', [0.3,0.3,0.3], 'ycolor', [0.3,0.3,0.3], 'color', 'none');
        ylabel('Frequency','Interpreter','latex')
        if i==1
            lgd=legend(' $~p=0$',' $~p=0.1$',' $~p=0.2$',' $~p=0.3$',' $~p=0.4$',' $~p=0.5$','Location','NorthWest',...
                'FontSize',0.6*fontSize,'Interpreter','latex','Box','off');
            lgd.Layout.Tile = 9;
            x1=25 ; x2=25;
            lgd.ItemTokenSize = [x1,x2];
        end

    end
    
    %lgd.Layout.Tile = 'east';
    sgtitle(['$c_h^\prime$, $b=',num2str(b),'$'],'Interpreter','latex','FontSize',ax.Title.FontSize)
   
    if 0
    exportgraphics(fig,[pathComp,'hist_HIV_p_0-1_b_',num2str(b),'.pdf'],'Resolution',300)
    exportgraphics(fig,[pathRecap,'hist_HIV_p_0-1_b_',num2str(b),'.pdf'],'Resolution',300)
    end
    %close all;
end



%% Comparaison des couts dans un groupe de kits par b
clear all; close all;
inf ='HIV';
vecP = 0:0.1:0.5;
filter = @(data) data(data>-0.6 & data<1000);
list_k = {'h','hc','hg','hcg'};
%list_k = {'hs','hsc','hsg','hscg'};
fontSize = 24;
set(groot,'defaultAxesTickLabelInterpreter','latex');  


path1   = '..\ParameterAnalysis\analysis_48\rds_1-320\';
costs1  = readtable([path1,'elimCosts_concatenated_rds_1-320_correction.txt']);
path3   = '..\ParameterAnalysis\analysis_59\rds_1-320\';
costs3  = readtable([path3,'elimCosts_concatenated_rds_1-320_correction.txt']);
path10  = '..\ParameterAnalysis\analysis_610\rds_1-320\';
costs10 = readtable([path10,'elimCosts_concatenated_rds_1-320_correction.txt']);


%pathRecap  = 'C:\Users\pepiot\Documents\RecapPhD\graphes\comparaison_b\';
pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\comparaison_b\';
pathComp = ['..\ParameterAnalysis\comparaison_b\'];

n = length(list_k);
fig = figure('Renderer', 'painters', 'Position', [10 10 1250 500])
tiledlayout(2,4);
ech1=[]; ech3=[]; ech10=[];
ntile=1;
for p=vecP
    %figure()
    for i=1:n
        k = list_k{i};
        ech1_k  = costs1(ismember(costs1.kit,k)   & abs(costs1.p-p)<1e-15 & costs1.f==1,:).(inf);
        ech1 = [ech1;ech1_k];
        ech3_k  = costs3(ismember(costs3.kit,k)   & abs(costs3.p-p)<1e-15 & costs3.f==3,:).(inf);
        ech3 = [ech3,ech3_k];
        ech10_k = costs10(ismember(costs10.kit,k) & abs(costs10.p-p)<1e-15 & costs10.f==10,:).(inf);
        ech10 = [ech10;ech10_k];
    end
    
    ech1_bis = filter(ech1(~isnan(ech1)));
    ech3_bis = filter(ech3(~isnan(ech3)));
    ech10_bis = filter(ech10(~isnan(ech10)));
    e1 = (max(ech1_bis)-min(ech1_bis))/40;
    e3 = (max(ech3_bis)-min(ech3_bis))/40;
    e10= (max(ech10_bis)-min(ech10_bis))/40;
    binW = max([e1,e3,e10]);
    
    if ntile==4
        nexttile(ntile+1,[1 1])
        ntile=ntile+2;
    else
        nexttile(ntile,[1 1])
        ntile=ntile+1;
    end
    if ~isnan(binW) & ~isempty(binW)
        if binW>0
            
            echtot = [ech1_bis;ech3_bis;ech10_bis];%

            low_bnd = -0.25; % quantile(echtot,0.004);
            up_bnd  = -0.0; % quantile(echtot,0.999);%max(echtot);
            binW = (up_bnd-low_bnd)./30;
            
            histogram(ech1_bis,'BinWidth',binW,'EdgeColor','none','Normalization','probability','FaceAlpha',0.5)
            hold on 
            histogram(ech3_bis,'BinWidth',binW,'EdgeColor','none','Normalization','probability','FaceAlpha',0.5)
            histogram(ech10_bis,'BinWidth',binW,'EdgeColor','none','Normalization','probability','FaceAlpha',0.5)

            xlim([low_bnd,up_bnd])
            xticks(round([(0.93*low_bnd+0.03*up_bnd),(low_bnd+up_bnd)/2,(0.03*low_bnd+0.97*up_bnd)],2))
            %h.Annotation.LegendInformation.IconDisplayStyle = 'on'; % Ensure the icon is displayed
        end
    
        ax=gca;
        ax.XAxis.FontSize  = 0.8*fontSize;
        ax.YAxis.FontSize  = 0.8*fontSize;
        ax.XLabel.FontSize = 0.8*fontSize;
        ax.YLabel.FontSize = 0.8*fontSize;
        ax.Title.FontSize  = 0.9*fontSize;
        ax.YAxis.TickLength = [0,0];
        %set(gca,'xtick',[])
        set(gca,'ytick',[])
        
        if (ntile-1>=5)
            xlabel('$c_h^\prime$','Interpreter','latex')
        else
            %set(gca,'xtick',[])
            xticklabels([])
        end
        if ntile-1==1 || ntile-1==5
            disp(ntile)
            ylabel('Frequency','Interpreter','latex')
        end
        title(['$p=',num2str(p),'$'],'Interpreter','latex')

        %juste mettre l'axe des abscisses
        %set(h7, 'box', 'off')
        set(gca, 'box', 'off', 'ytick', [], 'xcolor', [0.3,0.3,0.3], 'ycolor',[0.3,0.3,0.3], 'color', 'none');
          
        hold on
    end
end
lgd=legend(' $~b=1$',' $~b=3$',' $~b=10$','Location','NorthWest',...
    'FontSize',0.8*fontSize,'Interpreter','latex','Box','off');
lgd.Layout.Tile = 4;
x1=25 ; x2=25;
lgd.ItemTokenSize = [x1,x2];
%lgd.Layout.Tile = 'east';

%a=(strcat(list_k,','));
%sgtitle(['$c_h^\prime$ in $\{',strcat(a{:}),'\}$'],'Interpreter','latex','FontSize',20)
%sgtitle(['$c_{\texttt{k},h}^\prime$ for \texttt{k} in $\{h,hc,hg,hcg\}$'],'Interpreter','latex','FontSize',ax.Title.FontSize)
%sgtitle(['$c_{\texttt{k},h}^\prime$ for \texttt{k} in $\{hs,hsc,hsg,hscg\}$'],'Interpreter','latex','FontSize',ax.Title.FontSize)
warning('pas de sauvegarde: penser a sauvegarder manuellement')
if 0
    %path_fig = ['C:/Users/Moi/Documents/IPLESP/These/Rapports/PhDThesis/myPhDThesis/imgs/results/grps_hshschsghscg_b_p.pdf'];
    path_fig = ['C:/Users/Moi/Documents/IPLESP/These/Rapports/PhDThesis/myPhDThesis/imgs/results/grps_hhchghcg_b_p_light.pdf'];
    exportgraphics(fig,path_fig,'Resolution',300)
end

%% Comparaison des couts dans un groupe de kits par b par p

clear all; close all;
inf = 'HIV';
vecP = 0:0.1:0.5;
vecB = [1,3,10];
fontSize = 18;
filter = @(data) data(data>-0.8 & data<1000);
%list_k = {'h','hc','hg','hcg'};
list_k = {'hs','hsc','hsg','hscg'};

colors = [[0.6350 0.0780 0.1840];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];...
    [0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.4940 0.1840 0.5560]];

path1   = '..\ParameterAnalysis\analysis_48\rds_1-320\';
costs1  = readtable([path1,'elimCosts_concatenated_rds_1-320_correction.txt']);
path3   = '..\ParameterAnalysis\analysis_59\rds_1-320\';
costs3  = readtable([path3,'elimCosts_concatenated_rds_1-320_correction.txt']);
path10  = '..\ParameterAnalysis\analysis_610\rds_1-320\';
costs10 = readtable([path10,'elimCosts_concatenated_rds_1-320_correction.txt']);

cost_tot = [costs1;costs3;costs10];

%pathRecap  = 'C:\Users\pepiot\Documents\RecapPhD\graphes\comparaison_b\';
pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\comparaison_b\';
pathComp  = ['..\ParameterAnalysis\comparaison_b\'];

%1. histograms par kit
n = length(list_k);
fig = figure('Renderer', 'painters', 'Position', [10 10 800 500]);
tiledlayout(2,2);
for b=vecB
    ech2  = cost_tot(cost_tot.f==b & cost_tot.(inf)<1000,:);
    ax = nexttile;
    echtot=[];
    for j=1:length(vecP)
        ech=[];
        for i=1:n
            k = list_k{i};
            ech3 = ech2(ismember(ech2.kit,k),:);
            %binW = (max(filter(ech3.(inf)))-min(filter(ech3.(inf))))/50;
            binW = (quantile(filter(ech3.(inf)),0.999) - quantile(filter(ech3.(inf)),0.001))/10;
            %binW=0.001;
            %hA = axes(fig);
            
%             binW    =  0.0015*(b==1)+0.002*(b==3)+0.01*(b==10); 
%             low_bnd = -0.035*(b==1)+-0.07*(b==3)+-0.28*(b==10); 
%             up_bnd  = -0.005*(b==1)+-0.02*(b==3)+-0.06*(b==10); 

            p=vecP(j);
            ech = [ech;ech3];
        end
        ech = ech(abs(ech.p-p)<1e-15 ,:).(inf);
        ech_bis = filter(ech(~isnan(ech)));
        if ~isempty(binW) && ~isnan(binW)
            if binW>0
                histogram(ech_bis,'BinWidth',binW,'EdgeColor','none',...
                    'FaceAlpha',0.6,'FaceColor',colors(j,:),'Normalization','probability')
                hold on
            end
        end
        
        echtot = [echtot;ech_bis];
    end

    %box(ax,'off')

    low_bnd = quantile(filter(echtot),0.001);
    up_bnd = quantile(filter(echtot),0.999);
    xlim([low_bnd,up_bnd])
    xticks(unique(round([0.96*low_bnd+0.04*up_bnd, (low_bnd+up_bnd)/2, 0.03*low_bnd+0.97*up_bnd],2)))
    title(['$b=',num2str(b),'$'],'Interpreter','latex','FontSize',fontSize)

    xlabel('$c_h^\prime$','Interpreter','latex')
    ylabel('Frequency','Interpreter','latex')
    title(['$b=',num2str(b),'$'],'Interpreter','latex')

    ax=gca;
    ax.XAxis.FontSize  = 0.7*fontSize;
    ax.YAxis.FontSize  = 0.7*fontSize;
    ax.XLabel.FontSize = 0.8*fontSize;
    ax.YLabel.FontSize = 0.8*fontSize;
    ax.Title.FontSize  = 0.9*fontSize;
    ax.YAxis.TickLength = [0,0];
    %set(gca,'xtick',[])
    set(gca,'ytick',[])

    %juste mettre l'axe des abscisses
    %set(h7, 'box', 'off')
    set(gca, 'box', 'off', 'ytick', [], 'xcolor', [0.3,0.3,0.3], 'ycolor',[0.3,0.3,0.3], 'color', 'none');
        
    if b==1
        lgd=legend(' $~p=0$',' $~p=0.1$',' $~p=0.2$',' $~p=0.3$',' $~p=0.4$',' $~p=0.5$','Location','NorthWest',...
            'FontSize',0.7*fontSize,'Interpreter','latex','Box','off');
        lgd.Layout.Tile = 4;
        x1=25 ; x2=25;
        lgd.ItemTokenSize = [x1,x2];
    end

    
    %lgd.Layout.Tile = 'east';
    %a=(strcat(list_k,','));
    %sgtitle(['$c_h^\prime$ in $\{',strcat(a{:}),'\}$'],'Interpreter','latex','FontSize',20)
    %sgtitle(['$c_{\texttt{k},h}^\prime$ for \texttt{k} in $\{h,hc,hg,hcg\}$'],'Interpreter','latex','FontSize',ax.Title.FontSize)
    sgtitle(['$c_{\texttt{k},h}^\prime$ for \texttt{k} in $\{hs,hsc,hsg,hscg\}$'],'Interpreter','latex','FontSize',ax.Title.FontSize)
end

if(0)
    path_fig = ['C:/Users/Moi/Documents/IPLESP/These/Rapports/PhDThesis/myPhDThesis/imgs/results/grps_hshschsghscg_p_b.pdf'];
    %path_fig = ['C:/Users/Moi/Documents/IPLESP/These/Rapports/PhDThesis/myPhDThesis/imgs/results/grps_hhchghcg_p_b.pdf'];
    exportgraphics(fig,path_fig,'Resolution',400)
end
%%



%% Comparaison des couts groupe hs et groupe sans hs
clear all; close all;
inf ='HIV';
vecP = 0:0.1:0.5;
filter = @(data) data(data>-0.55 & data<1000);
list_k1 = {'h','hc','hg','hcg'};
list_k2 = {'hs','hsc','hsg','hscg'};

fontSize = 24;
set(groot,'defaultAxesTickLabelInterpreter','latex');  

path1   = '..\ParameterAnalysis\analysis_48\rds_1-320\';
costs1  = readtable([path1,'elimCosts_concatenated_rds_1-320_correction.txt']);

%pathRecap  = 'C:\Users\pepiot\Documents\RecapPhD\graphes\comparaison_b\';
pathRecap   = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Recap\graphes\comparaison_b\';
pathComp    = ['..\ParameterAnalysis\comparaison_b\'];

n1 = length(list_k1);
n2 = length(list_k2);

b = 1;

ech1 = costs1(ismember(costs1.kit,list_k1) & costs1.f==b,:);
ech2 = costs1(ismember(costs1.kit,list_k2) & costs1.f==b,:);

fig = figure('Renderer', 'painters', 'Position', [10 10 1250 500]);
tiledlayout(2,3);

ntile = 1;

for p=vecP
    ech1_1 = filter(ech1(abs(ech1.p-p)<1e-14,:).HIV);
    ech2_1 = filter(ech2(abs(ech1.p-p)<1e-14,:).HIV);

    binW = (max([ech1_1;ech2_1])-min([ech1_1;ech2_1]))/20;
    
%     if ntile==4
%         nexttile(ntile+1,[1 1])
%         ntile=ntile+2;
%     else
        nexttile(ntile,[1 1])
        ntile=ntile+1;
%     end
    if ~isnan(binW) & ~isempty(binW)
        histogram(ech1_1,'BinWidth',binW,'EdgeColor','none','Normalization','probability','FaceAlpha',1, 'FaceColor',[0.0, 0.58, 0.71])
        hold on
        histogram(ech2_1,'BinWidth',binW,'EdgeColor','none','Normalization','probability','FaceAlpha',1, 'FaceColor',[0.83,0.65,0.17])
        
        
%         if binW>0
%             
%             echtot = [ech1_bis];
% 
%             low_bnd = -0.25; % quantile(echtot,0.004);
%             up_bnd  = -0.0; % quantile(echtot,0.999);%max(echtot);
%             binW = (up_bnd-low_bnd)./30;
%             
%             histogram(ech1_bis,'BinWidth',binW,'EdgeColor','none','Normalization','probability','FaceAlpha',0.5)
% %             hold on 
% %             histogram(ech3_bis,'BinWidth',binW,'EdgeColor','none','Normalization','probability','FaceAlpha',0.5)
% %             histogram(ech10_bis,'BinWidth',binW,'EdgeColor','none','Normalization','probability','FaceAlpha',0.5)
% 
%             xlim([low_bnd,up_bnd])
%             xticks(round([(0.93*low_bnd+0.03*up_bnd),(low_bnd+up_bnd)/2,(0.03*low_bnd+0.97*up_bnd)],2))
%             %h.Annotation.LegendInformation.IconDisplayStyle = 'on'; % Ensure the icon is displayed
%         end
    
        ax=gca;
        ax.XAxis.FontSize  = 0.65*fontSize;
        ax.YAxis.FontSize  = 0.8*fontSize;
        ax.XLabel.FontSize = 0.8*fontSize;
        ax.YLabel.FontSize = 0.8*fontSize;
        ax.Title.FontSize  = 0.9*fontSize;
        ax.YAxis.TickLength = [0,0];
        %set(gca,'xtick',[])
        set(gca,'ytick',[])
        
%         if (ntile-1>=5)
             xlabel('cost of elimination, $c_h^\prime$','Interpreter','latex')
%         else
%             %set(gca,'xtick',[])
%             xticklabels([])
%         end
%         if ntile-1==1 || ntile-1==5
%             disp(ntile)
            ylabel('Frequency','Interpreter','latex')
%         end
        title(['$p=',num2str(p),'$'],'Interpreter','latex')

        %juste mettre l'axe des abscisses
        %set(h7, 'box', 'off')
        set(gca, 'box', 'off', 'ytick', [], 'xcolor', [0.3,0.3,0.3], 'ycolor',[0.3,0.3,0.3], 'color', 'none');
          
        hold on
    end
end
% lgd=legend(' $~b=1$',' $~b=3$',' $~b=10$','Location','NorthWest','FontSize',0.8*fontSize,'Interpreter','latex','Box','off');

x1=25 ; x2=25;
lgd.ItemTokenSize = [x1,x2];
%lgd.Layout.Tile = 'east';

%a=(strcat(list_k,','));
%sgtitle(['$c_h^\prime$ in $\{',strcat(a{:}),'\}$'],'Interpreter','latex','FontSize',20)
%sgtitle(['$c_{\texttt{k},h}^\prime$ for \texttt{k} in $\{h,hc,hg,hcg\}$'],'Interpreter','latex','FontSize',ax.Title.FontSize)
%sgtitle(['$c_{\texttt{k},h}^\prime$ for \texttt{k} in $\{hs,hsc,hsg,hscg\}$'],'Interpreter','latex','FontSize',ax.Title.FontSize)
warning('pas de sauvegarde: penser a sauvegarder manuellement')
if 0
    %path_fig = ['C:/Users/Moi/Documents/IPLESP/These/Rapports/PhDThesis/myPhDThesis/imgs/results/grps_hshschsghscg_b_p.pdf'];
    path_fig = ['C:/Users/Moi/Documents/IPLESP/These/Rapports/PhDThesis/myPhDThesis/imgs/results/grps_hhchghcg_b_p_light.pdf'];
    exportgraphics(fig,path_fig,'Resolution',300)
end




















































%%
%TO DO below

if (0)
    
    %% Analyse des resultats sur le cout d'elimination
    clear all; close all;
    %----------------%
    f=1; pHIV=0.5;
    folderRun='Run5';
    %----------------%
    path = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\_results'];
    
    sumup_tabAll = [];
    [status, msg, msgID] = mkdir([path,'\_histograms\']);
    path2 = [path,'\_histograms'];
    for strat=1:15
        tot = readtable([path,'\strat_',num2str(strat),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_elim_cost','_concatenated.txt']);
        n=size(tot,1);
        if (1)
            hCt = figure();
            histogram(tot.Ct,12); %Ct
            saveas(hCt,[path2,'\strat_',num2str(strat),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_hist_Ct_n_',num2str(n),'.png'])
            hNg = figure();
            histogram(tot.Ng,12) %Ng
            saveas(hNg,[path2,'\strat_',num2str(strat),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_hist_Ng_n_',num2str(n),'.png'])
            hHIV = figure();
            histogram(tot.HIV,12) %HIV
            saveas(hHIV,[path2,'\strat_',num2str(strat),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_hist_HIV_n_',num2str(n),'.png'])
            hS = figure();
            histogram(tot.syph,12) %syph
            saveas(hS,[path2,'\strat_',num2str(strat),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_hist_S_n_',num2str(n),'.png'])
        end
        close all;
        
        Infection =       ["Ct";               "Ng";               "HIV";              "Syph"];
        Moyenne   = round([mean(tot.Ct);       mean(tot.Ng);       mean(tot.HIV);       mean(tot.syph)],3);
        Mediane   = round([median(tot.Ct);     median(tot.Ng);     median(tot.HIV);     median(tot.syph)],3);
        Std       = round([std(tot.Ct);        std(tot.Ng);        std(tot.HIV);        std(tot.syph)],3); %standard deviation (t-student)
        N         = round([length(tot.Ct);     length(tot.Ng);     length(tot.HIV);     length(tot.syph)],3);
        Min       = round([min(tot.Ct);        min(tot.Ng);        min(tot.HIV);        min(tot.syph)],3);
        Max       = round([max(tot.Ct);        max(tot.Ng);        max(tot.HIV);        max(tot.syph)],3);
        low_95    = round([prctile(tot.Ct,05); prctile(tot.Ng,05); prctile(tot.HIV,05); prctile(tot.syph,05)],3);
        high_95   = round([prctile(tot.Ct,95); prctile(tot.Ng,95); prctile(tot.HIV,95); prctile(tot.syph,95)],3);
        low_90    = round([prctile(tot.Ct,10); prctile(tot.Ng,10); prctile(tot.HIV,10); prctile(tot.syph,10)],3);
        high_90   = round([prctile(tot.Ct,90); prctile(tot.Ng,90); prctile(tot.HIV,90); prctile(tot.syph,90)],3);
        numStrat  = [strat;strat;strat;strat];
        sumup_tab = table(numStrat,Infection,Moyenne,Mediane,Std,N,Min,Max,low_95,high_95,low_90,high_90);
        writetable(sumup_tab,[path,'\strat_',num2str(strat),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_summary.txt'])
        
        sumup_tabAll = [sumup_tabAll;sumup_tab];
    end
    writetable(sumup_tabAll,[path,'\_allStrat_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_summary.txt'])
    
    %% Reshape for LaTex
    clear all;
    %----------------%
    f=1; pHIV=0.5;
    folderRun='Run5';
    %----------------%
    path = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\_results'];
    tab  = readtable([path,'\_allStrat_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_summary.txt']);
    
    pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
    fileID = fopen([pathW,'recap-b-',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt'],'w');
    begin = {'\begin{longtable}{|l|l|l|l|l|l|l|}','\hline',...
        '\textbf{Strat.} & \textbf{Infection} & \textbf{N} & \textbf{Mean} & \textbf{Std} & \textbf{low\_95} & \textbf{high\_95} \\ \hline \hline'};
    
    for k=1:length(begin)
        fprintf(fileID,'%12s\r\n',[begin{k}]);
    end
    
    for k=1:15
        for i=1:4
            towrite = tab(4*(k-1)+i,{'numStrat','Infection','N','Moyenne','Std','low_95','high_95'});
            fprintf(fileID,'%7s',[num2str(k),' & ']);
            fprintf(fileID,'%7s',[towrite.Infection{:},' & ']);
            towrite2=strjoin(strsplit(num2str(table2array(towrite(:,3:end)))),' & ');
            fprintf(fileID,'%12s',[towrite2, '\\ \hline']);
        end
        if k<15
            fprintf(fileID,'%12s\r\n','\hline');
        end
    end
    fprintf(fileID,'%12s\r\n','\end{longtable}');
    
    fclose(fileID);
    
    %% Comparaison des moyennes des couts (Student test and others)
    % Distribution de reference : strategie [1], [2], [3], [4]
    clear all;
    %-----------------%
    f=1; pHIV=0.5;
    folderRun='Run5';
    %-----------------%
    pathInput = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\_results\'];
    
    createStatTests_v2(f,pHIV,pathInput,"Student")
    createStatTests_v2(f,pHIV,pathInput,"MWW")
    createStatTests_v2(f,pHIV,pathInput,"KS")
    createStatTests_v2(f,pHIV,pathInput,"StudentMatlab")
    
    %% Table for LaTex - matrix of p-values
    clear all;
    %----------------%
    f=1; pHIV=0.5;
    folderRun='Run5';
    pathInput = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\_results\_tests\'];
    %----------------%
    for infection = {'Ct','Ng','HIV','S'}
        for test={'KS','MWW','StudentMatlab','Student'}
            tab = readtable([pathInput,'_',test{:},'_p_values_',infection{:},'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt']);
            pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
            fileID = fopen([pathW,'pvalues-b-',num2str(f),'_p_',num2str(round(pHIV*100)),'-',test{:},'-',infection{:},'.txt'],'w');
            
            %Initialization
            %formatSpec = 'X is %4.2f meters or %8.3f mm\n';
            begin = {'\begin{sidewaystable}','\centering', ['\caption{',num2str(test{:}),' - b=',num2str(f),' - ',infection{:},'}'],...
                '\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|l|l|l|l|l|}','\hline'};
            for k=1:length(begin)
                fprintf(fileID,'%12s\r\n',[begin{k}]);
                %fprintf(fileID,'%6.2f %12.8f\n',A);
            end
            fprintf(fileID,'%12s\r\n','Strat & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 & 12 & 13 & 14 & 15\\ \hline');
            
            for k=1:15
                fprintf(fileID,'%4s',[num2str(k),' & ']);
                pvalues = tab(k,k:15);
                if k~=1
                    for i=1:(k-1)
                        fprintf(fileID,'%2s','~ &');
                    end
                end
                towrite = strjoin(strsplit(num2str(table2array(pvalues))),' & ');
                fprintf(fileID,'%5s \r\n',[towrite, '\\ \hline']);
            end
            
            fprintf(fileID,'%12s\r\n','\end{tabular}');
            fprintf(fileID,'%12s\r\n','\end{sidewaystable}');
            
            fclose(fileID);
        end
    end
    
    %% Histograms of costs strat A vs strat B
    clear all;
    %---------%
    pHIV=0.5;
    vecF=[1];
    folderRun = 'Run5';
    %---------%
    pathInput = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\_results\'];
    
    for f=vecF
        close all
        %maxNBins=15;
        minBinWidth=0.005;
        figure(1);
        
        %to define x-axis histograms :
        minCt=[];maxCt=[];
        minNg=[];maxNg=[];
        minHIV=[];maxHIV=[];
        minS=[];maxS=[];
        nStrat=15;
        
        fex=figure();
        Lf = 1000;%fex.Position(3); %largeur figure
        Ls = Lf/(nStrat+7); %largeur subplot
        a  = (Lf - nStrat*Ls)/(nStrat+1+1); %espace horizontal entre les subplots;
        
        hf = 700;%fex.Position(4);  %hauteur figure
        hs = hf/(nStrat+8); %hauteur subplot (10)
        b  = (hf - nStrat*hs)/(nStrat+3); %(10) %espace vertical entre les subplots;
        close all
        
        for stratA=1:nStrat
            path2 = [pathInput,'_histograms\'];
            stratAfile = readtable([pathInput,'strat_',num2str(stratA),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_elim_cost_concatenated.txt']);
            for stratB=stratA:nStrat
                plotPos = [round(a+(stratB-1)*(a+Ls),3),round(-2*b + hf-b-hs-(stratA-1)*(hs+b),3), Ls hs];
                if stratA==stratB
                    %if(0)
                    fig1 = figure(1);
                    fig1.Position = [10 10 Lf hf];
                    %h1 = subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    h1 = axes('Position', plotPos);
                    %nBins = min(length([min(stratAfile.Ct):0.002:max(stratAfile.Ct)]),maxNBins);
                    BinWidth=min([minBinWidth,(max(stratAfile.Ct)-min(stratAfile.Ct))./10]);
                    histogram(stratAfile.Ct,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
                    %minCt = min(minCt,min(stratAfile.Ct));
                    %maxCt = max(maxCt,max(stratAfile.Ct));
                    %set(gca,'xtick',[])
                    set(gca,'ytick',[])
                    set(gca,'fontsize',4)
                    set(h1,'Unit','pixels','Position',plotPos);
                    %hold on
                    %plot([0,0],[0;100],'k')
                    
                    fig2 = figure(2);
                    fig2.Position = [10 10 Lf hf];
                    %h2=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    h2 = axes('Position', plotPos);
                    %nBins = min(length([min(stratAfile.Ng):0.002:max(stratAfile.Ng)]),maxNBins);
                    BinWidth=min([minBinWidth,(max(stratAfile.Ng)-min(stratAfile.Ng))./10]);
                    histogram(stratAfile.Ng,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
                    minNg = min(minNg,min(stratAfile.Ng));
                    maxNg = max(maxNg,max(stratAfile.Ng));
                    %set(gca,'xtick',[])
                    set(h2,'Unit','pixels','Position',plotPos);
                    set(gca,'ytick',[])
                    set(gca,'fontsize',4)
                    %end
                    
                    fig3 = figure(3);
                    fig3.Position = [10 10 Lf hf];
                    %h3 = subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    h3 = axes('Position', plotPos);
                    %nBins = min(length([min(stratAfile.HIV):0.002:max(stratAfile.HIV)]),maxNBins);
                    BinWidth=min([minBinWidth,(max(stratAfile.HIV)-min(stratAfile.HIV))./10]);
                    histogram(stratAfile.HIV,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
                    minHIV = min(minHIV,min(stratAfile.HIV));
                    maxHIV = max(maxHIV,max(stratAfile.HIV));
                    %set(gca,'xtick',[])
                    set(h3,'Unit','pixels','Position',plotPos);
                    set(gca,'ytick',[])
                    set(gca,'fontsize',4)
                    
                    
                    fig4 = figure(4);
                    fig4.Position = [10 10 Lf hf];
                    %h4=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    h4=axes('Position', plotPos);
                    %nBins = min(length([min(stratAfile.syph):0.002:max(stratAfile.syph)]),maxNBins);
                    BinWidth=min([minBinWidth,(max(stratAfile.syph)-min(stratAfile.syph))./10]);
                    histogram(stratAfile.syph,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
                    minS = min(minS,min(stratAfile.syph));
                    maxS = max(maxS,max(stratAfile.syph));
                    %set(gca,'xtick',[])
                    set(h4,'Unit','pixels','Position',plotPos);
                    set(gca,'ytick',[])
                    set(gca,'fontsize',4)
                    
                else
                    stratBfile = readtable([pathInput,'strat_',num2str(stratB),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_elim_cost_concatenated.txt']);
                    
                    %Ct
                    figure(1);
                    %h5=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    h5 = axes('Position', plotPos);
                    %nBins = min(length([min(stratAfile.Ct):0.002:max(stratAfile.Ct)]),maxNBins);
                    %BinWidth=max([minBinWidth,min([(max(stratBfile.Ct)-min(stratBfile.Ct))./10,(max(stratAfile.Ct)-min(stratAfile.Ct))./10])]);
                    BinWidth= (max([stratBfile.Ct;stratAfile.Ct])-min([stratAfile.Ct;stratBfile.Ct]))/12;
                    histogram(h5,stratAfile.Ct,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
                    hold on
                    %nBins = min(length([min(stratBfile.Ct):0.002:max(stratBfile.Ct)]),maxNBins,'Normalization','probability');
                    histogram(h5,stratBfile.Ct,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth)
                    %axis off
                    set(gca,'xtick',[])
                    set(gca,'ytick',[])
                    set(h5,'Unit','pixels','Position',plotPos);
                    fprintf(['x=',num2str(a+(stratB-1)*(a+Ls)), ' y=',num2str(hf-b-hs-(stratA-1)*(hs+b)),'\n'])
                    %h=gca; h.XAxis.TickLength = [0 0];
                    %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                    %saveas(hCt,[path2,'hist_Ct_Strat_',num2str(stratA),'(',num2str(length(stratAfile.Ct)),')_vs_',...
                    %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.Ct)),')',...
                    %    '_b_',num2str(f),'.png'])
                    
                    figure(2);
                    %h6=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    h6 = axes('Position',plotPos);
                    %nBins = min(length([min(stratAfile.Ng):0.002:max(stratAfile.Ng)]),maxNBins);
                    %BinWidth=min([minBinWidth,(max(stratBfile.Ng)-min(stratBfile.Ng))./10,(max(stratAfile.Ng)-min(stratAfile.Ng))./10]);
                    BinWidth= (max([stratBfile.Ng;stratAfile.Ng])-min([stratAfile.Ng;stratBfile.Ng]))/12;
                    histogram(stratAfile.Ng,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
                    hold on
                    %nBins = min(length([min(stratBfile.Ng):0.002:max(stratBfile.Ng)]),maxNBins);
                    histogram(stratBfile.Ng,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth);
                    %axis off
                    set(h6,'Unit','pixels','Position',plotPos);
                    set(gca,'xtick',[])
                    set(gca,'ytick',[])
                    %set(gca,'xtick',[])
                    %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                    %saveas(hNg,[path2,'hist_Ng_Strat_',num2str(stratA),'(',num2str(length(stratAfile.Ng)),')_vs_',...
                    %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.Ng)),')',...
                    %    '_b_',num2str(f),'.png'])
                    
                    
                    figure(3);
                    %h7=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    h7=axes('Position', plotPos);
                    %BinWidth=min([minBinWidth,(max(stratBfile.HIV)-min(stratBfile.HIV))./10,(max(stratAfile.HIV)-min(stratAfile.HIV))./10]);
                    BinWidth= (max([stratBfile.HIV;stratAfile.HIV])-min([stratAfile.HIV;stratBfile.HIV]))/12;
                    %nBins = min(length([min(stratAfile.HIV):0.002:max(stratAfile.HIV)]),maxNBins);
                    histogram(h7,stratAfile.HIV,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
                    hold on
                    %nBins = min(length([min(stratBfile.HIV):0.002:max(stratBfile.HIV)]),maxNBins);
                    histogram(h7,stratBfile.HIV,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth);
                    %axis off
                    set(gca,'XTick',[])
                    set(gca,'ytick',[])
                    set(h7,'Unit','pixels','Position',plotPos);
                    %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                    %saveas(hHIV,[path2,'hist_HIV_Strat_',num2str(stratA),'(',num2str(length(stratAfile.HIV)),')_vs_',...
                    %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.HIV)),')',...
                    %    '_b_',num2str(f),'.png'])
                    
                    
                    figure(4);
                    %h8=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                    h8 = axes('Position', plotPos);
                    %BinWidth=min([minBinWidth,(max(stratBfile.syph)-min(stratBfile.syph))./10,(max(stratAfile.syph)-min(stratAfile.syph))./10]);
                    BinWidth=(max([stratBfile.syph;stratAfile.syph])-min([stratAfile.syph;stratBfile.syph]))/12;
                    %nBins = min(length([min(stratAfile.syph):0.002:max(stratAfile.syph)]),maxNBins);
                    histogram(stratAfile.syph,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
                    hold on
                    %nBins = min(length([min(stratBfile.syph):0.002:max(stratBfile.syph)]),maxNBins);
                    histogram(stratBfile.syph,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth);
                    %axis off
                    set(gca,'xtick',[])
                    set(gca,'ytick',[])
                    set(h8,'Unit','pixels','Position',plotPos);
                    %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                    %saveas(hS,[path2,'hist_S_Strat_',num2str(stratA),'(',num2str(length(stratAfile.syph)),')_vs_',...
                    %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.syph)),')',...
                    %    '_b_',num2str(f),'.png'])
                    
                    %close all;
                    %pause(1)
                    
                end
            end
        end
        
        %add number of the strategy
        for k=1:15
            for numFig=1:4
                figure(numFig)
                Ct = (k*a+(k-1)*Ls )/Lf;
                y1 = 1-(hf-b-hs-14*(hs+b))/hf;
                annotation('textbox', [Ct y1 Ls/Lf hs/hf], 'string', num2str(k),...
                    'VerticalAlignment','middle','HorizontalAlignment','center',...
                    'LineStyle','none')
                Ng = 1-2*a/Lf;
                y2 = 1-(4*b+b/3+k*b+(k-1)*hs)/hf;
                annotation('textbox', [Ng y2 Ls/Lf hs/hf], 'string', num2str(k),...
                    'VerticalAlignment','middle','HorizontalAlignment','center',...
                    'LineStyle','none')
            end
        end
        
        pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\graphes\';
        
        %Ct
        saveas(fig1,[path2,'hist_Ct_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.png'])
        set(fig1,'PaperOrientation','landscape');
        saveas(fig1,[pathRecap,'hist_Ct_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.pdf'])
        
        %Ng
        saveas(fig2,[path2,'hist_Ng_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.png'])
        set(fig2,'PaperOrientation','landscape');
        saveas(fig2,[pathRecap,'hist_Ng_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.pdf'])
        
        %HIV
        saveas(fig3,[path2,'hist_HIV_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.png'])
        set(fig3,'PaperOrientation','landscape');
        saveas(fig3,[pathRecap,'hist_HIV_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.pdf'])
        
        %Syphilis
        saveas(fig4,[path2,'hist_S_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.png'])
        set(fig4,'PaperOrientation','landscape');
        saveas(fig4,[pathRecap,'hist_S_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.pdf'])
    end
    
    
    %%
    clear all; close all;
    tabCorresp = zeros(15,15,4);
    %Ct
    tabCorresp(1,5,1)=1; tabCorresp(1,6,1)=1; tabCorresp(1,7,1)=1; tabCorresp(1,14,1)=1;
    tabCorresp(2,8,1)=1; tabCorresp(3,9,1)=1; tabCorresp(4,10,1)=1;tabCorresp(5,6,1)=1;
    tabCorresp(5,7,1)=1; tabCorresp(5,14,1)=1;tabCorresp(6,14,1)=1;tabCorresp(7,14,1)=1;
    
    %Ng
    tabCorresp(1,3,2)=1; tabCorresp(1,4,2)=1; tabCorresp(1,7,2)=1; tabCorresp(1,13,2)=1;
    tabCorresp(2,8,2)=1; tabCorresp(3,4,2)=1; tabCorresp(3,7,2)=1; tabCorresp(3,13,2)=1;
    tabCorresp(4,7,2)=1; tabCorresp(4,13,2)=1;tabCorresp(5,10,2)=1;tabCorresp(6,9,2)=1;
    tabCorresp(7,13,2)=1;
    
    %HIV
    tabCorresp(1,2,3)=1; tabCorresp(1,4,3)=1; tabCorresp(1,6,3)=1; tabCorresp(1,12,3)=1;
    tabCorresp(2,4,3)=1; tabCorresp(2,6,3)=1; tabCorresp(2,12,3)=1;tabCorresp(3,9,3)=1;
    tabCorresp(4,6,3)=1; tabCorresp(4,12,3)=1;tabCorresp(5,10,3)=1;tabCorresp(6,12,3)=1;
    tabCorresp(7,8,3)=1;
    
    %syphilis
    tabCorresp(1,2,4)=1; tabCorresp(1,3,4)=1; tabCorresp(1,5,4)=1; tabCorresp(1,11,4)=1;
    tabCorresp(2,3,4)=1; tabCorresp(2,5,4)=1; tabCorresp(2,11,4)=1;tabCorresp(3,5,4)=1;
    tabCorresp(3,11,4)=1; tabCorresp(4,10,4)=1;tabCorresp(5,11,4)=1;tabCorresp(6,9,4)=1;
    tabCorresp(7,8,4)=1;
    
    for i=1:4
        for k=1:15
            tabCorresp(k,k,i)=1;
            for m=(k+1):15
                tabCorresp(m,k,i)=tabCorresp(k,m,i); %symetrie par rapport a la diag
            end
        end
    end
    
    %% For each strategy : looking for other strategies that can be candidates to make a group of strategies (each line)
    % (resTest-b-**.txt)
    
    close all; clearvars -except tabCorresp
    %-------------------------%
    p=0.05; %p-value threshold
    f=1; test='KS';
    pHIV=0.5;
    folderRun = 'Run5';
    %-------------------------%
    
    pathInput = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\_results\_tests\'];
    pvCt_file = readtable([pathInput,'_',test,'_p_values_Ct_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt']);
    pvNg_file = readtable([pathInput,'_',test,'_p_values_Ng_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt']);
    pvHIV_file = readtable([pathInput,'_',test,'_p_values_HIV_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt']);
    pvS_file = readtable([pathInput,'_',test,'_p_values_S_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt']);
    
    pvCt  = table2array(pvCt_file(1:15,1:15))>p;
    pvNg  = table2array(pvNg_file(1:15,1:15))>p;
    pvHIV = table2array(pvHIV_file(1:15,1:15))>p;
    pvS   = table2array(pvS_file(1:15,1:15))>p;
    
    %Candidate groups of strategy (each line)
    listStrat = 1:15;
    for k=1:15
        gpCt{k} = listStrat(pvCt(k,1:15));
        gpNg{k} = listStrat(pvNg(k,1:15));
        gpHIV{k} = listStrat(pvHIV(k,1:15));
        gpS{k} = listStrat(pvS(k,1:15));
    end
    
    matP(:,:,1) = pvCt;
    matP(:,:,2) = pvNg;
    matP(:,:,3) = pvHIV;
    matP(:,:,4) = pvS;
    
    grps{1,:} = gpCt;
    grps{2,:} = gpNg;
    grps{3,:} = gpHIV;
    grps{4,:} = gpS;
    
    % Writing in LaTeX file
    test='KS';
    pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
    for infection={'Ct','Ng','HIV','S'}
        fileID = fopen([pathW,'resTest-b-',num2str(f),'_p_',num2str(round(pHIV*100)),'-',test,'-',infection{:},'.txt'],'w');
        numInfs=1:4; numInf = numInfs(strcmp(infection{:},{'Ct','Ng','HIV','S'}));
        pInf = matP(:,:,numInf);
        
        %Initialization
        begin = {'\begin{table}[h]','\centering', ['\caption{',test,' - b=',num2str(f),'p=',num2str(pHIV),' - ',infection{:},' - pvalue threshold=',num2str(p),'}'],...
            '\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|l|l|l|l|l|l|}','\hline'};
        for k=1:length(begin)
            fprintf(fileID,'%12s\r\n',[begin{k}]);
        end
        fprintf(fileID,'%12s\r\n','Strat & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 & 12 & 13 & 14 & 15 & group\\ \hline');
        
        for k=1:15
            fprintf(fileID,'%4s',[num2str(k),' & ']);
            pvalues = pInf(k,k:15);
            if k~=1
                for i=1:(k-1)
                    fprintf(fileID,'%2s','~ &');
                end
            end
            towrite = strjoin(strsplit(num2str(pvalues)),' & ');
            fprintf(fileID,'%5s ',towrite);
            fprintf(fileID,'%12s \r\n',['& [',num2str(grps{numInf}{k}),']', '\\ \hline']);
        end
        
        fprintf(fileID,'%12s\r\n','\end{tabular}');
        fprintf(fileID,'%12s\r\n','\end{table}');
        
        fclose(fileID);
    end
    
    
    %% Identifying problems (conflitctsTable.txt)
    % when the p-value should be big (same model used = 1), but is low instead
    %clear all
    clearvars -except tabCorresp
    %------------------------------------------------
    infections={'Ct','Ng','HIV','S'}; strategies=1:15;
    p = 0.05; %p-value threshold
    test = 'KS';
    vecpHIV = 0.5;
    vecF = 1;
    folderRun = 'Run5';
    %------------------------------------------------
    for f=vecF
        for pHIV=vecpHIV
            for test={'KS','MWW','Student'}
                i=1;
                pathInput = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\_results\_tests\'];
                pvCt_file = readtable([pathInput,'_',test{:},'_p_values_Ct_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt']); pMat(:,:,1)=table2array(pvCt_file(:,1:15));
                pvNg_file = readtable([pathInput,'_',test{:},'_p_values_Ng_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt']); pMat(:,:,2)=table2array(pvNg_file(:,1:15));
                pvHIV_file = readtable([pathInput,'_',test{:},'_p_values_HIV_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt']); pMat(:,:,3)=table2array(pvHIV_file(:,1:15));
                pvS_file = readtable([pathInput,'_',test{:},'_p_values_S_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt']); pMat(:,:,4)=table2array(pvS_file(:,1:15));
                
                for numInf=1:4
                    temp = (tabCorresp(:,:,numInf) - (pMat(:,:,numInf)>p));
                    %if 1, then the difference has been found significative whereas it
                    %should'nt be ((1-0).
                    
                    [x,y]=find(temp==1);
                    k=1;
                    while k<=length(x)/2
                        var_b(i) = f;
                        var_pHIV(i) = pHIV;
                        var_infection{i} = infections{numInf};
                        var_test{i}      = test;
                        var_strat1(i)    = strategies(x(k));
                        var_strat2(i)    = strategies(y(k));
                        var_pvalue(i)    = pMat(x(k),y(k),numInf);
                        i=i+1;
                        k=k+1;
                    end
                end
            end
        end
    end
    tab = table(var_b',var_pHIV',var_infection',var_test',var_strat1',var_strat2',var_pvalue');
    
    % For LaTex / conflict table
    pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
    fileID = fopen([pathW,'conflitctsTable.txt'],'w');
    
    %Initialization
    begin = {'\begin{table}[h]','\centering', ...
        '\begin{tabular}{|l|l|l|l|l|l|}','\hline'};
    for k=1:length(begin)
        fprintf(fileID,'%12s\r\n',begin{k});
    end
    fprintf(fileID,'%12s\r\n','b & p & infection & test & strat & strat & p-value \\ \hline \hline');
    
    for i=1:size(tab,1)
        fprintf(fileID,'%4s',[num2str(var_b(i)),' & ']);
        fprintf(fileID,'%4s',[num2str(var_pHIV(i)),' & ']);
        fprintf(fileID,'%10s',[var_infection{i},' & ']);
        fprintf(fileID,'%10s',[var_test{i}{:},' & ']);
        fprintf(fileID,'%4s',[num2str(var_strat1(i)),' & ']);
        fprintf(fileID,'%4s',[num2str(var_strat2(i)),' & ']);
        fprintf(fileID,'%5s \r\n',[num2str(var_pvalue(i)),' \\ \hline']);
    end
    
    fprintf(fileID,'%12s\r\n','\end{tabular}');
    fprintf(fileID,'%12s\r\n','\end{table}');
    
    fclose(fileID);
    
    %% ----------------------------------------
    
    %% Analysis by group (recap-per-group-b-**.txt)
    clear all; close all;
    %----------%
    % TO DO : adapter les groupes ci-dessous en fonction des resultats
    f=1;
    pHIV=0.5;
    folderRun='Run5';
    %----------%
    
    for inf=1:4
        if f==1%
            groups{1} = {[1 5 6 7 14],[2 8 12],[3 9 13],[4,10],[11,15]}; %Ct
            groups{2} = {[1,3,4,7,13],[2,8,11,12,15],[5,6,9,10,14]}; %Ng
            groups{3} = {[1,2,4,6,12],[3,9,13],[7,8],[5,10,14],[11],[15]}; %HIV
            groups{4} = {[1,2,3,5,11],[4,10],[6,9],[7,8],[12],[13],[14],[15]}; %S
        elseif f==3
            groups{1} = {[1 5 6 7 14],[2 8 12], [3 9 13],[4,10],[11,15]};
            groups{2} = {[1,3,4,7,13],[2,8,11,12,15],[5,6,9,10,14]};
            groups{3} = {[1,2,4,6,12],[3,9,13],[7,8],[5,10,14],[11,15]};
            groups{4} = {[1,2,3,5,11],[4,10],[6,9],[7,8],[12],[13],[14],[15]};
        elseif f==10
            groups{1} = {[1 5 6 7 14],[2 8 12],[3 9 13],[4,10],[11,15]};
            groups{2} = {[1,3,4,7,13],[2,8,11,12,15],[5,6,9,10,14]};
            groups{3} = {[1,2,4,6,7,8,12],[3,9,13],[5,10,14],[11,15]};
            groups{4} = {[1,2,3,5,11],[4,10],[6,9],[7,8],[12],[13],[14],[15]};
        end
        pathInput = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\_results\'];
        i=1; tab=[];
        
        for g=groups{inf}
            tot=[];
            for strat=g{:}
                tot = [tot;readtable([pathInput,'strat_',num2str(strat),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_elim_cost_concatenated.txt'])];
            end
            Infection =       ["Ct";               "Ng";               "HIV";              "Syph"];
            Moyenne   = round([mean(tot.Ct);       mean(tot.Ng);       mean(tot.HIV);       mean(tot.syph)],3);
            Mediane   = round([median(tot.Ct);     median(tot.Ng);     median(tot.HIV);     median(tot.syph)],3);
            Std       = round([std(tot.Ct);        std(tot.Ng);        std(tot.HIV);        std(tot.syph)],3); %standard deviation (t-student)
            N         = round([length(tot.Ct);     length(tot.Ng);     length(tot.HIV);     length(tot.syph)],3);
            Min       = round([min(tot.Ct);        min(tot.Ng);        min(tot.HIV);        min(tot.syph)],3);
            Max       = round([max(tot.Ct);        max(tot.Ng);        max(tot.HIV);        max(tot.syph)],3);
            low_95    = round([prctile(tot.Ct,2.5); prctile(tot.Ng,2.5); prctile(tot.HIV,2.5); prctile(tot.syph,2.5)],3);
            high_95   = round([prctile(tot.Ct,97.5); prctile(tot.Ng,97.5); prctile(tot.HIV,97.5); prctile(tot.syph,97.5)],3);
            low_90    = round([prctile(tot.Ct,05); prctile(tot.Ng,05); prctile(tot.HIV,10); prctile(tot.syph,10)],3);
            high_90   = round([prctile(tot.Ct,95); prctile(tot.Ng,95); prctile(tot.HIV,90); prctile(tot.syph,90)],3);
            group = string(repmat(['[',num2str(g{:}),']'],4,1));
            sumup_tab = table(group,Infection,Moyenne,Mediane,Std,N,Min,Max,low_95,high_95,low_90,high_90);
            tab = [tab;sumup_tab];
            i=i+1;
        end
        
        % For LaTex
        infections = {'Ct','Ng','HIV','S'};
        
        pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
        fileID = fopen([pathW,'recap-per-group-b-',num2str(f),'_p_',num2str(round(pHIV*100)),'-',num2str(infections{inf}),'.txt'],'w');
        begin = {'\begin{longtable}[H]{|l|l|l|l|l|l|l|}',['\caption{groups based on ',infections{inf},' analysis b=',num2str(f),' and p=',num2str(pHIV),' } \\'],'\hline',...
            '\textbf{Group} & \textbf{Infection} & \textbf{N} & \textbf{Mean} & \textbf{Std} & \textbf{low\_95} & \textbf{high\_95} \\ \hline \hline'};
        
        for k=1:length(begin)
            fprintf(fileID,'%12s\r\n',[begin{k}]);
        end
        
        for k=1:length(groups{inf})
            for i=1:4
                towrite = tab(4*(k-1)+i,{'group','Infection','N','Moyenne','Std','low_95','high_95'});
                fprintf(fileID,'%10s',[num2str(groups{inf}{k}),' & ']);
                fprintf(fileID,'%7s',[towrite.Infection{:},' & ']);
                towrite2=strjoin(strsplit(num2str(table2array(towrite(:,3:end)))),' & ');
                fprintf(fileID,'%12s \r\n',[towrite2, '\\ \hline']);
            end
            if k<length(groups{inf})
                fprintf(fileID,'%12s\r\n','\hline');
            end
        end
        fprintf(fileID,'%12s\r\n','\end{longtable}');
        
        fclose(fileID);
    end
    
    %% Analysis by group  (short version) (recap-per-group-b-**-**-short.txt)
    clearvars -except groups f
    close all;
    %----------%
    f=1;
    pHIV=0.5;
    folderRun='Run5';
    % adapter les groupes ci-dessous en fonction des resultats
    %----------%
    
    for inf=1:4
        pathInput = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\_results\'];
        
        i=1; tab=[];
        
        for g=groups{inf}
            tot=[];
            for strat=g{:}
                tot = [tot;readtable([pathInput,'strat_',num2str(strat),'_b_',num2str(f),'_p_',num2str(pHIV*100),'_elim_cost_concatenated.txt'])];
            end
            Infection =       ["Ct";               "Ng";               "HIV";              "Syph"];
            Moyenne   = round([mean(tot.Ct);       mean(tot.Ng);       mean(tot.HIV);       mean(tot.syph)],3);
            Mediane   = round([median(tot.Ct);     median(tot.Ng);     median(tot.HIV);     median(tot.syph)],3);
            N         = round([length(tot.Ct);     length(tot.Ng);     length(tot.HIV);     length(tot.syph)],3);
            low_95    = round([prctile(tot.Ct,2.5); prctile(tot.Ng,2.5); prctile(tot.HIV,2.5); prctile(tot.syph,2.5)],3);
            high_95   = round([prctile(tot.Ct,97.5); prctile(tot.Ng,97.5); prctile(tot.HIV,97.5); prctile(tot.syph,97.5)],3);
            group = string(repmat(['(',num2str(g{:}),')'],4,1));
            sumup_tab = table(group,Infection,Moyenne,Mediane,N,low_95,high_95);
            tab = [tab;sumup_tab];
            i=i+1;
        end
        
        % For LaTex
        infections = {'Ct','Ng','HIV','S'};
        
        pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
        fileID = fopen([pathW,'recap-per-group-b-',num2str(f),'_p_',num2str(pHIV*100),'-',num2str(infections{inf}),'-short.txt'],'w');
        begin = {'\begin{longtable}[H]{|l|l|l|}',['\caption{groups based on ',infections{inf},' analysis b=',num2str(f),' and p=',num2str(pHIV),' } \\'],'\hline',...
            '\textbf{group} & \textbf{n} & \textbf{mean}~\textbf{[95\% CI]} \\ \hline \hline'};
        
        for k=1:length(begin)
            fprintf(fileID,'%12s\r\n',[begin{k}]);
        end
        
        for k=1:length(groups{inf})
            for i=inf%1:4
                towrite = tab(4*(k-1)+i,{'group','Infection','N','Moyenne','low_95','high_95'});
                towrite2=['(',strjoin(strsplit(num2str(groups{inf}{k})),','),')'];
                %fprintf(fileID,'%10s',[num2str(groups{inf}{k}),' & ']);
                fprintf(fileID,'%10s',[towrite2,' & ']);
                fprintf(fileID,'%7s',[num2str(towrite.N),' & ']);
                fprintf(fileID,'%7s',[num2str(towrite.Moyenne),'~[',num2str(towrite.low_95),',',num2str(towrite.high_95),']']);
                fprintf(fileID,'%12s \r\n',['\\ \hline']);
            end
            %         if k<length(groups{inf})
            %             fprintf(fileID,'%12s\r\n','\hline');
            %         end
        end
        fprintf(fileID,'%12s\r\n','\end{longtable}');
        
        fclose(fileID);
    end
    
end %if(0)