%Analyse des résultats

%Creation des fichiers recap pour chaque strategie
% TO DO : ajouter les resultats des fichiers MAC
clear all
f=10; %beep off %last time : 16h le 19/04
pathtot = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v2\_results\_b_',num2str(f),'\'];
for strat=1:15   
    
    tot = [];
    
    %PC
    path  = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\_sauvegarde du 23-03-22\_b_',num2str(f),'\'];
    if f==1
        for round=[2,3,4,7,10]
            path2 = [path,'_round_',num2str(round),'_100','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_100','.txt'];
            A = readtable(path2,'ReadVariableNames', true);
            tot = [tot;A];
        end
        if sum([1,2,3,6,9,11,13]==strat)>0
            path2 = [path,'_round_17_100','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_100','.txt'];
            A = readtable(path2,'ReadVariableNames', true);
            tot = [tot;A];
        end
        if sum([2,6,7,13]==strat)>0
            path2 = [path,'_round_18_300','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_300','.txt'];
            A = readtable(path2,'ReadVariableNames', true);
            tot = [tot;A];
        end
        if sum([3,9,13]==strat)>0
            for num=28
                path2 = [path,'_round_',num2str(num),'_300','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_300','.txt'];
                A = readtable(path2,'ReadVariableNames', true);
                tot = [tot;A];
            end
        end
    end
    if f==3
        for round=[5,8,11,13]
            path2 = [path,'_round_',num2str(round),'_100','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_100','.txt'];
            A = readtable(path2,'ReadVariableNames', true);
            tot = [tot;A];
        end  
        if sum([1,7,8]==strat)>0
            path2 = [path,'_round_14_3000','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_3000','.txt'];
            A = readtable(path2,'ReadVariableNames', true);
            tot = [tot;A];
        end
        if sum([1,3,5,7]==strat)>0
            for num=19:23
                path2 = [path,'_round_',num2str(num),'_300','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_300','.txt'];
                A = readtable(path2,'ReadVariableNames', true);
                tot = [tot;A];
            end
        end
        if sum([1,5,14]==strat)>0
            for num=24:25
                path2 = [path,'_round_',num2str(num),'_300','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_300','.txt'];
                A = readtable(path2,'ReadVariableNames', true);
                tot = [tot;A];
            end
        end
        if sum([3,9,13]==strat)>0
            for num=26:27
                path2 = [path,'_round_',num2str(num),'_300','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_300','.txt'];
                A = readtable(path2,'ReadVariableNames', true);
                tot = [tot;A];
            end
        end
    end

    if f==10  
        for round=[6,9,12,14]
            path2 = [path,'_round_',num2str(round),'_100','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_100','.txt'];
            A = readtable(path2,'ReadVariableNames', true);
            tot = [tot;A];
        end
        if sum([1,7,8]==strat)>0
            path2 = [path,'_round_16_3000','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_3000','.txt'];
            A = readtable(path2,'ReadVariableNames', true);
            tot = [tot;A];
        end
        if strat==7
            path2 = [path,'_round_17_3000','\Strat_',num2str(strat),'_b_',num2str(f),'_tot_3000','.txt'];
            A = readtable(path2,'ReadVariableNames', true);
            tot = [tot;A];
        end
    end
    
    %MAC 
    if f==1 || f==3
        path3 = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesMAC\_resultats\_sauvegarde temporaire\_b_',num2str(f),'\'];
        nbSim = [1000,500];
        for k=1:2
            %path4 = [path3,'_round_',num2str(k),'_b',num2str(f),'_',num2str(nbSim(k)),'\','Strat_',num2str(strat),'_b_',num2str(f),'_from_1_to_',num2str(nbSim(k)),'.txt'];
            path4 = [path3,'_round_',num2str(k),'_b_',num2str(f),'\','Strat_',num2str(strat),'_b_',num2str(f),'_from_1_to_',num2str(nbSim(k)),'.txt'];
            B = readtable(path4,'ReadVariableNames', true);
            tot = [tot;B];
        end
    end
    if f==10
        if strat<=14
            path3 = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesMAC\_resultats\_sauvegarde temporaire\_b_',num2str(f),'\'];
            nbSim = [1000,500];
            for k=1:2
                path4 = [path3,'_round_',num2str(k),'_b_',num2str(f),'\','Strat_',num2str(strat),'_b_',num2str(f),'_from_1_to_',num2str(nbSim(k)),'.txt'];
                B = readtable(path4,'ReadVariableNames', true);
                tot = [tot;B];
            end
        elseif strat==15
            path4 = [path3,'_round_1_b_',num2str(f),'\','Strat_15_b_',num2str(f),'_from_1_to_500.txt'];
            B = readtable(path4,'ReadVariableNames', true);
            tot = [tot;B];
            path4 = [path3,'_round_1_b_',num2str(f),'\','Strat_15_b_',num2str(f),'_from_1_to_700.txt'];
            B = readtable(path4,'ReadVariableNames', true);
            tot = [tot;B];
            path4 = [path3,'_round_2_b_',num2str(f),'\','Strat_15_b_',num2str(f),'_from_1_to_500.txt'];
            B = readtable(path4,'ReadVariableNames', true);
            tot = [tot;B];
        end
    end
    
    if f==1 || f==3 || f==10
        for round=3:8
            path5 = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesMAC\_resultats\_round_',num2str(round)];
            path6 = [path5,'\','Strat_',num2str(strat),'_b_',num2str(f),'_from_1_to_500.txt'];
            B = readtable(path6,'ReadVariableNames', true);
            tot = [tot;B];
        end
    end
    
    
    % We take only 4500 simulations of each (last 4500)
    tot4500 = tot(size(tot,1)-4499:size(tot,1),:);
    
    writetable(tot4500,[pathtot,'Strat_',num2str(strat),'_b_',num2str(f),'_concatenated.txt'])
end

%% Creation des fichiers summary pour chaque strategie
clear all
f=10;
path = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v2\_results\_b_',num2str(f),'\'];
%path = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\_sauvegarde du 23-03-22\_b_',num2str(f),'\'];
sumup_tabAll = [];
path2 = [path,'_histograms\'];
for strat=1:15
    tot = readtable([path,'Strat_',num2str(strat),'_b_',num2str(f),'_concatenated.txt']);
    n=size(tot,1);
    if (1)
        hCt = figure();
        histogram(tot.x1,12); %Ct
        saveas(hCt,[path2,'Strat_',num2str(strat),'_b_',num2str(f),'_hist_Ct_n_',num2str(n),'.png'])
        hNg = figure();
        histogram(tot.x2,12) %Ng
        saveas(hNg,[path2,'Strat_',num2str(strat),'_b_',num2str(f),'_hist_Ng_n_',num2str(n),'.png'])
        hHIV = figure();
        histogram(tot.x3,12) %HIV
        saveas(hHIV,[path2,'Strat_',num2str(strat),'_b_',num2str(f),'_hist_HIV_n_',num2str(n),'.png'])
        hS = figure();
        histogram(tot.x4,12) %syph
        saveas(hS,[path2,'Strat_',num2str(strat),'_b_',num2str(f),'_hist_S_n_',num2str(n),'.png'])
    end
    close all;
    
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
    numStrat  = [strat;strat;strat;strat];
    sumup_tab = table(numStrat,Infection,Moyenne,Mediane,Std,N,Min,Max,low_95,high_95,low_90,high_90);
    writetable(sumup_tab,[path,'Strat_',num2str(strat),'_b_',num2str(f),'_summary.txt'])
    
    sumup_tabAll = [sumup_tabAll;sumup_tab];
end
writetable(sumup_tabAll,[path,'AllStrat_b_',num2str(f),'_summary.txt'])

%% Reshape for LaTex
clear all;
f=10;
path = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v2\_results\_b_',num2str(f),'\'];
tab  = readtable([path,'AllStrat_b_',num2str(f),'_summary.txt']);

pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
fileID = fopen([pathW,'recap-b-',num2str(f),'.txt'],'w');
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
f=10;
createStatTests(f,"Student")
createStatTests(f,"MWW")
createStatTests(f,"KS")
createStatTests(f,"StudentMatlab")

%% Table for LaTex - matrix of p-values
clear all;
f=10;
for infection = {'Ct','Ng','HIV','S'}
    for test={'KS','MWW','StudentMatlab','Student'}
        path = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v2\_results\_b_',num2str(f),'\'];
        tab = readtable([path,'_',test{:},'_p_values_',infection{:},'_b_',num2str(f),'.txt']);
        pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
        fileID = fopen([pathW,'pvalues-b-',num2str(f),'-',test{:},'-',infection{:},'.txt'],'w');
        
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
for f=[1,3,10]
    close all
    %maxNBins=15;
    minBinWidth=0.005;
    figure(1);
    
    %path  = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\_sauvegarde du 23-03-22\_b',num2str(f),'\'];
    path = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v2\_results\_b_',num2str(f),'\'];
    
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
        path2 = [path,'_histograms\'];
        stratAfile = readtable([path,'Strat_',num2str(stratA),'_b_',num2str(f),'_concatenated.txt']);
        for stratB=stratA:nStrat
            plotPos = [round(a+(stratB-1)*(a+Ls),3),round(-2*b + hf-b-hs-(stratA-1)*(hs+b),3), Ls hs];
            if stratA==stratB
                
                %if(0)
                fig1 = figure(1);
                fig1.Position = [10 10 Lf hf];
                %h1 = subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                h1 = axes('Position', plotPos);
                %nBins = min(length([min(stratAfile.x1):0.002:max(stratAfile.x1)]),maxNBins);
                BinWidth=min([minBinWidth,(max(stratAfile.x1)-min(stratAfile.x1))./10]);
                histogram(stratAfile.x1,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
                %minCt = min(minCt,min(stratAfile.x1));
                %maxCt = max(maxCt,max(stratAfile.x1));
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
                %nBins = min(length([min(stratAfile.x2):0.002:max(stratAfile.x2)]),maxNBins);
                BinWidth=min([minBinWidth,(max(stratAfile.x2)-min(stratAfile.x2))./10]);
                histogram(stratAfile.x2,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
                minNg = min(minNg,min(stratAfile.x2));
                maxNg = max(maxNg,max(stratAfile.x2));
                %set(gca,'xtick',[])
                set(h2,'Unit','pixels','Position',plotPos);
                set(gca,'ytick',[])
                set(gca,'fontsize',4)
                %end
                
                fig3 = figure(3);
                fig3.Position = [10 10 Lf hf];
                %h3 = subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                h3 = axes('Position', plotPos);
                %nBins = min(length([min(stratAfile.x3):0.002:max(stratAfile.x3)]),maxNBins);
                BinWidth=min([minBinWidth,(max(stratAfile.x3)-min(stratAfile.x3))./10]);
                histogram(stratAfile.x3,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
                minHIV = min(minHIV,min(stratAfile.x3));
                maxHIV = max(maxHIV,max(stratAfile.x3));
                %set(gca,'xtick',[])
                set(h3,'Unit','pixels','Position',plotPos);
                set(gca,'ytick',[])
                set(gca,'fontsize',4)
                
                
                fig4 = figure(4);
                fig4.Position = [10 10 Lf hf];
                %h4=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                h4=axes('Position', plotPos);
                %nBins = min(length([min(stratAfile.x4):0.002:max(stratAfile.x4)]),maxNBins);
                BinWidth=min([minBinWidth,(max(stratAfile.x4)-min(stratAfile.x4))./10]);
                histogram(stratAfile.x4,'FaceColor','black','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
                minS = min(minS,min(stratAfile.x4));
                maxS = max(maxS,max(stratAfile.x4));
                %set(gca,'xtick',[])
                set(h4,'Unit','pixels','Position',plotPos);
                set(gca,'ytick',[])
                set(gca,'fontsize',4)
                
            else
                stratBfile = readtable([path,'Strat_',num2str(stratB),'_b_',num2str(f),'_concatenated.txt']);
                %Ct
                figure(1);
                %h5=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                h5 = axes('Position', plotPos);
                %nBins = min(length([min(stratAfile.x1):0.002:max(stratAfile.x1)]),maxNBins);
                %BinWidth=max([minBinWidth,min([(max(stratBfile.x1)-min(stratBfile.x1))./10,(max(stratAfile.x1)-min(stratAfile.x1))./10])]);
                BinWidth= (max([stratBfile.x1;stratAfile.x1])-min([stratAfile.x1;stratBfile.x1]))/12;
                histogram(h5,stratAfile.x1,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
                hold on
                %nBins = min(length([min(stratBfile.x1):0.002:max(stratBfile.x1)]),maxNBins,'Normalization','probability');
                histogram(h5,stratBfile.x1,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth)
                %axis off
                set(gca,'xtick',[])
                set(gca,'ytick',[])
                set(h5,'Unit','pixels','Position',plotPos);
                fprintf(['x=',num2str(a+(stratB-1)*(a+Ls)), ' y=',num2str(hf-b-hs-(stratA-1)*(hs+b)),'\n'])
                %h=gca; h.XAxis.TickLength = [0 0];
                %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                %saveas(hCt,[path2,'hist_Ct_Strat_',num2str(stratA),'(',num2str(length(stratAfile.x1)),')_vs_',...
                %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.x1)),')',...
                %    '_b_',num2str(f),'.png'])
                
                figure(2);
                %h6=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                h6 = axes('Position',plotPos);
                %nBins = min(length([min(stratAfile.x2):0.002:max(stratAfile.x2)]),maxNBins);
                %BinWidth=min([minBinWidth,(max(stratBfile.x2)-min(stratBfile.x2))./10,(max(stratAfile.x2)-min(stratAfile.x2))./10]);
                BinWidth= (max([stratBfile.x2;stratAfile.x2])-min([stratAfile.x2;stratBfile.x2]))/12;
                histogram(stratAfile.x2,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
                hold on
                %nBins = min(length([min(stratBfile.x2):0.002:max(stratBfile.x2)]),maxNBins);
                histogram(stratBfile.x2,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth);
                %axis off
                set(h6,'Unit','pixels','Position',plotPos);
                set(gca,'xtick',[])
                set(gca,'ytick',[])
                %set(gca,'xtick',[])
                %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                %saveas(hNg,[path2,'hist_Ng_Strat_',num2str(stratA),'(',num2str(length(stratAfile.x2)),')_vs_',...
                %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.x2)),')',...
                %    '_b_',num2str(f),'.png'])
                
                
                figure(3);
                %h7=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                h7=axes('Position', plotPos);
                %BinWidth=min([minBinWidth,(max(stratBfile.x3)-min(stratBfile.x3))./10,(max(stratAfile.x3)-min(stratAfile.x3))./10]);
                BinWidth= (max([stratBfile.x3;stratAfile.x3])-min([stratAfile.x3;stratBfile.x3]))/12;
                %nBins = min(length([min(stratAfile.x3):0.002:max(stratAfile.x3)]),maxNBins);
                histogram(h7,stratAfile.x3,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
                hold on
                %nBins = min(length([min(stratBfile.x3):0.002:max(stratBfile.x3)]),maxNBins);
                histogram(h7,stratBfile.x3,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth);
                %axis off
                set(gca,'XTick',[])
                set(gca,'ytick',[])
                set(h7,'Unit','pixels','Position',plotPos);
                %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                %saveas(hHIV,[path2,'hist_HIV_Strat_',num2str(stratA),'(',num2str(length(stratAfile.x3)),')_vs_',...
                %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.x3)),')',...
                %    '_b_',num2str(f),'.png'])
                
                
                figure(4);
                %h8=subplot(nStrat,nStrat,(stratA-1)*nStrat+stratB);
                h8 = axes('Position', plotPos);
                %BinWidth=min([minBinWidth,(max(stratBfile.x4)-min(stratBfile.x4))./10,(max(stratAfile.x4)-min(stratAfile.x4))./10]);
                BinWidth=(max([stratBfile.x4;stratAfile.x4])-min([stratAfile.x4;stratBfile.x4]))/12;
                %nBins = min(length([min(stratAfile.x4):0.002:max(stratAfile.x4)]),maxNBins);
                histogram(stratAfile.x4,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
                hold on
                %nBins = min(length([min(stratBfile.x4):0.002:max(stratBfile.x4)]),maxNBins);
                histogram(stratBfile.x4,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth);
                %axis off
                set(gca,'xtick',[])
                set(gca,'ytick',[])
                set(h8,'Unit','pixels','Position',plotPos);
                %legend(['Strat ',num2str(stratA)],['Strat ',num2str(stratB)])
                %saveas(hS,[path2,'hist_S_Strat_',num2str(stratA),'(',num2str(length(stratAfile.x4)),')_vs_',...
                %    'Strat_',num2str(stratB),'(',num2str(length(stratBfile.x4)),')',...
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
            x1 = (k*a+(k-1)*Ls )/Lf;
            y1 = 1-(hf-b-hs-14*(hs+b))/hf;
            annotation('textbox', [x1 y1 Ls/Lf hs/hf], 'string', num2str(k),...
                'VerticalAlignment','middle','HorizontalAlignment','center',...
                'LineStyle','none')
            x2 = 1-2*a/Lf;
            y2 = 1-(4*b+b/3+k*b+(k-1)*hs)/hf;
            annotation('textbox', [x2 y2 Ls/Lf hs/hf], 'string', num2str(k),...
                'VerticalAlignment','middle','HorizontalAlignment','center',...
                'LineStyle','none')
        end
    end
    
    pathRecap = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\graphes\';
    
    %Ct
    saveas(fig1,[path2,'hist_Ct_b_',num2str(f),'.png'])
    set(fig1,'PaperOrientation','landscape');
    saveas(fig1,[pathRecap,'hist_Ct_b_',num2str(f),'.pdf'])
    
    %Ng
    saveas(fig2,[path2,'hist_Ng_b_',num2str(f),'.png'])
    set(fig2,'PaperOrientation','landscape');
    saveas(fig2,[pathRecap,'hist_Ng_b_',num2str(f),'.pdf'])
    
    %HIV
    saveas(fig3,[path2,'hist_HIV_b_',num2str(f),'.png'])
    set(fig3,'PaperOrientation','landscape');
    saveas(fig3,[pathRecap,'hist_HIV_b_',num2str(f),'.pdf'])
    
    %Syphilis
    saveas(fig4,[path2,'hist_S_b_',num2str(f),'.png'])
    set(fig4,'PaperOrientation','landscape');
    saveas(fig4,[pathRecap,'hist_S_b_',num2str(f),'.pdf'])    
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
% 
close all; clearvars -except tabCorresp 
p=0.05; %p-value threshold
f=3; test='KS';

path = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v2\_results\_b_',num2str(f),'\'];
pCt_file = readtable([path,'_',test,'_p_values_Ct_b_',num2str(f),'.txt']);
pNg_file = readtable([path,'_',test,'_p_values_Ng_b_',num2str(f),'.txt']);
pHIV_file = readtable([path,'_',test,'_p_values_HIV_b_',num2str(f),'.txt']);
pS_file = readtable([path,'_',test,'_p_values_S_b_',num2str(f),'.txt']);

pCt  = table2array(pCt_file(1:15,1:15))>p;
pNg  = table2array(pNg_file(1:15,1:15))>p;
pHIV = table2array(pHIV_file(1:15,1:15))>p;
pS   = table2array(pS_file(1:15,1:15))>p;

%Candidate groups of strategy (each line)
listStrat = 1:15;
for k=1:15
    gpCt{k} = listStrat(pCt(k,1:15));
    gpNg{k} = listStrat(pNg(k,1:15));
    gpHIV{k} = listStrat(pHIV(k,1:15));
    gpS{k} = listStrat(pS(k,1:15));
end

matP(:,:,1) = pCt;
matP(:,:,2) = pNg;
matP(:,:,3) = pHIV;
matP(:,:,4) = pS;

grps{1,:} = gpCt;
grps{2,:} = gpNg;
grps{3,:} = gpHIV;
grps{4,:} = gpS;

% Writing in LaTeX file
test='KS';
pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
for infection={'Ct','Ng','HIV','S'}
    fileID = fopen([pathW,'resTest-b-',num2str(f),'-',test,'-',infection{:},'.txt'],'w');
    numInfs=1:4; numInf = numInfs(strcmp(infection{:},{'Ct','Ng','HIV','S'}));
    pInf = matP(:,:,numInf);
    
    %Initialization
    begin = {'\begin{table}[h]','\centering', ['\caption{',test,' - b=',num2str(f),' - ',infection{:},' - pvalue threshold=',num2str(p),'}'],...
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
infections={'Ct','Ng','HIV','S'}; strategies=1:15;
p=0.05; %p-value threshold
test='KS';

for f=[1,3,10]
    for test={'KS','MWW','Student'}
        i=1;
        path = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v2\_results\_b_',num2str(f),'\'];
        pCt_file = readtable([path,'_',test{:},'_p_values_Ct_b_',num2str(f),'.txt']); pMat(:,:,1)=table2array(pCt_file(:,1:15));
        pNg_file = readtable([path,'_',test{:},'_p_values_Ng_b_',num2str(f),'.txt']); pMat(:,:,2)=table2array(pNg_file(:,1:15));
        pHIV_file = readtable([path,'_',test{:},'_p_values_HIV_b_',num2str(f),'.txt']); pMat(:,:,3)=table2array(pHIV_file(:,1:15));
        pS_file = readtable([path,'_',test{:},'_p_values_S_b_',num2str(f),'.txt']); pMat(:,:,4)=table2array(pS_file(:,1:15));
        
        for numInf=1:4
            temp = (tabCorresp(:,:,numInf) - (pMat(:,:,numInf)>p));
            %if 1, then the difference has been found significative whereas it
            %should'nt be ((1-0).
            
            [x,y]=find(temp==1);
            k=1;
            while k<=length(x)/2
                var_b(i) = f;
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
tab = table(var_b',var_infection',var_test',var_strat1',var_strat2',var_pvalue');

% For LaTex / conflict table
pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
fileID = fopen([pathW,'conflitctsTable.txt'],'w');

%Initialization
begin = {'\begin{table}[h]','\centering', ...
    '\begin{tabular}{|l|l|l|l|l|l|}','\hline'};
for k=1:length(begin)
    fprintf(fileID,'%12s\r\n',begin{k});
end
fprintf(fileID,'%12s\r\n','b & infection & test & strat & strat & p-value \\ \hline \hline');
    
for i=1:size(tab,1)    
    fprintf(fileID,'%4s',[num2str(var_b(i)),' & ']);
    fprintf(fileID,'%10s',[var_infection{i},' & ']);
    fprintf(fileID,'%10s',[var_test{i}{:},' & ']);
    fprintf(fileID,'%4s',[num2str(var_strat1(i)),' & ']);
    fprintf(fileID,'%4s',[num2str(var_strat2(i)),' & ']);
    fprintf(fileID,'%5s \r\n',[num2str(var_pvalue(i)),' \\ \hline']);
end

fprintf(fileID,'%12s\r\n','\end{tabular}');
fprintf(fileID,'%12s\r\n','\end{table}');

fclose(fileID);
    
%% Analysis by group (recap-per-group-b-**.txt)
clear all; close all;
f=10;
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
path = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v2\_results\_b_',num2str(f),'\'];
i=1; tab=[];

for g=groups{inf}
    tot=[];
    for strat=g{:}
        tot = [tot;readtable([path,'Strat_',num2str(strat),'_b_',num2str(f),'_concatenated.txt'])];
    end
    Infection =       ["Ct";               "Ng";               "HIV";              "Syph"];
    Moyenne   = round([mean(tot.x1);       mean(tot.x2);       mean(tot.x3);       mean(tot.x4)],3);
    Mediane   = round([median(tot.x1);     median(tot.x2);     median(tot.x3);     median(tot.x4)],3);
    Std       = round([std(tot.x1);        std(tot.x2);        std(tot.x3);        std(tot.x4)],3); %standard deviation (t-student)
    N         = round([length(tot.x1);     length(tot.x2);     length(tot.x3);     length(tot.x4)],3);
    Min       = round([min(tot.x1);        min(tot.x2);        min(tot.x3);        min(tot.x4)],3);
    Max       = round([max(tot.x1);        max(tot.x2);        max(tot.x3);        max(tot.x4)],3);
    low_95    = round([prctile(tot.x1,2.5); prctile(tot.x2,2.5); prctile(tot.x3,2.5); prctile(tot.x4,2.5)],3);
    high_95   = round([prctile(tot.x1,97.5); prctile(tot.x2,97.5); prctile(tot.x3,97.5); prctile(tot.x4,97.5)],3);
    low_90    = round([prctile(tot.x1,05); prctile(tot.x2,05); prctile(tot.x3,10); prctile(tot.x4,10)],3);
    high_90   = round([prctile(tot.x1,95); prctile(tot.x2,95); prctile(tot.x3,90); prctile(tot.x4,90)],3);
    group = string(repmat(['[',num2str(g{:}),']'],4,1));
    sumup_tab = table(group,Infection,Moyenne,Mediane,Std,N,Min,Max,low_95,high_95,low_90,high_90);
    tab = [tab;sumup_tab];
    i=i+1; 
end

% For LaTex
infections = {'Ct','Ng','HIV','S'};

pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
fileID = fopen([pathW,'recap-per-group-b-',num2str(f),'-',num2str(infections{inf}),'.txt'],'w');
begin = {'\begin{longtable}[H]{|l|l|l|l|l|l|l|}',['\caption{groups based on ',infections{inf},' analysis b=',num2str(f),' } \\'],'\hline',...
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

% Analysis by group  (short version) (recap-per-group-b-**-**-short.txt)
clearvars -except groups f
close all;
for inf=1:4
    path = ['C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\Article2_v2\_results\_b_',num2str(f),'\'];

    i=1; tab=[];
    
    for g=groups{inf}
        tot=[];
        for strat=g{:}
            tot = [tot;readtable([path,'Strat_',num2str(strat),'_b_',num2str(f),'_concatenated.txt'])];
        end
        Infection =       ["Ct";               "Ng";               "HIV";              "Syph"];
        Moyenne   = round([mean(tot.x1);       mean(tot.x2);       mean(tot.x3);       mean(tot.x4)],3);
        Mediane   = round([median(tot.x1);     median(tot.x2);     median(tot.x3);     median(tot.x4)],3);
        N         = round([length(tot.x1);     length(tot.x2);     length(tot.x3);     length(tot.x4)],3);
        low_95    = round([prctile(tot.x1,2.5); prctile(tot.x2,2.5); prctile(tot.x3,2.5); prctile(tot.x4,2.5)],3);
        high_95   = round([prctile(tot.x1,97.5); prctile(tot.x2,97.5); prctile(tot.x3,97.5); prctile(tot.x4,97.5)],3);
        group = string(repmat(['(',num2str(g{:}),')'],4,1));
        sumup_tab = table(group,Infection,Moyenne,Mediane,N,low_95,high_95);
        tab = [tab;sumup_tab];
        i=i+1;
    end
    
    % For LaTex
    infections = {'Ct','Ng','HIV','S'};
    
    pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
    fileID = fopen([pathW,'recap-per-group-b-',num2str(f),'-',num2str(infections{inf}),'-short.txt'],'w');
    begin = {'\begin{longtable}[H]{|l|l|l|}',['\caption{groups based on ',infections{inf},' analysis b=',num2str(f),' } \\'],'\hline',...
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

