% Results analysis of RHOHAT_v6
close all; clear all;
folderRun = 'Run5';
path  = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\'];
[status, msg, msgID] = mkdir([path,'\_results']);

%% Summary table creation
clear all; close all;
%The following section creates the concataned .txt files of the cost of
%elimination of each infection for each combination of [b,strategy,pHIV].
%-------------------------%
folderRun = 'Run5';
roundNos = 1;
stratNos = 1:15;
vecF     = [1];
vecpHIV  = [0.50];
nbSimPerRound = 50; %faire en sorte de lire n'importe quel fichier
%-------------------------%
path  = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',folderRun,'\'];

for f=vecF
    for strat=stratNos
        for pHIV=vecpHIV
            B = table('Size',[0,5],'VariableNames', ["Ct" "Ng","HIV","syph","IDsim"], 'VariableTypes', ["double","double","double","double","double"]);
            for roundNo=roundNos
                file = ['_round_',num2str(roundNo),'\strat_',num2str(strat),'_b_',num2str(f),...
                    '_p_',num2str(round(pHIV*100)),'_',num2str("1"),'-',num2str(nbSimPerRound),'_elim_cost.txt'];
                A = readtable([path,file],'ReadVariableNames', true);
                B = [B;A];
            end
            writetable(B,[path,'\_results','\strat_',num2str(strat),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),...
                '_elim_cost','_concatenated.txt'])
        end
    end
end


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

