% Code which solve the utlity maximization pb in a n-diseases model (n<=4
% for now)
%DIFFERENCE WITH RHOHAT_v3_3 : bestStrategy has been changed.

% one SICR : HIV
% one SEIIIS : syphilis
% one or two SEIIS : chlamydia, gonorrhea

%xx for each sim
%xx

%% Parameters

clear all;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ct=true; Ng=true; HIV=true; Syph=true; plotFigure=false;
nbSim = 500; sauvegardetousles = 25;

nSEIIS=Ct+Ng; nSICR=HIV; nSEIIIS=Syph;
N = nSICR+nSEIIS+nSEIIIS;
%STIChoice = 'Ct'; STIChoiceMini = 'Ct';%if one STI considered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = 5;
path  = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\numSim.txt';
path2 = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\';
%path = '_resultats/numSim.txt';
%path2 = '_resultats/';
%%
vecCnn=[];vecC0=[];
for biasFactor = 1
    bestStrat={};worstStrat={};tps1=tic;tabRecapStratCurrentSim=[];
    fileID = fopen(path,'r');
    formatSpec = '%d';
    A = fscanf(fileID,formatSpec);
    lastSim = A(end); %a recuperer dans un fichier
    fclose(fileID);
    for i=(lastSim+1):(nbSim+lastSim+1)
        %cond=0;t=0;
        %while(~cond & t<10000)
        i
        [paramTab,maxalpha] = sampleParameters(Ct,Ng,HIV,Syph,b)
        %%
        vecC = linspace(-0.06,0.09,10);
        optFindRhohat.sampleSize = length(vecC);%100;
        optFindRhohat.Tol = 10^(-5);
        [bestStrat{i-lastSim,1},worstStrat(i-lastSim),costStratTable,cmin,cmax,costStratTable_sorted,tab_refined,tab,costElimTable] =...
            findBestStrategy_v2(vecC(1),vecC(end),N,paramTab,mu,b,biasFactor,optFindRhohat);
        bestStrat{i-lastSim,2} = i;
        
        costStratTable_sorted.nbSim = i*ones(15,1);
        tabRecapStratCurrentSim = [tabRecapStratCurrentSim;costStratTable_sorted];
        
        %on sauvegarde les resultats tous les x simulations
        x = sauvegardetousles;
        if (mod(i,x)==0)
            fileID = fopen(path,'w');
            fprintf(fileID,'%d\n',i);
            fclose(fileID);
            fileName = ['Sim_',num2str(i-x+1),'_to_',num2str(i),'_b_',num2str(biasFactor)];
            %fileID2 = fopen(path2,'w');
            %fprintf(fileID2,'%s',tabRecapStratCurrentSim)
            writetable(tabRecapStratCurrentSim,[path2,fileName],'WriteVariableNames',false)
            writetable(cell2table(bestStrat,"VariableNames",["bestStrat","nSim"]),...
                [path2,'bestStrat_',num2str(lastSim+1),'_to_',num2str(i),'_b_',num2str(biasFactor),'.txt'])
            tabRecapStratCurrentSim=[]; % reinitializing to empty
            %bestStrat={};
        end
        vecCnn = [vecCnn,cmin];
        vecC0 = [vecC0,cmax];
    end
    
end
tps2 = toc(tps1);
bestStrat
worstStrat
%
% T=bestStrat;
% [g, T2] = findgroups(T);
% h = splitapply(@numel,T,g);
% B=table((T2'),(h'))
%1:'Chlam', 2:'Gono', 3:'HIV', 4:'Syphilis'
