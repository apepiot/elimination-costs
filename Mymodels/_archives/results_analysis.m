%Analyse des résultats

%path  = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\numSim.txt';
%path2 = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesMAC\_resultats\bestStrat_51_to_100_b_5.txt';
path2 = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesMAC\_resultats\b_1\bestStrat_1_to_1000.txt';

A = readtable(path2);
%[g, T2] = findgroups(T);
%h = splitapply(@numel,T,g);  
%B=table(T2,(h))
%1:'Chlam', 2:'Gono', 3:'HIV', 4:'Syphilis'
uV = unique(A(:,1))

x = size(uV,1);
for i=1:x
    nbOcc = sum(strcmp(uV.bestStrat(i),A.bestStrat));
    uV.n(i) = nbOcc;
end
%1:'Chlam', 2:'Gono', 3:'HIV', 4:'Syphilis'
uV.freq = round(uV.n/sum(uV.n)*100,1)


% N=numel(T.bestStrat);
% for i = 1:N
%     if isequal(T.bestStrat{i},'equality')
%         disp(num2str(i))
%     end
% end

%% with a loop on b
clear all;close all;
T=[];
for sim=[1:25:1476]+3000
%path2 = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesMAC\_resultats\Sim_',num2str(sim),'_to_',num2str(sim+49),'_b_10.txt'];
%path2 = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesMAC\_resultats\Sim_',num2str(sim),'_to_',num2str(sim+99),'_b_10.txt'];
path2 = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\Sim_',num2str(sim),'_to_',num2str(sim+24),'_b_10.txt'];

A = readtable(path2);
T = [T;A];
end
%for sim=[2101:50:2551]%[1601:50:2051,3101:50:3551,4601:50:5051] :1500 sim for b=10 %[1101:50:1551,2601:50:3051,4101:50:4551] : 1500 simulations for b=3
%    path2 = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesMAC\_resultats\Sim_',num2str(sim),'_to_',num2str(sim+49),'_b_1.txt'];
%    A = readtable(path2);
%    T = [T;A];
%end

%T.Properties.VariableNames = ["strategy","Ct","Ng","HIV","Syphilis","numSim"];
T.Properties.VariableNames = ["strategy","HIV","Syphilis","Ct","Ng","numSim"];

numSims = unique(T.numSim);
t=1;
for i=numSims'
    currentSim = T(T.numSim==i,:);
    costStratTable_sorted = sortrows(currentSim,[{'HIV'},{'Syphilis'},{'Ct'},{'Ng'}],'descend');
    costStratTable_sorted = costStratTable_sorted(:,[{'strategy','HIV'},{'Syphilis'},{'Ct'},{'Ng'}]);

    bestStrategy.strat(t) = costStratTable_sorted{1,1};
    bestStrategy.numSim(t) = i;
    t=t+1;
end
%1:'Chlam', 2:'Gono', 3:'HIV', 4:'Syphilis'
%uV = myUnique(A(:,1))'
uV(:,1)=unique(bestStrategy.strat);
x = size(uV,1);
for i=1:x
    nbOcc = sum(strcmp(uV(i,1),bestStrategy.strat),'omitnan');
    uV{i,2} = nbOcc;
end
%1:'Chlam', 2:'Gono', 3:'HIV', 4:'Syphilis'
uV
round([uV{:,2}]'/size(numSims,1)*100,1)
