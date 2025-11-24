function [bestStrategy,worstStrategy,costStratTable_sorted,costStratTable] = createCostElimStrat(costElimTable)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Odering strategies by cost
%e.g. 1x2 + 3 + 4 ou 1x2 + 3x4

costStratTable = array2table(NaN(15,5),'VariableNames',{'strategies', 'Chlam', 'Gono', 'HIV', 'Syphilis'});
% costStratTable = table('Size',NaN(15,5), 'VariableNames',{'strategies', 'Chlam', 'Gono', 'HIV', 'Syphilis'},...
%     'VariableTypes', {'string','double','double','double','double'});
costStratTable.strategies = {{'1','2','3','4'}, {'1', '2', '34'}, {'1', '234'},...
    {'1', '23', '4'}, {'1', '24', '3'},...
    {'12', '3', '4'}, {'12', '34'},...
    {'13', '2', '4'}, {'13', '24'},...
    {'14', '2', '3'}, {'14', '23'},...
    {'123','4'},{'124','3'},{'134','2'},...
    {'1234'}}';
strategies = costStratTable.strategies;

%costs associated to each strategy
nbStrat = size(strategies,1);
for i=1:nbStrat
    nbKits = size(strategies{i},2);
    for j=1:nbKits
        kitCost = costElimTable{costElimTable.model == strategies{i}{j},2:5};
        costStratTable{i,2:5} = sum([kitCost; costStratTable{i,2:5}],1,'omitnan') ;
    end
end

%detailing strategies
for i=1:nbStrat
    costStratTable.strategiesdetailed{i} = strjoin(costStratTable.strategies{i},'+');
    %costStratTable.strategiesdetailed{i};
end

%% sorting by costs (HIV, syphilis, chlamydia, gonorrhea)
costStratTable_sorted = sortrows(costStratTable,[{'HIV'},{'Syphilis'},{'Chlam'},{'Gono'}],'descend');
costStratTable_sorted = costStratTable_sorted(:,[{'strategiesdetailed','HIV'},{'Syphilis'},{'Chlam'},{'Gono'}]);

%writetable(costStratTable_sorted,'myData.txt','Delimiter',' ')

bestStrategy = costStratTable_sorted{1,1};

if(sum(table2array(costStratTable_sorted(1,2:5))==table2array(costStratTable_sorted(2,2:5)))==4)
    bestStrategy = {'equality '};
end

worstStrategy = costStratTable_sorted{nbStrat,1};
end

