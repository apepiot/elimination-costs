function [costElimTable,needToLowerC] = createCostElimTable(tab,tabcn,vecAlphas)
needToLowerC=0;
vecC = tab.c;
costElimTable = table('Size',[4^2-1,5], 'VariableNames',{'model', 'Chlam', 'Gono', 'HIV', 'Syphilis'},...
    'VariableTypes', {'string','double','double','double','double'});

costElimTable = standardizeMissing(costElimTable,0);

%1 disease model
for i=1:4
    costElimTable{i,1} = {num2str(i)};
    costElimTable{i,i+1} = round(tabcn.one(i),10,"significant");
end

%2 disease model
rhohat2d_reshape = reshape([tab.two(:).rhohat],[],6);
for k=1:6
    nameModel = num2str(tab.duos(k,:));
    costElimTable{k+4,1} = {nameModel(find(~isspace(nameModel)))};
    
    for dis=tab.duos(k,:)
        j=max(find(rhohat2d_reshape(:,k)>=vecAlphas(dis)));
        if ~isempty(j)
            costElimTable{k+4,dis+1} = round(vecC(j),10,"significant");
        else
            needToLowerC = 1;
        end
    end
end

%3 disease model
rhohat3d_reshape = reshape([tab.three(:).rhohat],[],4);

for k=1:4
    nameModel = num2str(tab.trios(k,:));
    costElimTable{k+4+6,1} = {nameModel(find(~isspace(nameModel)))};
    
    for dis=tab.trios(k,:)
        j=max(find(rhohat3d_reshape(:,k)>=vecAlphas(dis)));
        if ~isempty(j)
            costElimTable{k+4+6,dis+1} = round(vecC(j),10,"significant");
            else
            needToLowerC = 1;
        end
    end
end

%4 disease model
nameModel = num2str('1234');
costElimTable{4+6+4+1,1} = {nameModel};
for dis=1:4
    j=max(find(tab.rhohat.rhohat>=vecAlphas(dis)));
    if ~isempty(j)
        costElimTable{1+4+6+4,dis+1} = round(vecC(j),10,"significant");
    else
            needToLowerC = 1;
    end
end
end

