function [newTable] = replaceLineInTable(oldTable,ech)
newTable = oldTable;
VarNames = newTable.Properties.VariableNames;
    %File containing
    for i=1:size(ech,1)
        currSet = ech(i,:);
        ID_ech  = currSet.IDech; %tester cette simu avec hscg, round 9
        pHIV    = currSet.p;
        roundNo = currSet.roundNo;
        f       = currSet.f;
        kit     = currSet.kit{:};
        
        idx = find(newTable.IDech==ID_ech & newTable.p==pHIV & newTable.roundNo==roundNo & ismember(newTable.kit,kit) & newTable.f==f);
        rowToReplace = newTable(idx,:);
        if ~isempty(rowToReplace)
            newTable(idx,VarNames) = currSet(:,VarNames);
        end
        if mod(i,1000)==0
            disp(i)
        end
    end
end

