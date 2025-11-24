function [newtab] = sortTable(tab)
tabAsArray = struct2table(tab);
[~,order] = sort(tabAsArray.c);
tabAsArray2 = tabAsArray(order,:);
newtab = table2struct(tabAsArray2,"ToScalar",true);
end

