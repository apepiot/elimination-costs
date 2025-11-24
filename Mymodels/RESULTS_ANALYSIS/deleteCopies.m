function newTab = deleteCopies(tab)
newTab = tab;
rows_to_delete=[];
for i=1:length(tab.IDech)
    currSet = tab(i,:);
    %doublons = tab(tab.IDech==currSet.IDech & tab.p==currSet.p & tab.noEch==currSet.noEch,:);
    indexes =  find(tab.IDech==currSet.IDech & tab.p==currSet.p & tab.noEch==currSet.noEch);
    rows_to_delete = [rows_to_delete,indexes(2:end)'];
    %disp(i)
end

newTab(rows_to_delete,:) = [];    

end

