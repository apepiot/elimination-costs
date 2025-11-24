function res = isCharEqCell(cellArray,myChar)
    res=zeros(1,length(cellArray));
    for i=1:length(cellArray)
        res(i) = strcmp(cellArray{i},myChar); 
    end
    
    res = logical(res);
end

