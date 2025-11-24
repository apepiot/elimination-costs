function res = isCharInCell(cellArray,myChar)
    res=zeros(1,length(cellArray));
    for i=1:length(cellArray)
        %res(i) = ismember(myChar,cellArray{i}); 
        res(i) = contains(cellArray{i},myChar); 
    end
    
    res = logical(res);
end

