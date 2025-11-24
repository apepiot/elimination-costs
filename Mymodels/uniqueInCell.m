function myOutput = uniqueInCell(myCellofCells)

vecElim = [];
for j=2:length(myCellofCells)
    elt = myCellofCells{j};
    for i=1:(j-1)
        elt2 = myCellofCells{i};
        if isequal(elt,elt2)
            vecElim = [vecElim,j];
        end
        
    end
    j=j+1;
end

myOutput = {myCellofCells{setdiff(1:length(myCellofCells),vecElim)}};

end