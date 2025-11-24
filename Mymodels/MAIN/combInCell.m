function [res] = combInCell(myCell)
% gives all combinations of elements in the cell
j=1;
n=numel(myCell);
res={};
for k=1:n
    comb_k = nchoosek(1:n,k);
    for i=1:size(comb_k,1)
        res{j} = {myCell{comb_k(i,:)}};
        j=j+1;
    end
end

end

