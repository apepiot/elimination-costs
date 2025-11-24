function list_i = findModInStrats(mod,list)
n = length(list);
list_i = [];
for i=1:n
    if sum(isCharEqCell(list{i},mod))==1
        list_i = [list_i,i];
    end
end
end
