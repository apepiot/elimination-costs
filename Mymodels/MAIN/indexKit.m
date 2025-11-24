function k = indexKit(kit)
k = []; 

listIndexDis = {'h','s','c','g'};
listDis = {'HIV','syphilis','Ct','Ng'};
for j=1:4
    if ismember(listDis{j},kit)
        k=[k,listIndexDis{j}];
    end
end

end

