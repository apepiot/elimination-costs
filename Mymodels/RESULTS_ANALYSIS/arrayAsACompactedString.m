function res = arrayAsACompactedString(array)
%array : a vector containing integers

if 0
res = [];
array_sorted = sort(array);
diff = array_sorted(2:end)-array_sorted(1:(end-1));

debut = [array_sorted(1),array_sorted([0==1,diff(1:end)>1])];
fin   = [array_sorted([diff(1:end)>1,0==1]),array_sorted(end)];

for i=1:length(debut)
    res = [res, num2str(debut(i)),'-',num2str(fin(i))]; 
    if i~=length(debut)
        res = [res,'_'];
    end
end
else %short version
    res = [num2str(min(array)),'-',num2str(max(array))];
end


end

