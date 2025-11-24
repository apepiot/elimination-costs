function tabAlpha = createTabAlpha(list_kits_to_consider,n)

VarNames = {'IDech','roundNo','nbEch','p','timeCompil'};
j = length(VarNames);
for i=1:length(list_kits_to_consider)
    kit=list_kits_to_consider(i);
    %disp(kit{:});
    k = indexKit(kit{:});
    VarNames{i+j} = [k,'_',k];
    if ~contains(k,'h')
        VarNames{i+j+1} = [k,'_h',k];
        j=j+1;
    end
    
end
m=length(list_kits_to_consider)+j;
VarNames{m+1} = 'msg';

tabAlpha = table('Size',[n,m+1],...
    'VariableNames',VarNames, 'VariableTypes',[repmat({'double'},1,m),'string']);
end

