function tabAlpha = createTabAlpha2(list_kits_main,n)
list_kits_to_consider = {};
for kk=1:length(list_kits_main)
    list_kits_to_consider = [list_kits_to_consider,combInCell(list_kits_main{kk})];
end

VarNames_1 = {'IDech','roundNo','nbEch','p','timeCompil'};
VarNames_2 = {}; VarNames_3={}; j=0;
%j = length(VarNames);
for i=1:length(list_kits_to_consider)
    kit=list_kits_to_consider(i);
    %disp(kit{:});
    k = indexKit(kit{:});
    VarNames_2{i+j} = [k,'_',k];
    VarNames_3{i+j} = [k,'_',k,'_elim'];
    if ~contains(k,'h')
        VarNames_2{i+j+1} = [k,'_h',k];
        VarNames_3{i+j+1} = [k,'_h',k,'_elim'];
        j=j+1;
    end   
end
%m=length(list_kits_to_consider)+j;

VarNames_22 = unique(VarNames_2);
VarNames_33 = unique(VarNames_3);
m = length(VarNames_2);
VarNames_3{m+1} = 'msg';

VarNames = [VarNames_1,VarNames_2,VarNames_3];
tabAlpha = table('Size',[n,5+m+m+1],'VariableNames',VarNames,...
    'VariableTypes',[repmat({'double'},1,5+m),repmat({'string'},1,m+1)]);
end

