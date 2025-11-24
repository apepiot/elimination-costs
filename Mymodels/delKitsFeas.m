function list_kits_to_consider = delKitsFeas(list_kits_to_consider,msgSol)
indexModARetirer=[];
for dis=["HIV","syphilis","Ct","Ng"]
    if isequal(msgSol.(dis),'-2')
        for numKit=1:length(list_kits_to_consider)
            kit = list_kits_to_consider{numKit};
            if any(contains(kit,dis))
                indexModARetirer = [indexModARetirer,numKit];
            end
        end
    end
end
list_kits_to_consider=list_kits_to_consider(setdiff(1:length(list_kits_to_consider),indexModARetirer));
end

