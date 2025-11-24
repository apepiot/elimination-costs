function [res] = findCnn_inf_mod(temp,k,mod)

res = table(NaN,NaN,NaN,NaN, 'VariableNames', {'h','s', 'c', 'g'});


while ~isempty(k)
    k_mod = [k,'_',mod];
    
    if ismember({[k_mod,'_elim']},temp.Properties.VariableNames)
        inf_elim = temp.([k_mod,'_elim']){:};
        
        if ~isempty(inf_elim) && ~isequal(inf_elim,'0')
            
            for inf=inf_elim
                res.(inf) = temp.(k_mod);
                k   = erase(k,inf);
                mod = erase(mod,inf);
            end
            
        elseif isequal(inf_elim,'0') %one infection is already eliminated without VT
            for inf=mod
                if temp.(['elim_',inf])==1
                    res.(inf) = 0;
                    k   = erase(k,inf);
                    mod = erase(mod,inf);
                end
            end
        else
            k   = '';
        end
    else
        k='';
    end
end

end

