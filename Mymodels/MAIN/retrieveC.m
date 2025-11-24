function [res] = retrieveC(tabCosts,alphas,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

res = [];
nk  = length(k);
includeHIV = ~alphas.elim_h;
if includeHIV && ~contains(k,'h')
    mod = [k,'_','h',k];
else
    mod = [k,'_',k];
end

subMod = mod;
% while nk>1
    if ismember([subMod,'_elim'], alphas.Properties.VariableNames)
        infsElims = alphas.([subMod,'_elim']);
    else
        infsElims = {''};
    end
    infeliminated = infsElims{:};
    if ~isempty(infeliminated)
        nk = nk-length(infeliminated);
        for i=1:length(infeliminated)
            subMod = erase(subMod,infsElims{1}(i));
            k      = erase(k,infsElims{1}(i));
            subKit = kToKit(k);
        end
%    else
%         nk = 0;
     end
% end

cnew = table2array(tabCosts(tabCosts.kit==string(k),{'HIV','syphilis','Ct','Ng'}));
res  = cnew(~isnan(cnew));

end

