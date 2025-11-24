function alpha_max = getAlphaMax(alphas,mod,nk)
subMod = mod;
alpha_max = alphas.(mod);

while nk>=1
    if ismember([subMod,'_elim'], alphas.Properties.VariableNames)
        alpha_max = alphas.(subMod);
        
        %On passe au modele suivant
        infsElims = alphas.([subMod,'_elim']);
        nk = nk-length(infsElims{:});  %on retire au modele la ou les infections eliminee a alpha
        for i=1:length(infsElims{:})
            subMod = erase(subMod,infsElims{1}(i));
        end
    else
        infsElims = {''};
    end
end
% 
% while ~isempty(infsElims{:}) && nk>1
%     nk = nk-length(infsElims{:});  %on retire au modele la ou les infections eliminee a alpha
%     
%     for i=1:length(infsElims{:})
%         subMod = erase(subMod,infsElims{1}(i));
%         alpha_max = alphas.(subMod);
%     end
%     
%     if nk>1
%         if ismember([subMod,'_elim'], alphas.Properties.VariableNames)
%             infsElims = alphas.([subMod,'_elim']);
%         else
%             infsElims = {''};
%         end
%     end
% end

end

