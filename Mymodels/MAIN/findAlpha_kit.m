function [alpha,P,ES,msg,elim] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit,mod,afficherOutput,log_path,paramSolver,ampl_models_dir)
k = indexKit(kit);
if ~contains(k,'h')
    mod = ['h',k]; k_mod = [k,'_',mod];
else
    mod = k; k_mod = [k,'_',k];
end

elim = []; msg =[];
k_suivant       = k;
k_mod_suivant   = k_mod;
kit_suivant     = kToKit(k_suivant);
disp(k_mod_suivant);


while ~isempty(k_suivant)
    [alpha.(k_mod_suivant),P,ES,msgAlpha_mod,elim_i] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit_suivant,k_mod_suivant,afficherOutput,...
        log_path,paramSolver,ampl_models_dir);
    alpha.k = alpha.(k_mod_suivant);
    msg                  = [msg, ' ',msgAlpha_mod];
    elim.(k_mod_suivant) = elim_i;
    alpha.([k_mod_suivant,'_elim']) = {elim_i};
    for inf=elim_i
        k_suivant       = erase(k_suivant,inf);
        k_mod_suivant   = erase(k_mod_suivant,inf);
        kit_suivant     = kToKit(k_suivant);
    end
end
end

