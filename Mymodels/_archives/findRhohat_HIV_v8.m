function [rhohat,Cval,ES,msg_tot] = findRhohat_HIV_v8(paramTab,mu,b,f,paramRho,kit,c,alphas,verbose,mySeed,log_path,GlobalOpt)
%for any combination of SICTP x SEIIIS x SEIIS^m, this
%function finds hat_rho_k such that the prevalence of HIV is eliminated

%% a retravailler
msg_tot = '0';
includeHIV = ~alphas.elim_h;
Nk = length(kit); %nombre d'infections dans le kit

nk = Nk;
k = indexKit(kit);
if includeHIV && ~contains(k,'h')
    mod = [k,'_','h',k];
    error('the kit does not contain HIV')
else
    mod = [k,'_',k];
end

GlobalOpt.up_bnd_alpha = 0;
j=1; vecRhohat=[]; vecCval=[];
infsElims = alphas.([mod,'_elim']);
typeofInfsElim = class(infsElims);
if isequal(typeofInfsElim,'cell')
    if isempty(infsElims{:})
        alphas.(mod)=GlobalOpt.sup_bnd_alpha;
    end
else
    if isnan(infsElims)
        alphas.(mod)=GlobalOpt.sup_bnd_alpha;
        alphas.([mod,'_elim'])={''};
    elseif infsElims==0
        alphas.([mod,'_elim'])={''};
    end
end
while nk>0
    optSolver.inf_bnd_alpha = optSolver.up_bnd_alpha;
    if ismember(mod, alphas.Properties.VariableNames)
        optSolver.up_bnd_alpha = alphas.(mod);
    else
        optSolver.up_bnd_alpha = GlobalOpt.sup_bnd_alpha;
    end
    if ~contains(mod,'h')
        paramTab{3}.p=0;
    end
    [rhohat.(mod),Cval.(mod),ES,P,msg] = fminU_knitro_v7_bis(subKit,mod,paramTab,paramRho,mu,b,c,optSolver);
    vecRhohat(j) = rhohat.(mod);
    vecCval(j) = Cval.(mod);
    
    %Le nouveau modele depend de l'infection qui a ete eliminee.
    if ismember([mod,'_elim'], alphas.Properties.VariableNames)
        infsElims = alphas.([mod,'_elim']);
    else
        infsElims = {''};
    end
    
    if ~isempty(infsElims{:})
        nk = nk-length(infsElims{:});    
        for i=1:length(infsElims{:})
            mod = erase(mod,infsElims{1}(i));
            k = erase(k,infsElims{1}(i));
            subKit = kToKit(k);
        end
    else
        nk=0;
    end
    j=j+1;
    
    if ~isequal(msg,'0')
        msg_tot = msg;
    end
end
Cval.tot   = min(vecCval);
rhohat.tot = vecRhohat(vecCval==Cval.tot);
end

