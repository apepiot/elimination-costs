function [rhohat,Cval,ES,msg_tot,ampl_c] = findRhohat_kit_v8(paramTab,mu,b,f,paramRho,kit,c,alphas,verbose,mySeed,log_path,GlobalOpt,ampl_c,newmod)
%for any combination of SICTP x SEIIIS x SEIIS^m, this
%function finds the argmax of single diseases models, two-disease
%models,...,N diseases models

msg_tot='0';
Nk = length(kit); %nombre d'infections dans le kit

nk = Nk;
k = indexKit(kit);
includeHIV = ~alphas.elim_h;
if includeHIV && ~contains(k,'h')
    mod = [k,'_','h',k];
else
    mod = [k,'_',k];
end

subKit = kit;
optSolver = GlobalOpt;
optSolver.up_bnd_alpha = 0;
optSolver.seed = mySeed;
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
    [rhohat.(mod),Cval.(mod),ES.(mod),P,msg,ampl_c] = fmaxU_knitro_v7_bis(subKit,mod,paramTab,paramRho,mu,b,c,optSolver,ampl_c,newmod);
    vecRhohat(j) = rhohat.(mod);
    vecCval(j)   = Cval.(mod);
    tabES(j,:)   = ES.(mod);
    
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
            newmod=1;
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
ES.tot     = tabES(vecCval==Cval.tot,:);
end

