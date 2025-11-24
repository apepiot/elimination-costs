function [rhohat,Cval,ES,msg,ampl_c] = findRhohat_kit_v9(paramTab,mu,b,f,paramRho,kit,c,alphas,verbose,mySeed,log_path,GlobalOpt,ampl_c,newmod,ampl_models_dir)
%for any combination of SICTP x SEIIIS x SEIIS^m, this
%function finds the argmax of single diseases models, two-disease
%models,...,N diseases models
nbRelanceMax = 10; noRelance=0;
Nk = length(kit); %nombre d'infections dans le kit
msg = '0';
nk = Nk;
k = indexKit(kit);
includeHIV = ~alphas.elim_h;
if includeHIV && ~contains(k,'h')
    mod = [k,'_','h',k];
else
    mod = [k,'_',k];
end

subKit = kit;
subMod = mod;
optSolver = GlobalOpt;
optSolver.up_bnd_alpha = 0;
optSolver.seed = mySeed;
j=1; vecRhohat=[]; vecCval=[];
infsElims = alphas.([subMod,'_elim']);
typeofInfsElim = class(infsElims);
if isequal(typeofInfsElim,'cell')
    if isempty(infsElims{:})
        alphas.(subMod)=GlobalOpt.sup_bnd_alpha;
    end
else
    if isnan(infsElims)
        alphas.(subMod)=GlobalOpt.sup_bnd_alpha;
        alphas.([subMod,'_elim'])={''};
    elseif infsElims==0
        alphas.([subMod,'_elim'])={''};
    end
end
optSolver.inf_bnd_alpha = optSolver.up_bnd_alpha;
while nk>0 && noRelance<nbRelanceMax
    %optSolver.inf_bnd_alpha = optSolver.up_bnd_alpha; %#
    if ismember(subMod, alphas.Properties.VariableNames)
        optSolver.up_bnd_alpha = alphas.(subMod);
    else
        optSolver.up_bnd_alpha = GlobalOpt.sup_bnd_alpha;
    end
    if ~contains(subMod,'h')
        paramTab{3}.p=0;
    end
    
    if optSolver.up_bnd_alpha<optSolver.inf_bnd_alpha
        ES.tot     = nan(560,1);
        Cval.tot   = NaN;
        rhohat.tot = NaN;
        msg    = '299';
        warning('up_bnd_alpha<inf_bnd_alpha')
        break
    end
    
    [rhohat.(subMod),Cval.(subMod),ES.(subMod),~,msg,ampl_c] = fmaxU_knitro_v8(kit,subKit,mod,subMod,paramTab,paramRho,f,mu,b,c,optSolver,ampl_c,newmod,ampl_models_dir);
    
    if str2double(msg)==0
        vecRhohat(j) = rhohat.(subMod);
        vecCval(j)   = Cval.(subMod);
        tabES(j,:)   = ES.(subMod);
        
        %Le nouveau modele depend de l'infection qui a ete eliminee.
        if ismember([subMod,'_elim'], alphas.Properties.VariableNames)
            infsElims = alphas.([subMod,'_elim']);
        else
            infsElims = {''};
        end
        
        if ~isempty(infsElims{:})
            nk = nk-length(infsElims{:});  %on retire au modele la ou les infections eliminee a alpha
            for i=1:length(infsElims{:})
                subMod = erase(subMod,infsElims{1}(i));
                k = erase(k,infsElims{1}(i));
                subKit = kToKit(k);
            end
        else
            nk=0;
        end
        j=j+1;
        optSolver.inf_bnd_alpha = optSolver.up_bnd_alpha;
    else
        optSolver.seed = optSolver.seed+1;
        noRelance      = noRelance+1;
        disp('       relance dans findRhohat_v9')
    end
    newmod=0;
end

if exist('tabES','var')
    Cval.tot   = min(vecCval);
    rhohat.tot = unique(vecRhohat(vecCval==Cval.tot));
    EStot      = tabES(vecCval==Cval.tot,:); %pour eviter les doublons
    ES.tot     = EStot(1,:);
else
    ES.tot     = nan(560,1);
    Cval.tot   = NaN;
    rhohat.tot = NaN;
end

end

