function [rhohat,Cval,ES,msg_tot] = findRhohat_kit_v8_by_dis(paramTab,mu,b,f,paramRho,kit,bydis,c,alphas,verbose,mySeed,log_path,GlobalOpt)
%for any combination of SICTP x SEIIIS x SEIIS^m, this
%function finds the argmax of single diseases models, two-disease
%models,...,N diseases models

msg_tot='0';
Nk = length(kit); %nombre d'infections dont on cherche le rho_k qui maximise l'utilite
nk = Nk;
k = indexKit(kit);

includeHIV = ~alphas.elim_h;

for inf=bydis
    if ~contains(k,inf)
        error([' ',inf,' is not in the kit ', k])
    end
end


if ~includeHIV && isequal(bydis,'h')
    %faire le cas general pour n importe quelle infection
    warning(' on cherche a maximiser selon h, alors que h est deja eliminé')
    rhohat.h = 0;
    rhohat.tot = 0;
    msg_tot='-0';
    return
end

if includeHIV && ~contains(k,'h')
    mod_obj = [bydis,'_h',k];
    mod     = [k,'_h',k];
else
    mod_obj = [bydis,'_',k];
    mod     = [k,'_',k];
end

subKit = kit;
optSolver = GlobalOpt;
optSolver.up_bnd_alpha = 0;

j=1; vecRhohat=[]; vecCval=[];

% Which infection is eliminated first in the model mod
% infsElims = alphas.([mod,'_elim']); typeofInfsElim = class(infsElims);
% if isequal(typeofInfsElim,'cell')
%     if isempty(infsElims{:})
%         alphas.(mod)=GlobalOpt.sup_bnd_alpha;
%     end
% else
%     if isnan(infsElims)
%         alphas.(mod)=GlobalOpt.sup_bnd_alpha;
%         alphas.([mod,'_elim'])={''};
%     elseif infsElims==0
%         alphas.([mod,'_elim'])={''};
%     end
% end
% 
% optSolver.inf_bnd_alpha = optSolver.up_bnd_alpha;
% if ismember(mod, alphas.Properties.VariableNames)
%     optSolver.up_bnd_alpha = alphas.(mod);
% else
%     optSolver.up_bnd_alpha = GlobalOpt.sup_bnd_alpha;
% end
% if ~contains(mod,'h')
%     paramTab{3}.p=0;
% end
[rhohat.(mod_obj),Cval.(mod_obj),ES.(mod_obj),P,msg] = fminU_knitro_v7_bis(subKit,mod_obj,paramTab,paramRho,mu,b,c,optSolver);
% vecRhohat(j) = rhohat.(mod);
% vecCval(j) = Cval.(mod);
% tabES(j,:) = ES.(mod);
% 
% %Le nouveau modele depend de l'infection qui a ete eliminee.
% if ismember([mod,'_elim'], alphas.Properties.VariableNames)
%     infsElims = alphas.([mod,'_elim']);
% else
%     infsElims = {''};
% end
% 
% if ~isempty(infsElims{:})
%     nk = nk-length(infsElims{:});
%     for i=1:length(infsElims{:})
%         mod = erase(mod,infsElims{1}(i));
%         k = erase(k,infsElims{1}(i));
%         subKit = kToKit(k);
%     end
% else
%     nk=0;
% end
% j=j+1;
% 
% if ~isequal(msg,'0')
%     msg_tot = msg;
% end
% 
% Cval.tot   = min(vecCval);
% rhohat.tot = vecRhohat(vecCval==Cval.tot);
% ES.tot     = tabES(vecCval==Cval.tot,:);



