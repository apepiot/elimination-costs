function [msgSolver,msgSol] = checkFeasability(paramTab,paramRho,b,mu,f,kit,verbose,myTol)
%Chercher rhohat_k' a-t-il un sens ? Non, si pour rho_k=0, les infections
%sont déjà éliminées
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';

[ES,~,~,~,~,~,msgSolver] = ...
    P1234_SICTPSEIIISSEIIS2_v6(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,f,'knitroampl',20,verbose);
if ~isequal(msgSolver,'0')
    disp('Le solver n a pas convergé')
end

N = length(kit);
for i=1:N
    dis = kit{i};
    Pi = P_kit(ES,{dis});
    if (Pi<myTol)
        disp(['L infection ', dis, ' est deja eliminee sans depistage volontaire']);
        msgSol.(dis) = '-2';
    else
        msgSol.(dis) = '0';
    end
end
end