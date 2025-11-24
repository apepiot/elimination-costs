function [infElim,msgSolver,msgSol] = checkFeasability_v2(paramTab,paramRho,b,mu,f,mod,verbose,paramSolver,log_path,ampl_models_dir)
%Chercher rhohat_k' a-t-il un sens ? Non, si pour rho_k=0, les infections
%sont déjà éliminées
myTol = paramSolver.tolP0;
%addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';

%[ES,~,~,~,~,~,msgSolver] = ...
%    P1234_SICTPSEIIISSEIIS2_v6(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,f,'knitroampl',20,verbose);
[P,ES,msgSolver,ampl] = P_mod_v8(paramTab,paramRho,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,0);
ampl.close();

disp(P);

if ~isequal(msgSolver,0)
    disp('Le solver n a pas convergé')
    msgSol.mod = msgSolver;
end

N = length(mod);
for i=1:N
    dis = mod(i);
    Pi = P.(dis);
    if (Pi<myTol)
        disp(['L infection ', dis, ' est deja eliminee sans depistage volontaire']);
        msgSol.(dis) = '-2';
        infElim.(dis) = 1;
    else
        msgSol.(dis) = '0';
        infElim.(dis) = 0;
    end
end
end