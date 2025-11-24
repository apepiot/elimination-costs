function [U,outputP] = U1234_SICTPSEIIISSEIIS2_v5(param1,param2,param3,param4,mu,b,paramRho,c,f,kit,solvingMethod)
% this code computes the utility for some rho_S where S is a set of
% infections 
% U(rho_k,c) = rho_k*(Pi_k(rho_k)-c)
% param1 : HIV parameters
% param2 : syphilis parameters
% param3 : Ct parameter
% param4 : Ng parameters
% kit : is a cell containing the infections considered in the utility
% paramRho : structure containing all the rho, eta rates



% if (paramRho.rho_h ~= param1.rhob || paramRho.rho_s ~= param2.rhob || paramRho.rho_c ~= param3.rhob || paramRho.rho_g ~=param4.rhob)
%     disp('erreur dans les parametres rho')
%     return;
% end

% Theoretical total population at the equilibrium
N=b/mu; 

%Computation of ES
[ES,~,~,iterNo,restart,changeSolver] = P1234_SICTPSEIIISSEIIS2_v5(param1,param2,param3,param4,mu,b,paramRho,f,solvingMethod,25);

outputP.iterNo = iterNo;
outputP.restart = restart;
outputP.changeSolver = changeSolver;

[HIV,syphilis,Ct,Ng,tot,N] = assigningPrevalence(ES);

%On considère tous les infectés/asymptomatiques et on enlève au fur et à
%mesure ceux qui ne doivent pas être dans P_k

infected_and_asymptomatic_HIV    = HIV.asymptomaticTotal;

if ismember('HIV',kit) && ismember('syphilis',kit) && ismember('Ct',kit) && ismember('Ng',kit)
    infected_and_asymptomatic_others = setdiff(tot.asymptomatic,infected_and_asymptomatic_HIV);
    Ph = ES(HIV.asymptomaticTotal)/N;
    Pscg_h = ES(infected_and_asymptomatic_others)/N;
    P_hscg = max(f*Ph + Pscg_h,1);
    U = rho_hscg*(P_hscg - c);
end

if ismember('HIV',kit) && ismember('syphilis',kit) && ismember('Ct',kit) && ~ismember('Ng',kit)
    infected_and_asymptomatic_tot    = unique([HIV.asymptomaticTotal,syphilis.asymptomatic,Ct.asymptomatic]); %voir notes du 8/11
    infected_and_asymptomatic_others = setdiff(infected_and_asymptomatic_tot,infected_and_asymptomatic_HIV);
    Ph = ES(infected_and_asymptomatic_HIV)/N;
    Psc_h = S(infected_and_asymptomatic_others)/N;
    P_hscg = max(f*Ph + Psc_h,1);
    U = paramRho.rho_hs*(P_hs - c);
end

%faire les autres cas

end