function [cres] = findAllCswitch(vecC,tabRhohat,N,vecAlpha, pres)
 %precision : 1% of vecC interval, to find cswitch
%this function finds c switch, i.e. cost values for which the utility
%switch from a curve to another 
% input : vecC, tab with rhohats, N the total number of diseases

%first to do : check if cswitch does not appear outside vecC

%alpha's sorted
[alphaSorted, orderEl] = sort(vecAlpha); %i.e. disease elimination order

vecRhohat = tabRhohat.rhohat;
%rhohat jumps from Uijk to Ujk to Uk

presABS = pres*abs(vecC(1)-vecC(end));
if(N==2) 
    %U = U12 or Uj such that alphaj>alphai
    rhohatj  = tabRhohat.one(:,orderEl(end));
    rhohatij = tabRhohat.nm(:); %argmaxUij (not Uixj)
    [cres] = findCswitch(vecRhohat,rhohatj,rhohatij,vecC,presABS); 
end


end

