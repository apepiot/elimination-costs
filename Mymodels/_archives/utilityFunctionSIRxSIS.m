function [U,vecPrev, vecR10,vecR20,P1,P2,P12, rho12max] = utilityFunctionSIRxSIS(beta1,beta2,gamma1,gamma2,vecDelta1, vecDelta2, s1, s2, p,mu)
    % Utility function of testing for the SISxSIS model
    % case were e1 = beta1, e2 = beta2; version2
    % U12 = rho*P12
    
    l       = length(vecDelta1);
    N       = p/mu;
    vecGamma1 = gamma1 + s1*vecDelta1;
    vecGamma2 = gamma2 + s2*vecDelta2;
    
    vecR10  = beta1./(vecGamma1 + mu);
    vecR20  = beta2./(vecGamma2 + mu);
    

    
    %% How to choose the equilibrium ?
    % On consid�re ici que les maladies ne s'influencent pas l'une et
    % l'autre, donc la persistence de 1 est uniquement d�termin�e par R1
    % (idem pour 2).
    
    %case of disease free equilibrium 
    %S0 = N*ones(1,l);
    P0 = zeros(1,l);
    
    %case of single infection equilibrium I1
    %S1 = N./vecR10.*ones(1,l);
    P1 = mu/beta1.*(vecR10-1); %divise par N

    %case of single infection equilibrium I2
    %S2 = N./vecR20.*ones(1,l);
    P2 = 1 - 1./vecR20; %divise par N

    %case of the coinfection equilibrium
    %S12 = N.*(gamma1./vecR20 + gamma2./vecR10 + mu + s.*vecDelta)./(beta1+beta2 - mu-s.*vecDelta);
    P12 = mu*(vecR10 -1)./(beta2+vecGamma1).*(vecGamma2./beta1 + 1./vecR10.*(mu.*(vecR10-1)+vecGamma2+mu)./(mu.*(vecR10-1)+beta2)) + (1-1./vecR20);
    [U12max,i12max] = max(vecDelta1.*P12);
    rho12max = vecDelta1(i12max);
    
    P12 = P12.*(vecR10>1 & vecR20>1);
    
    %S12  a revoir
    %Sequ = (vecR10<=1 & vecR20<=1).*S0 + (vecR10>1 & vecR20<=1).*S1 + (vecR10<=1 & vecR20>1).*S2 + (vecR10>1 & vecR20>1).*S12; %a revoir
    
    vecPrev = (vecR10<=1 & vecR20<=1).*P0 +...
        (vecR10>1 & vecR20<=1).*P1 + ...
        (vecR10<=1 & vecR20>1).*P2 + ...
        (vecR10>1 & vecR20>1).*P12; 
    %vecPrev = max(P12;
    
    
    %conditionPositivite = (vecPrev>0);
    U = vecDelta1.*vecPrev;
    
end

