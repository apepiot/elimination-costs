function [U,vecPrev, vecR10,vecR20] = utilityFunctionSIS2_V2(beta1,beta2,gamma1,gamma2,s, vecDelta, p,mu,e1,e2)
    % Utility function of testing for the SISxSIS model
    % case were e1 = beta1, e2 = beta2;
    
    if e1 ~= beta1
        warning('e1 and beta1 are different: the result is not garanted')
    end
    
    if e2 ~= beta2
        warning('e2 and beta2 are different: the result is not garanted')
    end
    
    l       = length(vecDelta);
    N       = p/mu;

    vecR10  = beta1./(gamma1 + s*vecDelta +mu);
    vecR20  = beta2./(gamma2 + s*vecDelta +mu);
    
    %% How to choose the equilibrium ?
    % On consid�re ici que les maladies ne s'influencent pas l'une et
    % l'autre, donc la persistence de 1 est uniquement d�termin�e par R1
    % (idem pour 2).
    
    %case of disease free equilibrium 
    S0 = N*ones(1,l);

    %case of single infection equilibrium I1
    S1 = N./vecR10;

    %case of single infection equilibrium I2
    S2 = N./vecR20;

    %case of the coinfection equilibrium
    S12 = N.*(gamma1./vecR20 + gamma2./vecR10 + mu + s.*vecDelta)./(beta1+beta2 - mu-s.*vecDelta);

    
    Sequ = (vecR10<=1 & vecR20<=1).*S0 + (vecR10>1 & vecR20<=1).*S1 + (vecR10<=1 & vecR20>1).*S2 + (vecR10>1 & vecR20>1).*S12;
    %Sequ = S12;
    
    %prevalence
    vecPrev = max(N - Sequ,0);
    
    %conditionPositivite = (vecPrev>0);
    %U = vecGamma1.*vecPrev./N;
    U = vecDelta.*vecPrev./N;
    %Gamma th�orique qui maximise la fonction

end

