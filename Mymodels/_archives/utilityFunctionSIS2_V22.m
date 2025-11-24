function [U,P, R1p,R2p,P1,P2,P12] = utilityFunctionSIS2_V22(beta1,beta2,gamma1,gamma2,s1,s2,vecRho,p,mu)
    % Utility function of testing for the SISxSIS model , s1 diff of s2 and
    % I12-> S
    
    l       = length(gamma1);
    N       = p/mu;
    gamma1p = gamma1 + s1*vecRho;
    gamma2p = gamma2 + s2*vecRho;
    R1p  = beta1./(gamma1p+mu);
    R2p  = beta2./(gamma2p+mu);
    
    %% How to choose the equilibrium ?
    % On consid�re ici que les maladies ne s'influencent pas l'une et
    % l'autre, donc la persistence de 1 est uniquement d�termin�e par R1
    % (idem pour 2).
   
    %case of disease free equilibrium 
    s0 = ones(1,l);

    %case of single infection equilibrium I1
    s1 = 1./R1p;
    P1 = 1-s1;

    %case of single infection equilibrium I2
    s2 = 1./R2p;
    P2 = 1-s2;

    %case of the coinfection equilibrium
    s12 = (gamma2p./R1p + gamma1p./R2p + mu + (s1.*s2.*vecRho).*(1-1./R1p-1./R2p))./(beta1+beta2-mu-s1.*s2.*vecRho);
    P12 = 1-s12;
    
    %sequ = (R1p<=1 & R2p<=1).*s0 + (R1p>1 & R2p<=1).*s1 + (R1p<=1 & R2p>1).*s2 + (R1p>1 & R2p>1).*s12;
    
    %prevalence
    P = (R1p<=1 & R2p<=1).*0 + (R1p>1 & R2p<=1).*P1 + (R1p<=1 & R2p>1).*P2 + (R1p>1 & R2p>1).*P12;
    
    %utility
    U = vecRho.*P;

end

