function [U,P] = U_SIS3(rho,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c)
    gamma1p = gamma1+s1*rho;
    gamma2p = gamma2+s2*rho;
    gamma3p = gamma3+s3*rho;
    R1p = beta1./(gamma1p+mu);
    R2p = beta2./(gamma2p+mu);
    R3p = beta3./(gamma3p+mu);
    
    %one disease prevalence
    P1 = 1-1./R1p;
    P2 = 1-1./R2p;
    P3 = 1-1./R3p;

    %two diseases prevalence
    [U12,P12] = U_SIS2(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    [U13,P13] = U_SIS2(rho,beta1,beta3,gamma1,gamma3,s1,s3,b,mu,c);
    [U23,P23] = U_SIS2(rho,beta2,beta3,gamma2,gamma3,s2,s3,b,mu,c);

    %three disease prevalence
    [U123,P123] = U123_SIS3(rho,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c);
        
    %prevalence
    P = (R1p<=1 & R2p<=1).*0 +...
        (R1p>1 & R2p<=1 & R2p<=1).*P1 + (R1p<=1 & R2p>1 & R3p<=1).*P2 + (R1p<=1 & R2p<=1 & R3p>1).*P3 +...
        (R1p>1 & R2p>1 & R3p<=1).*P12 + (R1p>1 & R2p<=1 & R3p>1).*P13 + (R1p<=1 & R2p>1 & R3p>1).*P23 +...
        (R1p>1 & R2p>1 & R3p>1).*P123;
    
    %utility
    U = rho.*(P-c);
end

