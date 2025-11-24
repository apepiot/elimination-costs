function [U,P] = U_SIRSIS2(rho,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c)
    R1     = beta1./(gamma1+mu);
    R2     = beta2./(gamma2+mu);
    R3     = beta3./(gamma3+mu);
    alpha1 = beta1/s1*(1-1/R1);
    alpha2 = beta2/s2*(1-1/R2);
    alpha3 = beta3/s3*(1-1/R3);
    
    gamma1p=gamma1+s1*rho;
    gamma2p=gamma2+s2*rho;
    gamma3p=gamma3+s3*rho;

    R1p = beta1./(gamma1p+mu);
    R2p = beta2./(gamma2p+mu);
    R3p = beta3./(gamma3p+mu);

    P1 = mu/beta1*(R1p-1);
    P2 = 1-1./R2p;
    P3 = 1-1./R3p; 
    %note : I could have also computed P12 directly instead of P1 and P2
    %separatly.

    [U12,P12] = U_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    [U13,P13] = U_SIRSIS7(rho,beta1,beta3,gamma1,gamma3,s1,s3,b,mu,c);
    [U23,P23] = U_SIS2(rho,beta2,beta3,gamma2,gamma3,s2,s3,b,mu,c);
    
    [U123,P123] = U123_SIRSIS2(rho,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c);
    P = P123.*(rho<alpha1 & rho<alpha2 & rho<alpha3) +...
        P1.*(rho<alpha1 & rho>=alpha2 & rho>=alpha3)+... %SIR
        P2.*(rho>=alpha1 & rho<alpha2 & rho>=alpha3)+... %SIS
        P3.*(rho>=alpha1 & rho>=alpha2 & rho<alpha3)+... %SIS
        P12.*(rho<alpha1 & rho<alpha2 & rho>=alpha3)+... %SIRxSIS
        P23.*(rho>=alpha1 & rho<alpha2 & rho<alpha3)+... %SISxSIS
        P13.*(rho<alpha1 & rho>=alpha2 & rho<alpha3); %SIRxSIS*
    U = rho.*(P-c);
end