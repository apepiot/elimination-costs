function [U,P] = U_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c)
% utility function of the SIRXSIS model (combined testing)
    R1     = beta1./(gamma1+mu);
    R2     = beta2./(gamma2+mu);
    alpha1 = beta1/s1*(1-1/R1);
    alpha2 = beta2/s2*(1-1/R2);
    gamma1p = gamma1+s1*rho;
    gamma2p = gamma2+s2*rho;
    gamma12t = s1*s2*rho;
    gamma2t = gamma2p-s1*s2*rho;
    R1p     = beta1./(gamma1p+mu);
    R2p     = beta2./(gamma2p+mu);
    Lambda1 = mu*(R1p-1);
    SES12   = b/mu*(Lambda1 + gamma2p+mu)./(Lambda1+beta2)./R1p;
    I1ES12  = b/mu*Lambda1./(beta2+gamma1p-gamma12t).*(gamma2t./beta1 + mu/b.*SES12);
    
    P1 = max(mu/beta1.*(R1p-1),0);
    P2 = max(1-1./R2p,0);
    P12 = P2 + (I1ES12*mu/b);
    P = P12.*(rho<alpha1 & rho<alpha2) + P1.*(rho<alpha1 & rho>=alpha2) + P2.*(rho>=alpha1 & rho<alpha2);
    U  = rho.*(P - c);
    moinsU = -U; 
end
