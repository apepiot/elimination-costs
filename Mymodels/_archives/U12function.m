function [U12,s] = U12function(beta1,beta2,gamma1,gamma2,s1,s2,mu,rho1,rho2)
    gamma1p = gamma1+s1.*rho1;
    gamma2p = gamma2+s2.*rho2;
    R1p = beta1./(gamma1p+mu);
    R2p = beta2./(gamma2p+mu);
    s = (gamma1p./R2p + gamma2p./R1p + mu)./(beta1+beta2-mu); %s/N
    P12 = 1 - s;
    U12 = (rho1+rho2).*P12;
    
end

