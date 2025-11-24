function [moinsU1] = U1_SIRSIS7(rho,beta1,gamma1,s1,mu,c)
% utility function of the SIRXSIS model (combined testing)
    gamma1p = gamma1+s1*rho;
    R1p     = beta1./(gamma1p+mu);
    P1 = max(mu/beta1.*(R1p-1),0);
    U1  = rho.*(P1 - c);
    moinsU1 = -U1; 
end
