function [moinsU2] = U2_SIRSIS7(rho,beta2,gamma2,s2,mu,c)
% utility function of the SIRXSIS model (combined testing)
    gamma2p = gamma2+s2*rho;
    R2p     = beta2./(gamma2p+mu);
    P2 = max(1-1./R2p,0);
    U2  = rho.*(P2 - c);
    moinsU2 = -U2; 
end
