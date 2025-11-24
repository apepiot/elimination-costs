function [U12] = U12_SISSIS(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c)
    Nequ = b./mu;
    gamma2p = gamma2 + s2*rho;
	gamma1p = gamma1 + s1*rho;
    R1p = beta1/(gamma1p + mu);
	R2p = beta2/(gamma2p + mu);
    gamma1t = gamma1p-s1*s2*rho;
    gamma2t = gamma2p-s1*s2*rho;
    num     = gamma1t/R2p+gamma2t/R1p+mu+s1*s2*rho;
    denum   = beta1+beta2-mu-s1*s2*rho;
    S12     = Nequ*num/denum;
    
    P12 = max(1-S12/Nequ,0);

    U12 = rho.*(P12-c);
end

