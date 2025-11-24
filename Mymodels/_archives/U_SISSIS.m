function [U] = U_SISSIS(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c)
    R10 = beta1/(gamma1+mu);
    R20=beta2/(gamma2+mu);
    alpha1 = beta1*(1-1/R10)/s1;
    alpha2 = beta2*(1-1/R20)/s2;
    
    R1p = beta1./(gamma1+s1*rho+mu);
    R2p = beta2./(gamma2+s2*rho+mu);
    P1 = 1-1./R1p; U1 = rho.*(P1-c);
    P2 = 1-1./R2p; U2 = rho.*(P2-c);
    
    U12 = U12_SISSIS(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    
    if alpha2<alpha1
        U = (rho<=alpha2).*U12 + (rho>alpha2).*U1 ;
    end
    
    if alpha1<alpha2 
        U = (rho<=alpha1).*U12 + (rho>alpha1).*U2 ;
    end
    
end

