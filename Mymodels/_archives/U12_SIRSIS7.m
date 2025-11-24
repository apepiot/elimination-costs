function [U12,P12,dU12] = U12_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c)
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
    P2 = 1-1./R2p;
    P12 = P2 + (I1ES12*mu/b);
    U12  = rho.*(P12 - c);
    %moinsU12 = -U12; 
    dU12 = (mu*((gamma2 + rho*s2 - rho*s1*s2)/beta1 + ((gamma1 + mu + rho*s1)*(gamma2 + mu + rho*s2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1))))*(beta1/(gamma1 + mu + rho*s1) - 1))/(beta2 + gamma1 + rho*s1 - rho*s1*s2) - (gamma2 + mu + rho*s2)/beta2 - rho*(s2/beta2 - (mu*(beta1/(gamma1 + mu + rho*s1) - 1)*((s2 - s1*s2)/beta1 + (s1*(gamma2 + mu + rho*s2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1))) + ((s2 - (beta1*mu*s1)/(gamma1 + mu + rho*s1)^2)*(gamma1 + mu + rho*s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1))) + (mu*s1*(gamma2 + mu + rho*s2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))/((beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1))^2*(gamma1 + mu + rho*s1))))/(beta2 + gamma1 + rho*s1 - rho*s1*s2) + (mu*((gamma2 + rho*s2 - rho*s1*s2)/beta1 + ((gamma1 + mu + rho*s1)*(gamma2 + mu + rho*s2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1))))*(s1 - s1*s2)*(beta1/(gamma1 + mu + rho*s1) - 1))/(beta2 + gamma1 + rho*s1 - rho*s1*s2)^2 + (beta1*mu*s1*((gamma2 + rho*s2 - rho*s1*s2)/beta1 + ((gamma1 + mu + rho*s1)*(gamma2 + mu + rho*s2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))))/((gamma1 + mu + rho*s1)^2*(beta2 + gamma1 + rho*s1 - rho*s1*s2))) - c + 1;
end
