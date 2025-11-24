function [S,I1,I2,I12,IR2,R1] = ES_SIRSIS_V4(b,beta1,beta2,s1,s2,rho,gamma1,gamma2,mu)
%Theoretical calculations of ES for the SIRxSIS model with combined testing
%strategy.

    %I12
    gamma1p = gamma1+s1.*rho;
    R1p = beta1/(gamma1p+mu);
    gamma2p = gamma2+s2.*rho;
    R2p = beta2/(gamma2p+mu);

    Lambda2 = beta2*(1-1/R2p); 
    gamma12t = s1*s2*rho; 
    gamma1t = gamma1p-gamma12t;
    souslaracine = Lambda2^2*R1p^2*beta2^2*gamma12t^2 - 2*Lambda2*R1p^2*beta1*beta2^2*gamma12t*mu - 2*Lambda2*R1p^2*beta1*beta2*gamma1t*gamma12t*mu - 2*Lambda2*R1p*beta1*beta2^3*gamma12t - 2*Lambda2*R1p*beta1*beta2^2*gamma1t*gamma12t + 2*Lambda2*R1p*beta1*beta2^2*gamma12t^2 + 2*Lambda2*R1p*beta1*beta2^2*gamma12t*mu - 2*Lambda2*R1p*beta1*beta2*gamma2p*gamma12t^2 + 2*Lambda2*R1p*beta1*beta2*gamma1t*gamma12t*mu - 2*Lambda2*R1p*beta1*beta2*gamma12t^2*mu + R1p^2*beta1^2*beta2^2*mu^2 + 2*R1p^2*beta1^2*beta2*gamma1t*mu^2 + R1p^2*beta1^2*gamma1t^2*mu^2 + 2*R1p*beta1^2*beta2^3*mu + 4*R1p*beta1^2*beta2^2*gamma1t*mu + 2*R1p*beta1^2*beta2^2*gamma12t*mu - 2*R1p*beta1^2*beta2^2*mu^2 - 2*R1p*beta1^2*beta2*gamma2p*gamma12t*mu + 2*R1p*beta1^2*beta2*gamma1t^2*mu + 2*R1p*beta1^2*beta2*gamma1t*gamma12t*mu - 4*R1p*beta1^2*beta2*gamma1t*mu^2 - 2*R1p*beta1^2*beta2*gamma12t*mu^2 - 2*R1p*beta1^2*gamma2p*gamma1t*gamma12t*mu - 2*R1p*beta1^2*gamma1t^2*mu^2 - 2*R1p*beta1^2*gamma1t*gamma12t*mu^2 + beta1^2*beta2^4 + 2*beta1^2*beta2^3*gamma1t - 2*beta1^2*beta2^3*gamma12t - 2*beta1^2*beta2^3*mu + 2*beta1^2*beta2^2*gamma2p*gamma12t + beta1^2*beta2^2*gamma1t^2 - 2*beta1^2*beta2^2*gamma1t*gamma12t - 4*beta1^2*beta2^2*gamma1t*mu + beta1^2*beta2^2*gamma12t^2 + beta1^2*beta2^2*mu^2 + 2*beta1^2*beta2*gamma2p*gamma1t*gamma12t - 2*beta1^2*beta2*gamma2p*gamma12t^2 + 2*beta1^2*beta2*gamma2p*gamma12t*mu - 2*beta1^2*beta2*gamma1t^2*mu + 2*beta1^2*beta2*gamma1t*mu^2 - 2*beta1^2*beta2*gamma12t^2*mu + 2*beta1^2*beta2*gamma12t*mu^2 + beta1^2*gamma2p^2*gamma12t^2 + 2*beta1^2*gamma2p*gamma1t*gamma12t*mu + 2*beta1^2*gamma2p*gamma12t^2*mu + beta1^2*gamma1t^2*mu^2 + 2*beta1^2*gamma1t*gamma12t*mu^2 + beta1^2*gamma12t^2*mu^2;
    I12 = -(b*beta1*beta2^2 - b*(souslaracine)^(1/2) + b*beta1*beta2*gamma1t - b*beta1*beta2*gamma12t + b*beta1*gamma2p*gamma12t - b*beta1*beta2*mu - b*beta1*gamma1t*mu + b*beta1*gamma12t*mu - Lambda2*R1p*b*beta2*gamma12t + 2*Lambda2*R1p*b*gamma12t*mu + R1p*b*beta1*beta2*mu + R1p*b*beta1*gamma1t*mu - 2*Lambda2*R1p^2*b*gamma12t*mu)/(2*mu*(R1p*beta1*beta2*gamma12t - Lambda2*R1p^2*gamma12t^2 + R1p*beta1*gamma1t*gamma12t));    
    
    %Lambda1
    Lambda1 = mu*(R1p-1)+mu/b*R1p*gamma12t*I12;
    
    %S
    S  = b/mu*(mu+gamma2p/R1p+gamma12t*mu/b*I12)/(Lambda1+beta2);   
    %I2
    I2 = b/(mu*R1p)-S;
    %I1 
    I1 = b/beta1*(R1p-1)+(R1p/beta1*gamma12t-1)*I12;
    %R1
    R1 = b/(mu*R2p) - I1 -S;
    %IR2 
    IR2= b/mu*(1-1/R2p-1/R1p)+S-I12;

end

