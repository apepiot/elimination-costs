function [dY] = ODE_V2(t,Y,b,beta1,beta2,s1,s2,gamma1,gamma2,mu,rho1,rho2,type)
%Each disease has its own voluntary testing rate. Disease-driven strategy
%(strategy 1)
    
    S=Y(1);I1=Y(2);I2=Y(3);I12=Y(4);IR2=Y(5);R1=Y(6);
    
    if strcmp(type,'frequency')
        N = sum(Y);
    end
    if strcmp(type,'density')
        N = 1;
    end
    
    gamma1p = gamma1+s1*rho1;
    gamma2p = gamma2+s2*rho2;
    lambda1 = beta1.*(I1 + I12)/N;
    lambda2 = beta2.*(IR2 + I2 + I12)/N;

    dS = b - (lambda1 + lambda2)*S + gamma2p*I2 - mu*S; %S
    dI1 = lambda1*S - lambda2*I1 + gamma2p*I12 - (gamma1p + mu)*I1; %I1
    dI2 = lambda2*S - lambda1*I2 - (gamma2p + mu)*I2; %I2
    dI12 = lambda2*I1 + lambda1*I2 - (gamma1p + gamma2p + mu)*I12;%I12
    dIR2 = gamma1p*I12 + lambda2*R1  - gamma2p*IR2 - mu*IR2; %IR2
    dR1 = gamma1p*I1 - lambda2*R1 + gamma2p*IR2 - mu*R1; % R1
    
    dY = [dS; dI1; dI2; dI12; dIR2; dR1];   
end

