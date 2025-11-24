function [dY] = ODE_V3(t,Y,b,beta1,beta2,s,s1,s2,rho,gamma1,gamma2,mu,type)
%we add a diagnosis rate between I12 and S, now the diagnosis rates of I12
%to I1 and I12 to I2 change in consequences
%Parameters : s is the common sensitivity of 1 and 2, and should be between
% s1+s2-1 and min(s1,s2) whith s1 and s2 > 0.5.
    
    S=Y(1);I1=Y(2);I2=Y(3);I12=Y(4);IR2=Y(5);R1=Y(6);
    
    if strcmp(type,'frequency')
        N = sum(Y);
    end
    if strcmp(type,'density')
        N = 1;
    end
    
    lambda1 = beta1*(I1 +  I12)/N;
    lambda2 = beta2*(IR2 + I2 + I12)/N;
    gamma1p = gamma1+s1*rho;
    gamma2p = gamma2+s2*rho;

    dS      = b - (lambda1 + lambda2)*S + gamma2p*I2 +s*rho*I12- mu*S; %S
    dI1     = lambda1*S - lambda2*I1 + gamma2p*I12 - s*rho*I12 - (gamma1p + mu)*I1; %I1
    dI2     = lambda2*S - lambda1*I2 - (gamma2p + mu)*I2; %I2
    dI12    = lambda2*I1 + lambda1*I2 - (gamma1p + gamma2p + mu)*I12 + s*rho*I12;%I12
    dIR2    = gamma1p*I12 + lambda2*R1  - gamma2p*IR2 - mu*IR2 - s*rho*I12; %IR2
    dR1     = gamma1p*I1 - lambda2*R1 + gamma2p*IR2 - mu*R1; % R1
    
    dY = [dS; dI1; dI2; dI12; dIR2; dR1];   
end

