function [dU] = dU_SIRSIS(beta1,beta2,gamma1,gamma2,delta,s1, s2, mu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    gamma1p = gamma1 + s1*delta;
    Gamma1 = gamma1p + mu;
    gamma2p = gamma2 + s2*delta;
    R1p = beta1/Gamma1;
    R1 = beta1/(gamma1+mu);
    %R2p = beta2/(gamma2p+mu);
    R2 = beta2/(gamma2+mu);
    
    dC = (-s1/beta1*R1p^2*(beta2+gamma1p)-(R1p-1)*s1)/(beta2+gamma1p)^2;
    dA = mu*(R1p-1)/(beta2+gamma1p) + delta*mu*dC;
    B = gamma2p/beta1 + 1/R1p*(mu*(R1p-1)+gamma2p+mu)/(mu*(R1p-1)+beta2);
    C = (R1p-1)/(beta2+gamma1p);
    A = delta*mu*C;
    U = mu + gamma2p/R1p;
    dU = (gamma2*s1 + s2*(gamma1+mu) + 2*s1*s2*delta)/beta1;
    V = mu*(R1p-1) + beta2;
    dV = -mu*s1/beta1*R1p^2;
    dB = s2/beta1 + (dU*V - U*dV)/V^2; 
    %dB = s2/beta1 + (((gamma2*s1+2*s1*s2*delta)/beta1 + s2/R1)*(mu*R1p-mu+beta2) + (mu+gamma2p/R1p)*(mu*s1/beta1*R1p^2))/(mu*R1p-mu+beta2)^2;

    dUI1 = dA*B + dB*A;
    dUI2 = 1-1/R2-2*delta*s2/beta2;

    dU = dUI1 + dUI2; %est-ce que l'expression est bonne ?
    
end

