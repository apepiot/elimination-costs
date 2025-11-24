function [I] = i_WA(t, p, mu, beta, gamma, N0, I0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    c1 = exp((N0 - p/mu)*beta/mu)/I0;
    c2 = p/mu - (gamma+mu)/beta;
    %c1=0;
    %c1=1000;
    param = [p, mu, beta, gamma, N0, I0];
    
    Z1 = 1./I_WA_sousfct(t, param);
    
    INT = integral(@(u) I_WA_sousfct(u,param) ,0,t);
    Z2 = (c1 + beta* INT);
    
    I = 1./(Z1*Z2) +c2;

    
end

