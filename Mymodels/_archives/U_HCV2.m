function [moinsU] = U_HCV2(rho, beta, gamma, s, sigma, omega, mu, c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    gammap = gamma+s*rho;
    Rp = beta.*(omega*sigma+mu+gammap)./((sigma+gammap+mu).*(mu+gammap));

    %prevalence 
    P = mu/beta.*(Rp-1);
    
    %utility
    U = rho.*(P-c);
    moinsU = -U;
end

