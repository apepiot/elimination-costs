function [Y] = I_WA_sousfct(t, param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    p       = param(1);
    mu      = param(2);
    beta    = param(3);
    gamma   = param(4);
    N0      = param(5);
    I0      = param(6);
    alpha = N0 - p/mu;
    Y = exp((alpha*beta.*exp(-mu*t)/mu) - (p/mu*beta - gamma - mu)*t);
end

