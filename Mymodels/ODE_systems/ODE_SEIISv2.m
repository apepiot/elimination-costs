function [dY] = ODE_SEIISv2(t,Y,b,beta,nu,p,sigma,gamma0,rho,mu)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    dY=zeros(4,1);
    lambda = beta*(Y(3)+Y(4))/sum(Y);
    dY(1)  = b - lambda*Y(1) + rho*Y(2) + (gamma0+nu)*Y(4) + (rho+nu)*Y(3) - mu*Y(1);
    dY(2)  = lambda*Y(1) - (sigma+rho)*Y(2) - mu*Y(2);
    dY(3)  = (1-p)*sigma*Y(2) - (rho+nu)*Y(3) - mu*Y(3); %3 : asymptomatic IA
    dY(4)  = p*sigma*Y(2) - (gamma0+nu)*Y(4) - mu*Y(4);  %4 : symptomatic IS
end