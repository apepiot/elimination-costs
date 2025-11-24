function [dY] = ODE_SEIIS_v4(t,Y,beta,nu,gamma0,sigma,p,rho,mu,b)
%ODE to write
    %S,E,I1,I2,I3
    S=Y(1);E=Y(2);IA=Y(3);IS=Y(4);
    N=sum(Y);
    lambda = beta*(IS+IA)./N;
    dS  = b - lambda*S - mu*S + rho*E +(nu+rho)*IA+(gamma0+nu)*IS;
    dE  = lambda*S - (sigma+rho+mu)*E;
    dIA = (1-p)*sigma*E - (nu+rho+mu)*IA;
    dIS = p*sigma*E - (gamma0+nu+mu)*IS;
    
    dY = [dS;dE;dIA;dIS];
end