function [dY] = ODE_SEIIIS_v4(t,Y,beta,sigma, tau, theta, gamma10,gamma30,nu,rho,mu,b)
    %ODE to write
    %S,E,I1,I2,I3
    S=Y(1);E=Y(2);I1=Y(3);I2=Y(4);I3=Y(5);
    N=sum(Y);
    lambda = beta*(I1+I2+I3)./N;
    dS  = b - lambda*S - mu*S + rho*E +(gamma10+rho)*I1+rho*I2+(nu+gamma30+rho)*I3;
    dE  = lambda*S - (sigma+rho+mu)*E;
    dI1 = sigma*E - (gamma10+rho+tau+mu)*I1;
    dI2 = tau*I1 - (theta+rho+mu)*I2;
    dI3 = theta*I2 - (nu+gamma30+rho+mu)*I3;

    dY = [dS;dE;dI1;dI2;dI3];
end