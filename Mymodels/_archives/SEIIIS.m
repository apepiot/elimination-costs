function dY = SEIIIS(t,Y,b,beta, sigma, gamma10, gamma30, tau, theta, nu, mu,rho)
    % SEI1I2I3 model, syphilis
    S=Y(1);E=Y(2);I1=Y(3);I2=Y(4);I3=Y(5);

    N = b/mu;
    lambda = beta*(I1+I2)/N;

    dS  = b - lambda*S - mu*S + rho*E +(gamma10+rho)*I1+rho*I2+(nu+gamma30+rho)*I3;
    dE  = lambda*S - (sigma+rho+mu)*E;
    dI1 = sigma*E -(gamma10+rho+tau+mu)*I1;
    dI2 = tau*I1 - (theta+rho+mu)*I2;
    dI3 = theta*I2 - (nu+gamma30+rho+mu)*I3;
    
    dY = [dS; dE; dI1; dI2; dI3];   
end

