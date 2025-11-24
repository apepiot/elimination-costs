function [U,dU] = U12_SEIIIS(b, beta, sigma, gamma10, gamma30, tau, theta, nu, mu, rho, c)  
    tspan = 1:100;
    odefun = ODE_SEIIS2_2(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                 beta2,nu2,eps2,sigma2,gamma2,mu,b,rho,...
                                 indexInf1,indexInf2)
    P = 1-1/Rp;
    U = max(rho.*(P-c),0)*(Rp>=1);
    dU = 1 - rho*(((mu + rho + sigma)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) + ((mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) + ((mu + rho + sigma)*(mu + rho + theta))/(beta*sigma*(mu + rho + tau + theta)) - ((mu + rho + sigma)*(mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)^2)) - ((mu + rho + sigma)*(mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) - c;
end