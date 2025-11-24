function [U,dU] = U_SEIIIS(b, beta, sigma, gamma10, gamma30, tau, theta, nu, mu, rho, c,alpha,f)  
    Rp = sigma*beta*(tau+theta+rho+mu)/((theta+rho+mu)*(gamma10+rho+tau+mu)*(sigma+rho+mu));
    P = 1-1/Rp;
    U = rho.*(P-c).*(rho<alpha) + (-c*alpha)*(rho>=alpha);
    dU = 1 - rho*(((mu + rho + sigma)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) + ((mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) + ((mu + rho + sigma)*(mu + rho + theta))/(beta*sigma*(mu + rho + tau + theta)) - ((mu + rho + sigma)*(mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)^2)) - ((mu + rho + sigma)*(mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) - c;
end