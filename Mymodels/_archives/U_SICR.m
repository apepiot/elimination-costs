function [U,dU] = U_SICR(betaI, betaC, sigma, gamma0, theta0, mu, rho, c, alpha, f)
    gamma = gamma0+rho; theta = theta0+rho;
    Rp = (betaI*(theta+mu) + betaC*sigma)/((mu+theta)*(sigma+gamma+mu));
    P = f*max(mu*((theta+mu+sigma)./(betaI*(theta+mu) + betaC*sigma)).*(Rp-1),0);

    %alpha = (betaI-gamma0-2*mu-sigma-theta0)/2+((betaI-gamma0-sigma+theta0)^2+4*betaC*sigma)^(1/2)/2;
    %rho = min(rho, alpha);
    %U = max(rho.*(P-c),0);
    %f = 10;
    %disp(['attention: there is a factor ',num2str(f), ' before prevalence of HIV']);
    U = rho.*(P-c)*(rho<alpha)+ (-alpha*c).*(rho>=alpha);
    dU = (f*mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1)*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0)) - rho*((f*mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)^2) - betaI/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) + (betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)^2*(gamma0 + mu + rho + sigma)))*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0)) - (f*mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1))/(betaC*sigma + betaI*(mu + rho + theta0)) + (betaI*f*mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1)*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0))^2) - c;
    %without f : 
    %dU = (mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1)*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0)) - rho*((mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)^2) - betaI/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) + (betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)^2*(gamma0 + mu + rho + sigma)))*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0)) - (mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1))/(betaC*sigma + betaI*(mu + rho + theta0)) + (betaI*mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1)*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0))^2) - c;
end