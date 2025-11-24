function [U,dU] = U_SIJSv2(p, beta, sigma, gamma0, mu, nu, rho, c)
    gamma = gamma0+rho;
    R = beta.*(gamma+nu+mu+p*sigma+rho)./((gamma+nu+mu+rho).*(p*sigma+(1-p)*nu+rho+mu));
    P = 1-1./R;
    alpha = (beta-(gamma0+2*mu+(1-p)*nu+nu+p*sigma))/2 + sqrt((gamma0+nu*p+beta+p*sigma)^2-4*p*sigma*(gamma0+nu*p))/2;
    rho = min(rho, alpha);
    U = max(rho.*(P-c),0);
    dU = 1 - rho*((gamma0 + mu + nu + rho)/(beta*(gamma0 + mu + nu + rho + p*sigma)) + (mu + rho + p*sigma - nu*(p - 1))/(beta*(gamma0 + mu + nu + rho + p*sigma)) - ((gamma0 + mu + nu + rho)*(mu + rho + p*sigma - nu*(p - 1)))/(beta*(gamma0 + mu + nu + rho + p*sigma)^2)) - ((gamma0 + mu + nu + rho)*(mu + rho + p*sigma - nu*(p - 1)))/(beta*(gamma0 + mu + nu + rho + p*sigma)) - c;
end