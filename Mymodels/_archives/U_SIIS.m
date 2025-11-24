function [U,dU] = U_SIIS(p, beta, sigma, gamma0, mu, nu, rho, c)
    gamma = gamma0+rho;
    R = beta.*(p*sigma + gamma+mu+nu)./((gamma+mu+nu).*(p*sigma+mu+(1-p).*(nu+rho)));
    P = 1-1./R;
    alpha = - nu - (beta - gamma0 - 2*mu + gamma0*p + mu*p - p*sigma +...
        (gamma0^2*p^2 - 2*gamma0^2*p + gamma0^2 + 2*gamma0*mu*p^2 - 2*gamma0*mu*p +...
        2*gamma0*p^2*sigma - 2*gamma0*p*sigma - 2*beta*gamma0*p + 2*beta*gamma0 +....
        mu^2*p^2 + 2*mu*p^2*sigma - 2*beta*mu*p + p^2*sigma^2 - 4*beta*p^2*sigma + 2*beta*p*sigma + beta^2)^(1/2))/(2*(p - 1));
    rho = min(rho, alpha);
    U = max(rho.*(P-c),0);
    dU = rho*(((mu + p*sigma - (nu + rho)*(p - 1))*(gamma0 + mu + nu + rho))/(beta*(gamma0 + mu + nu + rho + p*sigma)^2) - (mu + p*sigma - (nu + rho)*(p - 1))/(beta*(gamma0 + mu + nu + rho + p*sigma)) + ((p - 1)*(gamma0 + mu + nu + rho))/(beta*(gamma0 + mu + nu + rho + p*sigma))) - c - ((mu + p*sigma - (nu + rho)*(p - 1))*(gamma0 + mu + nu + rho))/(beta*(gamma0 + mu + nu + rho + p*sigma)) + 1;
end