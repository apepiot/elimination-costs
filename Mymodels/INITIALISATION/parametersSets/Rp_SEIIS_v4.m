function [Rp,Lambdap,alpha,Pas] = Rp_SEIIS_v4(beta,nu,p,sigma,gamma0,mu,b,rho)
%%Case lambda = beta*(IS+IA)/N
Rp = (beta*sigma*(gamma0*(1-p) + mu + nu + p*rho))./((mu + sigma + rho).*(gamma0 + mu + nu).*(mu + nu + rho));   
Lambdap = ((Rp-1).*beta.*(sigma+rho+mu)./(beta+(sigma+rho+mu).*Rp));
alpha1 = -(2*gamma0*mu + gamma0*nu + gamma0*sigma + 3*mu*nu + mu*sigma + nu*sigma + (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2)^(1/2) + 2*mu^2 + nu^2 - beta*p*sigma)/(2*(gamma0 + mu + nu));
alpha2 = -(2*gamma0*mu + gamma0*nu + gamma0*sigma + 3*mu*nu + mu*sigma + nu*sigma - (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2)^(1/2) + 2*mu^2 + nu^2 - beta*p*sigma)/(2*(gamma0 + mu + nu));
alpha  = max(alpha1,alpha2); %a ecrire plus proprement, alpha2

S  = b/(mu*Rp);
E  = Lambdap*S./(sigma+rho+mu);
IA = (1-p)*sigma*E./(nu+rho+mu);
IS = p*sigma*E./(gamma0+nu+mu);
N  = b./mu;
Pas = (E+IA)./N;
end

