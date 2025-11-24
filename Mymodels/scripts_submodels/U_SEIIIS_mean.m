function [U] = U_SEIIIS_mean(rho,c)
    mu          = 1/30.6; b=10000;
    P_VT        = 6.6/100;                            
    R_VT        = 1/(1-P_VT); 
    rhob        = 1/2.48;      
    sigma       = 365/25; 
    tau         = 365/46;
    theta       = 12/3.6; 
    gamma10     = 0;
    gamma1      = gamma10+rhob;
    gamma30     = 1/20;
    gamma3      = gamma30+rhob;
    nu          = 0; 
    beta        = R_VT/((sigma.*((gamma3+mu+nu).*(mu+rhob+tau+theta)+tau*theta))./((mu+rhob+sigma).*(mu+rhob+theta).*(gamma3+mu+nu).*(gamma1+mu+tau)));

    S = (b*mu^4 + b*rho^4 + 4*b*mu*rho^3 + 4*b*mu^3*rho + b*nu*rho^3 + b*mu^3*sigma + b*mu^3*tau + b*mu^3*theta + b*rho^3*sigma + b*rho^3*tau + b*rho^3*theta + 6*b*mu^2*rho^2 + b*gamma10*mu^3 + b*gamma30*mu^3 + b*gamma10*rho^3 + b*gamma30*rho^3 + b*mu^3*nu + b*gamma10*gamma30*mu^2 + b*gamma10*gamma30*rho^2 + b*gamma10*mu^2*nu + 3*b*gamma10*mu*rho^2 + 3*b*gamma10*mu^2*rho + 3*b*gamma30*mu*rho^2 + 3*b*gamma30*mu^2*rho + b*gamma10*nu*rho^2 + b*gamma10*mu^2*sigma + b*gamma30*mu^2*sigma + b*gamma30*mu^2*tau + b*gamma10*mu^2*theta + b*gamma30*mu^2*theta + b*gamma10*rho^2*sigma + b*gamma30*rho^2*sigma + b*gamma30*rho^2*tau + b*gamma10*rho^2*theta + b*gamma30*rho^2*theta + 3*b*mu*nu*rho^2 + 3*b*mu^2*nu*rho + b*mu^2*nu*sigma + b*mu^2*nu*tau + b*mu^2*nu*theta + 3*b*mu*rho^2*sigma + 3*b*mu^2*rho*sigma + 3*b*mu*rho^2*tau + 3*b*mu^2*rho*tau + 3*b*mu*rho^2*theta + 3*b*mu^2*rho*theta + b*nu*rho^2*sigma + b*nu*rho^2*tau + b*mu^2*sigma*tau + b*nu*rho^2*theta + b*mu^2*sigma*theta + b*mu^2*tau*theta + b*rho^2*sigma*tau + b*rho^2*sigma*theta + b*rho^2*tau*theta + 2*b*gamma10*gamma30*mu*rho + b*gamma10*gamma30*mu*sigma + b*gamma10*gamma30*mu*theta + b*gamma10*gamma30*rho*sigma + b*gamma10*gamma30*rho*theta + 2*b*gamma10*mu*nu*rho + b*gamma10*gamma30*sigma*theta + b*gamma10*mu*nu*sigma + b*gamma10*mu*nu*theta + 2*b*gamma10*mu*rho*sigma + 2*b*gamma30*mu*rho*sigma + 2*b*gamma30*mu*rho*tau + 2*b*gamma10*mu*rho*theta + 2*b*gamma30*mu*rho*theta + b*gamma10*nu*rho*sigma + b*gamma30*mu*sigma*tau + b*gamma10*nu*rho*theta + b*gamma10*mu*sigma*theta + b*gamma30*mu*sigma*theta + b*gamma30*mu*tau*theta + b*gamma10*nu*sigma*theta + b*gamma30*rho*sigma*tau + b*gamma10*rho*sigma*theta + b*gamma30*rho*sigma*theta + 2*b*mu*nu*rho*sigma + b*gamma30*rho*tau*theta + 2*b*mu*nu*rho*tau + 2*b*mu*nu*rho*theta + b*gamma30*sigma*tau*theta + b*mu*nu*sigma*tau + b*mu*nu*sigma*theta + b*mu*nu*tau*theta + 2*b*mu*rho*sigma*tau + 2*b*mu*rho*sigma*theta + 2*b*mu*rho*tau*theta + b*nu*rho*sigma*tau + b*nu*rho*sigma*theta + b*nu*rho*tau*theta + b*mu*sigma*tau*theta + b*nu*sigma*tau*theta + b*rho*sigma*tau*theta)/(mu*sigma*(beta*mu^2 + beta*rho^2 + beta*gamma30*mu + beta*gamma30*rho + beta*gamma30*tau + beta*gamma30*theta + beta*mu*nu + 2*beta*mu*rho + beta*nu*rho + beta*mu*tau + beta*mu*theta + beta*nu*tau + beta*nu*theta + beta*rho*tau + beta*rho*theta + beta*tau*theta));
    P = 1-S/(b/mu);
    U = rho*(P-c);
end

