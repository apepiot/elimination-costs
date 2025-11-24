function [U,P,dU] = U_SIAS(rho,beta,gamma0,ksi0,omega,s,b,mu,c)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    gamma=gamma0+s*rho;
    ksi=ksi0+s*rho;
    R = beta.*(gamma.*(1-omega)+ksi*omega+mu)./((gamma+mu).*(ksi+mu));
    P = 1-1./R;

    U = max(rho.*(P-c),0);
    dU  = 1 - rho*((s*(gamma0 + mu + rho*s))/(beta*(mu - (gamma0 + rho*s)*(omega - 1) + omega*(ksi0 + rho*s))) + (s*(ksi0 + mu + rho*s))/(beta*(mu - (gamma0 + rho*s)*(omega - 1) + omega*(ksi0 + rho*s))) - ((omega*s - s*(omega - 1))*(gamma0 + mu + rho*s)*(ksi0 + mu + rho*s))/(beta*(mu - (gamma0 + rho*s)*(omega - 1) + omega*(ksi0 + rho*s))^2)) - ((gamma0 + mu + rho*s)*(ksi0 + mu + rho*s))/(beta*(mu - (gamma0 + rho*s)*(omega - 1) + omega*(ksi0 + rho*s))) - c;

end

