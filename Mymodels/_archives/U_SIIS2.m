function [U] = U_SIIS2(b,beta1,beta2,gamma10,gamma20,eps1,eps2,sigma1,sigma2,nu1,nu2,mu,vecRho,c)
 %utility function of the SIISxSIIS model %numerically
    alpha1 = - nu1 - (beta1 - gamma10 - 2*mu + gamma10*eps1 + mu*eps1 - eps1*sigma1 +...
        (gamma10^2*eps1^2 - 2*gamma10^2*eps1 + gamma10^2 + 2*gamma10*mu*eps1^2 - 2*gamma10*mu*eps1 +...
        2*gamma10*eps1^2*sigma1- 2*gamma10*eps1*sigma1 - 2*beta1*gamma10*eps1 + 2*beta1*gamma10 +...
        mu^2*eps1^2 + 2*mu*eps1^2*sigma1 - 2*beta1*mu*eps1 + eps1^2*sigma1^2 -...
        4*beta1*eps1^2*sigma1 + 2*beta1*eps1*sigma1 + beta1^2)^(1/2))/(2*(eps1 - 1));
    alpha2 = - nu2 - (beta2 - gamma20 - 2*mu + gamma20*eps2 + mu*eps2 - eps2*sigma2 + (gamma20^2*eps2^2 - 2*gamma20^2*eps2 + gamma20^2 + 2*gamma20*mu*eps2^2 - 2*gamma20*mu*eps2 +...
        2*gamma20*eps2^2*sigma2- 2*gamma20*eps2*sigma2 - 2*beta2*gamma20*eps2 + 2*beta2*gamma20 + mu^2*eps2^2 + 2*mu*eps2^2*sigma2 - 2*beta2*mu*eps2 + eps2^2*sigma2^2 -...
        4*beta2*eps2^2*sigma2 + 2*beta2*eps2*sigma2 + beta2^2)^(1/2))/(2*(eps2 - 1));

    tspan = 0:1:500; U=zeros(1,length(vecRho)); i=1;
    for rho=vecRho
        if(rho<alpha1 && rho<alpha2)
            Y0 = [99; 1; 1 ; 1 ; 1 ; 1 ; 1 ;1 ; 1];
            options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
            [~,Ys] = ode45(@(t,Y) SIIS2(t,Y,b,beta1,beta2,gamma10,gamma20,rho,eps1,eps2,sigma1,sigma2,nu1,nu2,mu,'frequency'),tspan,Y0, options);
            N = b/mu;
            P = sum(Ys(end,2:9))/N;
            U(i) = rho*(P-c);
        elseif (rho<alpha1 && rho>=alpha2)
           U(i) = U_SIIS(eps1, beta1, sigma1, gamma10, mu, nu1, rho, c);
        elseif (rho<alpha2 && rho>=alpha1)
           U(i) = U_SIIS(eps2, beta2, sigma2, gamma20, mu, nu2, rho, c);
        else
           U(i) = 0;
        end
        i=i+1;
    end  
end