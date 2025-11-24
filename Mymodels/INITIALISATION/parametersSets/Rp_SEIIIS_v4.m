function [Rp,Lambdap,alphaS] = Rp_SEIIIS_v4(beta,sigma,tau,nu,gamma10,theta,gamma30,mu,b,rho)
    %%Case lambda = beta*(I1+I2+I3)/N
    Rp = (beta*sigma).*((gamma30+mu+nu+rho).*(mu+rho+tau+theta)+tau*theta)./((mu + rho + sigma).*(mu + rho + theta).*(gamma30 + mu + nu + rho).*(gamma10 + mu + rho + tau));
    
    Lambdap = max(beta*(Rp-1)./(Rp+beta./(sigma+rho+mu)),0);
    Rpfun  = @(rho) (((beta*sigma)*((gamma30+mu+nu+rho)*(mu+rho+tau+theta)+tau*theta)./...
        ((mu + rho + sigma)*(mu + rho + theta)*(gamma30 + mu + nu + rho)*(gamma10 + mu + rho + tau)))-1);
    alphaS=-1;
    %if (rho==0)
    iterMax=100;startingPoint=0;iter=0;
    while (alphaS<0 && iter<iterMax)
        alphaS = fzero(Rpfun, startingPoint); 
        iter=iter+1; startingPoint = startingPoint+1;
        %iter;
    end
    if (alphaS<0)
        error('alphaS<0')
    end
    %end
end

