function [U,dU] = U1_SEIIIS_v4(param,mu,b,rho,c,f)  
    sigma=param.sigma;  tau=param.tau;  theta=param.theta;  
    gamma10=param.gamma1;gamma30=param.gamma3;nu = param.nu;
    beta=param.beta; %alpha=param.alpha;
    %Rp = sigma*beta*(tau+theta+rho+mu)/((theta+rho+mu)*(gamma10+rho+tau+mu)*(sigma+rho+mu));
    %P = 1-1/Rp;
    [Rp,~,~] = Rp_SEIIIS_v4(beta,sigma,tau,nu,gamma10,theta,gamma30,mu,b,rho);
    P = 1-1./Rp;
    U = rho.*(P-c);
    %dU = 1 - rho*(((mu + rho + sigma)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) + ((mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) + ((mu + rho + sigma)*(mu + rho + theta))/(beta*sigma*(mu + rho + tau + theta)) - ((mu + rho + sigma)*(mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)^2)) - ((mu + rho + sigma)*(mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) - c;
    
    vecRho = rho; dU=zeros(1,length(vecRho));i=1;
    for rho=vecRho
        dU(i) = 1 - rho*(((mu + rho + sigma)*(mu + rho + theta)*(gamma30 + mu + nu + rho))/(beta*sigma*(tau*theta + (gamma30 + mu + nu + rho)*(mu + rho + tau + theta))) + ((mu + rho + sigma)*(mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(tau*theta + (gamma30 + mu + nu + rho)*(mu + rho + tau + theta))) + ((mu + rho + sigma)*(gamma30 + mu + nu + rho)*(gamma10 + mu + rho + tau))/(beta*sigma*(tau*theta + (gamma30 + mu + nu + rho)*(mu + rho + tau + theta))) + ((mu + rho + theta)*(gamma30 + mu + nu + rho)*(gamma10 + mu + rho + tau))/(beta*sigma*(tau*theta + (gamma30 + mu + nu + rho)*(mu + rho + tau + theta))) - ((mu + rho + sigma)*(mu + rho + theta)*(gamma30 + mu + nu + rho)*(gamma10 + mu + rho + tau)*(gamma30 + 2*mu + nu + 2*rho + tau + theta))/(beta*sigma*(tau*theta + (gamma30 + mu + nu + rho)*(mu + rho + tau + theta))^2)) - ((mu + rho + sigma)*(mu + rho + theta)*(gamma30 + mu + nu + rho)*(gamma10 + mu + rho + tau))/(beta*sigma*(tau*theta + (gamma30 + mu + nu + rho)*(mu + rho + tau + theta))) - c;
        i=i+1;
    end
end