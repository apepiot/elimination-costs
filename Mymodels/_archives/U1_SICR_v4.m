function [U,dU] = U1_SICR_v4(param,mu,b,vecRho,c,f)
    gamma0=param.gamma;theta0=param.theta;betaI=param.betaI;betaC=param.betaC;
    sigma = param.sigma;
    %gamma = gamma0+vecRho; 
    theta = theta0+vecRho;
    [Rp,~,~] = Rp_SICR_v4(betaI,betaC,theta0,sigma,gamma0,mu,b,vecRho);
    P = min(f*(mu*((theta+mu+sigma)./(betaI*(theta+mu) + betaC*sigma)).*(Rp-1)),1);
    U = vecRho.*(P-c);
    dU = zeros(length(vecRho),1);
  
    for i=1:length(vecRho)
      rho = vecRho(i);
      dU(i) = (f*mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1)*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0)) - rho*((f*mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)^2) - betaI/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) + (betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)^2*(gamma0 + mu + rho + sigma)))*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0)) - (f*mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1))/(betaC*sigma + betaI*(mu + rho + theta0)) + (betaI*f*mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1)*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0))^2) - c;
    end
    %without f : 
    %dU = (mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1)*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0)) - rho*((mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)^2) - betaI/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) + (betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)^2*(gamma0 + mu + rho + sigma)))*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0)) - (mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1))/(betaC*sigma + betaI*(mu + rho + theta0)) + (betaI*mu*((betaC*sigma + betaI*(mu + rho + theta0))/((mu + rho + theta0)*(gamma0 + mu + rho + sigma)) - 1)*(mu + rho + sigma + theta0))/(betaC*sigma + betaI*(mu + rho + theta0))^2) - c;
end