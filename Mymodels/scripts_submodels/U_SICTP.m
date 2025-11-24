function [U] = U_SICTP(param,mu,b,vecRho,c,f)
    theta0 = param.theta0; betaI = param.betaI;   betaC = param.betaC;
    sigma  = param.sigma; eta   = param.eta;
    zeta   = param.zeta;  p = param.p;
    gamma0 = param.gamma0;
    [~,lambda,~] = Rp_SICTP(betaI,betaC,theta0,gamma0,sigma,zeta,eta,p,mu,b,vecRho);

    Ith  = -(b*lambda*(p - 1))./((lambda + mu)*(mu + vecRho + sigma));
    %Ipth = -(b*lambda*p*(zeta - 1))./((lambda + mu - lambda*zeta)*(eta + mu + sigma));
    Cth  = -(b*lambda*sigma*(p - 1))./((lambda + mu)*(mu + vecRho + sigma)*(mu + vecRho + theta0));
    %Cpth = -(b*lambda*p*sigma*(zeta - 1))./((lambda + mu - lambda*zeta)*(eta + mu + sigma)*(eta + mu + theta0));
    
    %P  = min(f*mu/b*(Ith+Ipth+Cth+Cpth),1);
    P = f*(Ith+Cth)/(b/mu);
    U  = vecRho.*(P-c);
end