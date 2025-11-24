function [U,dU] = U_SEIISSICR(beta1,gamma1,nu1,sigma1,eps1,betaI2,betaC2,gamma2,sigma2,theta2,mu,b,vecRho,c,alpha1,alpha2,f)
    U = zeros(1,length(vecRho));opt.TolFun=1e-5;
    %souslaracine = @(gamma0,beta,nu,p,sigma) (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2);
    %alpha1 = (beta1*eps1*sigma1 + sqrt(souslaracine(gamma1,beta1,nu1,eps1,sigma1)))/(2*(gamma1+mu+nu1)) - (2*mu+nu1+sigma1)/2;
    %alpha2 = (beta2*eps2*sigma2 + sqrt(souslaracine(gamma2,beta2,nu2,eps2,sigma2)))/(2*(gamma2+mu+nu2)) - (2*mu+nu2+sigma2)/2;
    %R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*vecRho))./((mu + sigma1+vecRho).*(gamma1 + mu + nu1).*(mu + nu1 + vecRho)); 
    %R2 = (betaI2.*(theta2+vecRho+mu)+betaC2*sigma2)./((sigma2+gamma2+vecRho+mu).*(theta2+vecRho+mu));
    i=1;
    for rho=vecRho
        if rho<alpha1 && rho<alpha2 %R1(i)>=1 && R2(i)>=1
            U(i) = U12_SEIISSICR(beta1,gamma1,nu1,sigma1,eps1,betaI2,betaC2,gamma2,sigma2,theta2,mu,b,rho,c,alpha1,alpha2,f,'fsolve',opt);
        elseif rho>=alpha1 && rho<alpha2 %R1(i)<1 && R2(i)>=1
            U(i) = U_SICR(betaI2, betaC2, sigma2, gamma2, theta2, mu, rho, c ,alpha2,f);
        elseif rho>=alpha2 && rho<alpha1 %R1(i)>=1 && R2(i)<1
            U(i) = U_SEIISv2(eps1, beta1, sigma1, gamma1, mu, nu1, rho, c, alpha1, f);
        elseif rho>=alpha2 && rho>=alpha1
            U(i) = -c*max(alpha1,alpha2); %instead of -c*rho ?
        end
        i=i+1;
    end
    %P=0;
    dU  = 0;
end