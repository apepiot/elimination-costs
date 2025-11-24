function [U,dU] = U_SEIISSEIIIS(beta1,gamma1,nu1,sigma1,eps1,beta2,sigma2,tau2,theta2,gamma12,gamma32,nu2,mu,b,vecRho,c,alpha1,alpha2,f)
    U = zeros(1,length(vecRho));opt.TolFun=1e-5;
    %R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*vecRho))./((mu + sigma1+vecRho).*(gamma1 + mu + nu1).*(mu + nu1 + vecRho)); 
    %R2 = sigma2*beta2*(tau2+theta2+vecRho+mu)./((theta2+vecRho+mu).*(gamma12+vecRho+tau2+mu).*(sigma2+vecRho+mu));
    i  = 1;
    for rho=vecRho
        if rho<alpha1 && rho<alpha2 %R1(i)>=1 && R2(i)>=1
            U(i) = U12_SEIISSEIIIS(beta1,gamma1,nu1,sigma1,eps1,beta2,sigma2,tau2,theta2,gamma12,gamma32,nu2,mu,b,rho,c,alpha1,alpha2,f,'fsolve',opt);
        elseif rho>=alpha1 && rho<alpha2 %R1(i)<1 && R2(i)>=1
            U(i) = U_SEIIIS(b, beta2, sigma2, gamma12, gamma32, tau2, theta2, nu2, mu, rho, c, alpha2,f)  ;
        elseif rho>=alpha2 && rho<alpha1 %R1(i)>=1 && R2(i)<1
            U(i) = U_SEIISv2(eps1, beta1, sigma1, gamma1, mu, nu1, rho, c, alpha1 ,f);
        else
            U(i) = -c*rho;
        end
        i=i+1;
    end
    %P   = 0;
    dU  = 0;
end