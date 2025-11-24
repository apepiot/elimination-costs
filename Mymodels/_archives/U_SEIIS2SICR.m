function [U,dU] = U_SEIIS2SICR(beta1,gamma1,nu1,sigma1,eps1,...
                                     beta2,gamma2,nu2,sigma2,eps2,...
                                     betaI3,betaC3,gamma3,sigma3,theta3,...
                                     mu,b,vecRho,c,alpha1,alpha2,alpha3,f)
    U  = zeros(1,length(vecRho));
%     R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*vecRho))./((mu + sigma1+vecRho).*(gamma1 + mu + nu1).*(mu + nu1 + vecRho)); 
%     R2 = (beta2*sigma2*(gamma2*(1-eps2) + mu + nu2 + eps2*vecRho))./((mu + sigma2+vecRho).*(gamma2 + mu + nu2).*(mu + nu2 + vecRho)); 
%     R3 = (betaI3*(theta3+mu+vecRho) + betaC3*sigma3)./((mu+theta3+vecRho).*(sigma3+gamma3+vecRho+mu));    
    i = 1;opt.TolFun=1e-5;
    for rho=vecRho
        if rho<min([alpha1,alpha2,alpha3])%R1(i)>1 && R2(i)>1 && R3(i)>1
            U(i) = U123_SEIIS2SICR(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,...
                                     betaI3,betaC3,gamma3,sigma3,theta3,mu,b,rho,c,alpha1,alpha2,alpha3,f,'fsolve',opt);
        elseif rho>=alpha1 %R1(i)<=1 %R1(i)<=1 && R2(i)>1 && R3(i)>1
            U(i) = U_SEIISSICR(beta2,gamma2,nu2,sigma2,eps2,betaI3,betaC3,gamma3,sigma3,theta3,mu,b,rho,c,alpha2,alpha3,f);
        elseif rho>=alpha2 %R2(i)<=1 %R1(i)>1 && R2(i)<=1 && R3(i)>1
            U(i) = U_SEIISSICR(beta1,gamma1,nu1,sigma1,eps1,betaI3,betaC3,gamma3,sigma3,theta3,mu,b,rho,c,alpha1,alpha3,f);
        elseif rho>=alpha3 %R3(i)<=1 %R1(i)>1 && R2(i)>1 && R3(i)<=1
            U(i) = U_SEIIS2_v2(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,rho,c,alpha1,alpha2,f);
%         elseif R1(i)>1 && R2(i)<=1 && R3(i)<=1
%             U(i) = U_SEIISv2(eps1, beta1, sigma1, gamma1, mu, nu1, rho, c);
%         elseif R1(i)<=1 && R2(i)>1 && R3(i)<=1
%             U(i) = U_SEIISv2(eps2, beta2, sigma2, gamma2, mu, nu2, rho, c);
%         elseif R1(i)<=1 && R2(i)<=1 && R3(i)>1
%             U(i) = U_SICR(betaI3, betaC3, sigma3, gamma3, theta3, mu, rho, c);
        end
        i=i+1;
    end
    P = 0;
    dU  = 0;
end