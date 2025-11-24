function [U,dU] = U_SEIIS2SEIIIS(beta1,gamma1,nu1,sigma1,eps1,...
                                     beta2,gamma2,nu2,sigma2,eps2,...
                                     beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,...
                                     mu,b,vecRho,c,alpha1,alpha2,alpha3,f)
    U  = zeros(1,length(vecRho));opt.TolFun=1e-5;
    %R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*vecRho))./((mu + sigma1+vecRho).*(gamma1 + mu + nu1).*(mu + nu1 + vecRho)); 
    %R2 = (beta2*sigma2*(gamma2*(1-eps2) + mu + nu2 + eps2*vecRho))./((mu + sigma2+vecRho).*(gamma2 + mu + nu2).*(mu + nu2 + vecRho)); 
    %R3 = sigma3*beta3*(tau3+theta3+vecRho+mu)./((theta3+vecRho+mu).*(gamma13+vecRho+tau3+mu).*(sigma3+vecRho+mu));
    i = 1;
    for rho=vecRho
        if rho<min([alpha1,alpha2,alpha3])%R1(i)>1 && R2(i)>=1 && R3(i)>=1
            U(i) = U123_SEIIS2SEIIIS(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,...
                                     beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,mu,b,rho,c,alpha1,alpha2,alpha3,f,'fsolve',opt);
        elseif rho>=alpha1 %R1(i)<=1 && R2(i)>1 && R3(i)>1
            U(i) = U_SEIISSEIIIS(beta2,gamma2,nu2,sigma2,eps2,beta3,sigma3,tau3,theta3,gamma13,gamma33,nu3,mu,b,rho,c,alpha2,alpha3,f);
        elseif rho>=alpha2 %R1(i)>1 && R2(i)<=1 && R3(i)>1
            U(i) = U_SEIISSEIIIS(beta1,gamma1,nu1,sigma1,eps1,beta3,sigma3,tau3,theta3,gamma13,gamma33,nu3,mu,b,rho,c,alpha1,alpha3,f);
        elseif rho>=alpha3 %R1(i)>1 && R2(i)>1 && R3(i)<=1
            U(i) = U12_SEIIS2_v2(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,b,mu,rho,c,alpha1,alpha2,f,'fsolve',opt);
%         elseif R1(i)>1 && R2(i)<=1 && R3(i)<=1
%             U(i) = U_SEIISv2(eps1, beta1, sigma1, gamma1, mu, nu1, rho, c);
%         elseif R1(i)<=1 && R2(i)>1 && R3(i)<=1
%             U(i) = U_SEIISv2(eps2, beta2, sigma2, gamma2, mu, nu2, rho, c);
%         elseif R1(i)<=1 && R2(i)<=1 && R3(i)>1
%             U(i) = U_SEIIIS(b, beta3, sigma3, gamma13, gamma33, tau3, theta3, nu3, mu, rho, c); 
        end
        i=i+1;
    end
    P = 0;
    dU  = 0;
end