function [U,dU] = U_SEIIS2SICRSEIIIS(beta1,gamma1,nu1,sigma1,eps1,...
                                    beta2,gamma2,nu2,sigma2,eps2,...
                                    betaI3,betaC3,gamma3,sigma3,theta3,...
                                    beta4,sigma4,tau4,gamma14,theta4,gamma34,nu4,...
                                    mu,b,vecRho,c,alpha1,alpha2,alpha3,alpha4,f)
    U  = zeros(1,length(vecRho));
    %R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*vecRho))./((mu + sigma1+vecRho).*(gamma1 + mu + nu1).*(mu + nu1 + vecRho)); 
    %R2 = (beta2*sigma2*(gamma2*(1-eps2) + mu + nu2 + eps2*vecRho))./((mu + sigma2+vecRho).*(gamma2 + mu + nu2).*(mu + nu2 + vecRho)); 
    %R3 = (betaI3*(theta3+vecRho+mu) + betaC3*sigma3)./((mu+theta3+vecRho).*(sigma3+gamma3+vecRho+mu));
    %R4 = sigma4*beta4*(tau4+theta4+vecRho+mu)./((theta4+vecRho+mu).*(gamma14+vecRho+tau4+mu).*(sigma4+vecRho+mu));

    i = 1;opt.TolFun=1e-5;
    for rho=vecRho
        %rho;
        if rho<min([alpha1,alpha2,alpha3,alpha4]) %R1(i)>=1 && R2(i)>=1 && R3(i)>=1 && R4(i)>=1
            U(i) = U1234_SEIIS2SICRSEIIIS(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,...
                                            betaI3,betaC3,gamma3,sigma3,theta3,...
                                            beta4,sigma4,tau4,gamma14,theta4,gamma34,nu4,mu,b,rho,c,...
                                            alpha1,alpha2,alpha3,alpha4,f,'fsolve',opt);
        %else
        elseif rho>=alpha1
            U(i) = U_SEIISSICRSEIIIS(beta2,gamma2,nu2,sigma2,eps2,betaI3,betaC3,gamma3,sigma3,theta3,...
                                 beta4,sigma4,tau4,nu4,gamma14,theta4,gamma34,mu,b,rho,c,alpha2,alpha3,alpha4,f);
        elseif rho>=alpha2%R2(i)<1
            U(i) = U_SEIISSICRSEIIIS(beta1,gamma1,nu1,sigma1,eps1,betaI3,betaC3,gamma3,sigma3,theta3,...
                                 beta4,sigma4,tau4,nu4,gamma14,theta4,gamma34,mu,b,rho,c,alpha1,alpha3,alpha4,f);
        elseif rho>=alpha3%R3(i)<1 
            U(i) = U_SEIIS2SEIIIS(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,...
                                 beta4,sigma4,tau4,nu4,gamma14,theta4,gamma34,mu,b,rho,c,alpha1,alpha2,alpha4,f);
        elseif rho>=alpha4%R4(i)<1
            U(i) = U_SEIIS2SICR(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,...
                                 betaI3,betaC3,gamma3,sigma3,theta3,mu,b,rho,c,alpha1,alpha2,alpha3,f);
        end
        %end
        i=i+1;
%         elseif R1(i)<=1 && R2(i)>1 && R3(i)<1 && R4(i)>1
%             U(i) = U12_SEIISSEIIIS(beta2,gamma2,nu2,sigma2,eps2,beta3,sigma3,tau3,theta3,gamma13,gamma33,nu3,mu,b,rho,c);
%         elseif R1(i)>1 && R2(i)<=1 && R3(i)<1 && R4(i)>1
%             U(i) = U12_SEIISSEIIIS(beta1,gamma1,nu1,sigma1,eps1,beta3,sigma3,tau3,theta3,gamma13,gamma33,nu3,mu,b,rho,c);
%         elseif R1(i)>1 && R2(i)>1 && R3(i)<=1 && R4(i)<=1
%             U(i) = U12_SEIIS2_v2(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,rho,c);
%         elseif R1(i)>1 && R2(i)<=1 && R3(i)<=1 && R4(i)<=1
%             U(i) = U_SEIISv2(eps1, beta1, sigma1, gamma1, mu, nu1, rho, c);
%         elseif R1(i)<=1 && R2(i)>1 && R3(i)<=1 && R4(i)<=1
%             U(i) = U_SEIISv2(eps2, beta2, sigma2, gamma2, mu, nu2, rho, c);
%         elseif R1(i)<=1 && R2(i)<=1 && R3(i)<=1 && R4(i)>1
%             U(i) = U_SEIIIS(b, beta3, sigma3, gamma13, gamma33, tau3, theta3, nu3, mu, rho, c);
%         elseif R1(i)<=1 && R2(i)<=1 && R3(i)>1 && R4(i)<=1
        
   
    end
    P = 0;
    dU  = 0;
end