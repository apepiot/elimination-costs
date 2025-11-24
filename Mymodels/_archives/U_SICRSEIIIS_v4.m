function [U,dU] = U_SICRSEIIIS_v4(param1,param2,mu,b,vecRho,c,f)
    alpha1=param1.alpha;alpha2=param2.alpha;
    U = zeros(1,length(vecRho));opt.TolFun=1e-5;
    %R1 = (betaI1.*(theta1+vecRho+mu)+betaC1*sigma1)./((sigma1+gamma1+vecRho+mu).*(theta1+vecRho+mu)); %SICR
    %R2 = sigma2*beta2*(tau2+theta2+vecRho+mu)./((theta2+vecRho+mu).*(gamma12+vecRho+tau2+mu).*(sigma2+vecRho+mu)); %SEIIIS
    i=1;
    for rho=vecRho       
        if rho<alpha1 && rho<alpha2 %R1(i)>1 && R2(i)>1
            U(i) = U12_SICRSEIIIS_v4(param1,param2,mu,b,rho,c,f,'fsolve',opt);
        elseif rho>=alpha1 && rho<alpha2 %R1(i)<=1 && R2(i)>1
            U(i) = U_SEIIIS_v4(param2,mu,b,rho,c,f);
        elseif rho>=alpha2 && rho<alpha1 %R1(i)>1 && R2(i)<=1
            U(i) = U_SICR_v4(param1,mu,b,rho,c,f);
        else
            U(i) = -c*rho;
        end
        i=i+1;
    end
    %P = 0;
    dU  = 0;
end