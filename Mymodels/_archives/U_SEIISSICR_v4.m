function [U,dU] = U_SEIISSICR_v4(param1,param2,mu,b,vecRho,c,f)
    U = zeros(1,length(vecRho));opt.TolFun=1e-5;
    alpha1=param1.alpha;alpha2=param2.alpha;
    i=1;
    for rho=vecRho
        if rho<alpha1 && rho<alpha2 %R1(i)>=1 && R2(i)>=1
            U(i) = U12_SEIISSICR_v4(param1,param2,mu,b,rho,c,f,'fsolve',opt);
        elseif rho>=alpha1 && rho<alpha2 %R1(i)<1 && R2(i)>=1
            U(i) = U_SICR_v4(param2,mu,b,rho,c,f);
        elseif rho>=alpha2 && rho<alpha1 %R1(i)>=1 && R2(i)<1
            U(i) = U_SEIIS_v4(param1,mu,b,rho,c,f);
        elseif rho>=alpha2 && rho>=alpha1
            U(i) = -c*rho;%max(alpha1,alpha2); %instead of -c*rho ?
        end
        i=i+1;
    end
    %P=0;
    dU  = 0;
end