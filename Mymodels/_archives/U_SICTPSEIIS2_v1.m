function [U] = U_SICTPSEIIS2_v1(param1,param2,param3,mu,b,vecRho,c,f)                                 
    U  = zeros(1,length(vecRho)); opt.TolFun=1e-5;
    alpha1=param1.alpha;alpha2=param2.alpha;alpha3=param3.alpha;
    i = 1;
    for rho=vecRho
        if rho<alpha1 && rho<alpha2 && rho<alpha3 %R1(i)>1 && R2(i)>1 && R3(i)>1
            U(i) = U123_SICTPSEIIS2_v1(param1,param2,param3,mu,b,rho,c,f,'fsolve',opt);
        elseif rho>=alpha1 %R1(i)<=1 && R2(i)>1 && R3(i)>1
            U(i) = U_SEIIS2_v4(param2,param3,mu,b,rho,c,f);
        elseif rho>=alpha2 %R1(i)>1 && R2(i)<=1 && R3(i)>1
            U(i) = U_SICTPSEIIS_v1(param1,param3,mu,b,rho,c,f);
        elseif rho>=alpha3 %R1(i)>1 && R2(i)>1 && R3(i)<=1
            U(i) = U_SICTPSEIIS_v1(param1,param2,mu,b,rho,c,f); 
        end
        i=i+1;
    end
end