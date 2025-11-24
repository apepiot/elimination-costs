function [U] = U_SEIIS2SEIIIS_v4(param1,param2,param3,mu,b,vecRho,c,f)   
    alpha1=param1.alpha;
    alpha2=param2.alpha;
    alpha3=param3.alpha;
    i = 1;
    U = zeros(1,length(vecRho));opt.TolFun=1e-5;
    for rho=vecRho
        if rho<min([alpha1,alpha2,alpha3])%R1(i)>1 && R2(i)>=1 && R3(i)>=1
            U(i) = U123_SEIIS2SEIIIS_v4(param1,param2,param3,mu,b,rho,c,f,'fsolve',opt);
        elseif rho>=alpha1 %R1(i)<=1 && R2(i)>1 && R3(i)>1
            U(i) = U_SEIISSEIIIS_v4(param2,param3,mu,b,rho,c,f);
        elseif rho>=alpha2 %R1(i)>1 && R2(i)<=1 && R3(i)>1
            U(i) = U_SEIISSEIIIS_v4(param1,param3,mu,b,rho,c,f);
        elseif rho>=alpha3 %R1(i)>1 && R2(i)>1 && R3(i)<=1
            U(i) = U12_SEIIS2_v4(param1,param2,mu,b,rho,c,f,'fsolve',opt);
        end
        i=i+1;
    end
end