function [U] = U_SICTPSEIIISSEIIS2_v2(param1,param2,param3,param4,...
                                    mu,b,paramRho,c,f)
%     U = zeros(1,length(vecRho));
%     alpha1=param1.alpha;alpha2=param2.alpha;
%     alpha3=param3.alpha;alpha4=param4.alpha;
    i = 1; opt.TolFun=1e-5;
%     for rho=vecRho
%         if rho<min([alpha1,alpha2,alpha3,alpha4])
            U(i) = U1234_SICTPSEIIISSEIIS2_v1(param1,param2,param3,param4,...
                                          mu,b,paramRho,c,f,'fsolve',opt);
%         elseif rho>=alpha1
%             U(i) = U_SEIIS2SEIIIS_v4(param3,param4,param2,mu,b,rho,c,f);
%         elseif rho>=alpha2
%             U(i) = U_SICTPSEIIS2_v1(param1,param3,param4,mu,b,rho,c,f);
%         elseif rho>=alpha3%R3(i)<1 
%             U(i) = U_SICTPSEIIISSEIIS_v1(param1,param2,param4,mu,b,rho,c,f);
%         elseif rho>=alpha4%R4(i)<1
%             U(i) = U_SICTPSEIIISSEIIS_v2(param1,param2,param3,mu,b,rho,c,f);
%         end
%         i=i+1;
%     end
end