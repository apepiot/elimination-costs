function [rhohat, rhohat3d] = findRhohat3d_v3(c,type1,type2,type3,...
                                                rhohat1d,rhohat2d,param1,param2,param3,mu,b,...
                                                alpha1,alpha2,alpha3,f)
    opt.TolFun=1e-5;
    % rhohat1 : rhohat for 1disease model (1xN) of cost c,
    % rhohat2 : rhohats of 2-disease models (1xN2)
    %This function finds the argmax of U in a 3-diseases model. 
    options = optimset('Display','off'); %options for minsearch  
    %SEIISXSICR model
    if strcmp(type1,"SEIIS") && strcmp(type2,"SEIIS") && strcmp(type3,"SICR")
        beta1=param1.beta;gamma1=param1.gamma;nu1=param1.nu;sigma1=param1.sigma;eps1=param1.eps;
        beta2=param2.beta;gamma2=param2.gamma;nu2=param2.nu;sigma2=param2.sigma;eps2=param2.eps;
        betaI3=param3.betaI;betaC3=param3.betaC;gamma3=param3.gamma;sigma3=param3.sigma;theta3=param3.theta;
       
        fun123 = @(rho) -U123_SEIIS2SICR(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,...
            betaI3,betaC3,gamma3,sigma3,theta3,mu,b,rho,c,alpha1,alpha2,alpha3,f,'fsolve',opt);
        rhohat123 = min(max(fminsearch(fun123,0,options),0),max([alpha1,alpha2,alpha3]));
        list_rhohat123 = max([rhohat1d,rhohat2d,rhohat123],0);
        U_list = U_SEIIS2SICR(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,...
            betaI3,betaC3,gamma3,sigma3,theta3,mu,b,list_rhohat123,c,alpha1,alpha2,alpha3,f);
        [~,imax] = max(U_list);
        rhohat = list_rhohat123(imax);
        rhohat3d=rhohat123;
    %SIS^2xSIR model
    elseif strcmp(type1,"SEIIS") && strcmp(type2,"SEIIS") && strcmp(type3,"SEIIIS")
        beta1=param1.beta;gamma1=param1.gamma;nu1=param1.nu;sigma1=param1.sigma;eps1=param1.eps;
        beta2=param2.beta;gamma2=param2.gamma;nu2=param2.nu;sigma2=param2.sigma;eps2=param2.eps;
        beta3=param3.beta;sigma3=param3.sigma;tau3=param3.tau;gamma13=param3.gamma1;theta3=param3.theta;gamma33=param3.gamma3;nu3=param3.nu;
        fun = @(rho) -U123_SEIIS2SEIIIS(beta1,gamma1,nu1,sigma1,eps1,...
                                     beta2,gamma2,nu2,sigma2,eps2,...
                                     beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,...
                                     mu,b,rho,c,alpha1,alpha2,alpha3,f,'fsolve',opt);
        rhohat123 = min(max(fminsearch(fun,0,options),0),max([alpha1,alpha2,alpha3])); %argmaxU12
        list_rhohat123 = max([rhohat1d,rhohat2d,rhohat123],0);
        U_list = U_SEIIS2SEIIIS(beta1,gamma1,nu1,sigma1,eps1,...
                                     beta2,gamma2,nu2,sigma2,eps2,...
                                     beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,...
                                     mu,b,list_rhohat123,c,alpha1,alpha2,alpha3,f);
        [~,imax] = max(U_list);
        rhohat = list_rhohat123(imax);   
        rhohat3d=rhohat123;
        
    elseif strcmp(type1,"SEIIS") && strcmp(type2,"SICR") && strcmp(type3,"SEIIIS")  
        beta1=param1.beta;gamma1=param1.gamma;nu1=param1.nu;sigma1=param1.sigma;eps1=param1.eps;
        betaI2=param2.betaI;betaC2=param2.betaC;gamma2=param2.gamma;sigma2=param2.sigma;theta2=param2.theta;
        beta3=param3.beta;sigma3=param3.sigma;tau3=param3.tau;gamma13=param3.gamma1;theta3=param3.theta;gamma33=param3.gamma3;nu3=param3.nu;
       
        fun = @(rho) -U123_SEIISSICRSEIIIS(beta1,gamma1,nu1,sigma1,eps1,...
                                     betaI2,betaC2,gamma2,sigma2,theta2,...
                                     beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,...
                                     mu,b,rho,c,alpha1,alpha2,alpha3,f,'fsolve',opt);
        rhohat123 = min(max(fminsearch(fun,0,options),0),max([alpha1,alpha2,alpha3])); %argmaxU12
        list_rhohat123 = max([rhohat1d,rhohat2d,rhohat123],0); %candidates for Umax
        U_list = U_SEIISSICRSEIIIS(beta1,gamma1,nu1,sigma1,eps1,...
                                     betaI2,betaC2,gamma2,sigma2,theta2,...
                                     beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,...
                                     mu,b,list_rhohat123,c,alpha1,alpha2,alpha3,f);
        [~,imax] = max(U_list);
        rhohat   = list_rhohat123(imax);   
        rhohat3d = rhohat123;                         
    else 
        error('Error. Models sorting. Cannot find any combination of models')
    end
end

