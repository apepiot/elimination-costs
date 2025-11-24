function [rhohat, rhohat3d] = findRhohat3d_v2(c,type1,type2,type3,...
                                                rhohat1d,rhohat2d,param1,param2,param3,mu,b,...
                                                alpha1,alpha2,alpha3,f)
    % rhohat1 : rhohat for 1disease model (1xN) of cost c,
    % rhohat2 : rhohats of 2-disease models (1xN2)
    %This function finds the argmax of U in a 3-diseases model. 
    options = optimset('Display','off'); %options for minsearch  
    %SEIISXSICR model
    if strcmp(type1,"SEIIS") && strcmp(type2,"SEIIS") && strcmp(type3,"SICR")
        beta1=param1(1);gamma1=param1(2);nu1=param1(3);sigma1=param1(4);eps1=param1(5);
        beta2=param2(1);gamma2=param2(2);nu2=param2(3);sigma2=param2(4);eps2=param2(5);
        betaI3=param3(1);betaC3=param3(2);gamma3=param3(3);sigma3=param3(4);theta3=param3(5);
       
        fun123 = @(rho) -U123_SEIIS2SICR(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,...
            betaI3,betaC3,gamma3,sigma3,theta3,mu,b,rho,c,alpha1,alpha2,alpha3,f);
        rhohat123 = min(max(fminsearch(fun123,0,options),0),max([alpha1,alpha2,alpha3]));
        list_rhohat123 = max([rhohat1d,rhohat2d,rhohat123],0);
        U_list = U_SEIIS2SICR(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,...
            betaI3,betaC3,gamma3,sigma3,theta3,mu,b,list_rhohat123,c,alpha1,alpha2,alpha3,f);
        [~,imax] = max(U_list);
        rhohat = list_rhohat123(imax);
        rhohat3d=rhohat123;
    %SIS^2xSIR model
    elseif strcmp(type1,"SEIIS") && strcmp(type2,"SEIIS") && strcmp(type3,"SEIIIS")
        beta1=param1(1);gamma1=param1(2);nu1=param1(3);sigma1=param1(4);eps1=param1(5);
        beta2=param2(1);gamma2=param2(2);nu2=param2(3);sigma2=param2(4);eps2=param2(5);
        beta3=param3(1);sigma3=param3(2);tau3=param3(3);gamma13=param3(4);theta3=param3(5);gamma33=param3(6);nu3=param3(7);
        fun = @(rho) -U123_SEIIS2SEIIIS(beta1,gamma1,nu1,sigma1,eps1,...
                                     beta2,gamma2,nu2,sigma2,eps2,...
                                     beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,...
                                     mu,b,rho,c,alpha1,alpha2,alpha3,f);
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
        beta1=param1(1);gamma1=param1(2);nu1=param1(3);sigma1=param1(4);eps1=param1(5);
        betaI2=param2(1);betaC2=param2(2);gamma2=param2(3);sigma2=param2(4);theta2=param2(5);
        beta3=param3(1);sigma3=param3(2);tau3=param3(3);gamma13=param3(4);theta3=param3(5);gamma33=param3(6);nu3=param3(7);
        
        fun = @(rho) -U123_SEIISSICRSEIIIS(beta1,gamma1,nu1,sigma1,eps1,...
                                     betaI2,betaC2,gamma2,sigma2,theta2,...
                                     beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,...
                                     mu,b,rho,c,alpha1,alpha2,alpha3,f);
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

