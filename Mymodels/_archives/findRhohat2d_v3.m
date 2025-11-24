function [rhohat, rhohat2d] = findRhohat2d_v3(c,type1,type2,rhohat1,rhohat2,param1,param2,mu,b,alpha1,alpha2,f)
%findRhohat2d(c,modelTypes(n),modelTypes(m),vecRhohat1d(i,n),vecRhohat1d(i,m),listParam{n},listParam{m},mu,b,alphan,alpham);
%This function finds the argmax of U in a 2-diseases model.  
options = optimset('Display','off'); %options for minsearch
    opt.TolFun=1e-5;
    %SEIIS^2 model
    if strcmp(type1,"SEIIS") && strcmp(type2,"SEIIS")
        beta1=param1.beta;gamma1=param1.gamma;nu1=param1.nu;sigma1=param1.sigma;eps1=param1.eps;
        beta2=param2.beta;gamma2=param2.gamma;nu2=param2.nu;sigma2=param2.sigma;eps2=param2.eps;
        fun12 = @(rho) -U12_SEIIS2_v2(beta1,gamma1,nu1,sigma1,eps1,...
                                   beta2,gamma2,nu2,sigma2,eps2,b,mu,rho,c,alpha1,alpha2,f,'fsolve',opt);
        rhohat12 = min(max(fminsearch(fun12,0,options),0),max(alpha1,alpha2));
        list_rhohat12 = max([rhohat1,rhohat2,rhohat12],0);
        U_list = U_SEIIS2_v2(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,b,mu,list_rhohat12,c,alpha1,alpha2,f);        
        [~,imax] = max(U_list);
        rhohatSEIIS2 = list_rhohat12(imax);
        rhohat = rhohatSEIIS2; rhohat2d=rhohat12;
        
    %SEIISxSICR model   
    elseif strcmp(type1,"SEIIS") && strcmp(type2,"SICR")
        beta1=param1.beta;gamma1=param1.gamma;nu1=param1.nu;sigma1=param1.sigma;eps1=param1.eps;
        betaI2=param2.betaI;betaC2=param2.betaC;gamma2=param2.gamma;sigma2=param2.sigma;theta2=param2.theta;        
        fun21 = @(rho) -U12_SEIISSICR(beta1,gamma1,nu1,sigma1,eps1,betaI2,betaC2,gamma2,sigma2,theta2,mu,b,rho,c,alpha1,alpha2,f,'fsolve',opt);
        rhohat21 = min(max(fminsearch(fun21,0,options),0),max(alpha1,alpha2)); %argmaxU12        
        list_rhohat21 = max([rhohat1,rhohat2,rhohat21],0);
        U_list = U_SEIISSICR(beta1,gamma1,nu1,sigma1,eps1,betaI2,betaC2,gamma2,sigma2,theta2,mu,b,list_rhohat21,c,alpha1,alpha2,f);
        [~,imax] = max(U_list);
        rhohatSEIISSICR = list_rhohat21(imax);   
        rhohat=rhohatSEIISSICR;rhohat2d=rhohat21;
        
    %SEIISxSEIIIS model
    elseif strcmp(type1,"SEIIS") && strcmp(type2,"SEIIIS")
        beta1=param1.beta;gamma1=param1.gamma;nu1=param1.nu;sigma1=param1.sigma;eps1=param1.eps;
        beta2=param2.beta;sigma2=param2.sigma;tau2=param2.tau;gamma12=param2.gamma1;theta2=param2.theta;gamma32=param2.gamma3;nu2=param2.nu;
        fun21 = @(rho) -U12_SEIISSEIIIS(beta1,gamma1,nu1,sigma1,eps1,beta2,sigma2,tau2,theta2,gamma12,gamma32,nu2,mu,b,rho,c,alpha1,alpha2,f,'fsolve',opt);
        rhohat21 = min(max(fminsearch(fun21,0,options),0),max(alpha1,alpha2)); %argmaxU12
        list_rhohat21 = max([rhohat1,rhohat2,rhohat21],0);
        U_list   = U_SEIISSEIIIS(beta1,gamma1,nu1,sigma1,eps1,beta2,sigma2,tau2,theta2,gamma12,gamma32,nu2,mu,b,list_rhohat21,c,alpha1,alpha2,f);
        [~,imax] = max(U_list);
        rhohat   = list_rhohat21(imax);   
        rhohat2d = rhohat21;
    %SICRxSEIIIS
    elseif strcmp(type1,"SICR") && strcmp(type2,"SEIIIS")
        betaI1=param1.betaI;betaC1=param1.betaC;gamma1=param1.gamma;sigma1=param1.sigma;theta1=param1.theta;
        beta2=param2.beta;sigma2=param2.sigma;tau2=param2.tau;gamma12=param2.gamma1;theta2=param2.theta;gamma32=param2.gamma3;nu2=param2.nu;
        
        fun21 = @(rho) -U12_SICRSEIIIS(betaI1,betaC1,gamma1,sigma1,theta1,...
                                       beta2,sigma2,tau2,nu2,gamma12,theta2,gamma32,...
                                       mu,b,rho,c,alpha1,alpha2,f,'fsolve',opt);%,...
                                       %eta1,omega1,eta2,omega2);
        rhohat21 = min(max(fminsearch(fun21,0,options),0),max(alpha1,alpha2)); %argmaxU12
        list_rhohat21 = max([rhohat1,rhohat2,rhohat21],0);
        U_list   = U_SICRSEIIIS(betaI1,betaC1,gamma1,sigma1,theta1,beta2,sigma2,tau2,nu2,gamma12,theta2,gamma32,mu,b,...
                                list_rhohat21,c,alpha1,alpha2,f);%,eta1,omega1,eta2,omega2);
        [~,imax] = max(U_list);
        rhohat   = list_rhohat21(imax);   
        rhohat2d = rhohat21;   
    else 
        error('Error. Models sorting. Cannot find any combination of models')
    end
end

