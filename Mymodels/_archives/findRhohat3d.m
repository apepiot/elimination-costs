function [rhohat, rhohat3d] = findRhohat3d(c,type1,type2,type3,rhohat1d,rhohat2d,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,mu,b,omega1,omega2,sigma1,sigma2,alpha1,alpha2,alpha3)
% rhohat1 : rhohat for 1disease model (1xN) of cost c,
% rhohat2 : rhohats of 2-disease models (1xN2)
%This function finds the argmax of U in a 2-diseases model. 
    options = optimset('Display','off'); %options for minsearch
    
    %SIS^3 model
    if strcmp(type1,"SIS") && strcmp(type2,"SIS") && strcmp(type3,"SIS")
        fun123 = @(rho) -U123_SIS3(rho,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c);
        rhohat123 = min(max(fminsearch(fun123,0,options),0),max([alpha1,alpha2,alpha3]));
        list_rhohat123 = max([rhohat1d,rhohat2d,rhohat123],0);
        U_list = U_SIS3(list_rhohat123,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c);
        [Umax,imax] = max(U_list);
        rhohatSIS3 = list_rhohat123(imax);
        rhohat=rhohatSIS3;rhohat3d=rhohat123;

    %SIS^2xSIR model
    elseif strcmp(type1,"SIS") && strcmp(type2,"SIS") && strcmp(type3,"SIR")
        fun = @(rho) -U123_SIRSIS2(rho,beta3,beta1,beta2,gamma3,gamma1,gamma2,s3,s1,s2,b,mu,c); %attention a l'ordre des parametres
        rhohat123 = min(max(fminsearch(fun,0,options),0),max([alpha1,alpha2,alpha3])); %argmaxU12
        list_rhohat123 = max([rhohat1d,rhohat2d,rhohat123],0);
        U_list = U_SIRSIS2(list_rhohat123,beta3,beta1,beta2,gamma3,gamma1,gamma2,s3,s1,s2,b,mu,c);
        [Umax,imax] = max(U_list);
        rhohatSIS2SIR = list_rhohat123(imax);   
        rhohat=rhohatSIS2SIR;rhohat3d=rhohat123;
    else 
        error('Error. Models sorting. Cannot find any combination of models')
    end
end

