function [rhohat, rhohat2d] = findRhohat2d(c,type1,type2,rhohat1,rhohat2,param1,param2,mu,b,alphan,alpham);
%findRhohat2d(c,modelTypes(n),modelTypes(m),vecRhohat1d(i,n),vecRhohat1d(i,m),listParam{n},listParam{m},mu,b,alphan,alpham);

%retrouver les bons parametres, la foncitona  été écrasée



%This function finds the argmax of U in a 2-diseases model. 
    options = optimset('Display','off'); %options for minsearch
   
    %to do : reordonner l'ordre des diseases : si type1=SIR et type 2=SIS
    %par exemple
    
    %SIS^2 model
    if strcmp(type1,"SIS") && strcmp(type2,"SIS")
        if (s1==1 && s2==1)
            %%
            r = beta1*beta2/(gamma1*gamma2);
            rhohat12 = (beta1+beta2-mu)-sqrt((beta1+beta2-mu)*(beta1+beta2).*((beta1/gamma2+beta2/gamma1+r+1)./(beta1/gamma2+beta2/gamma1+(2-c)*r)));
        else
            fun12 = @(rho) -U12_SIS2(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
            rhohat12 = min(max(fminsearch(fun12,0,options),0),max(alpha1,alpha2));
        end
        list_rhohat12 = max([rhohat1,rhohat2,rhohat12],0);
        U_list = U_SIS2(list_rhohat12,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
        [Umax,imax] = max(U_list);
        rhohatSIS2 = list_rhohat12(imax);
        rhohat=rhohatSIS2;rhohat2d=rhohat12;
        
    %SISxSIR model
    elseif strcmp(type1,"SIS") && strcmp(type2,"SIR")
        fun21 = @(rho) -U12_SIRSIS7(rho,beta2,beta1,gamma2,gamma1,s2,s1,b,mu,c); %attention a l'ordre des parametres
        rhohat21 = min(max(fminsearch(fun21,0,options),0),max(alpha1,alpha2)); %argmaxU12
        list_rhohat21 = max([rhohat1,rhohat2,rhohat21],0);
        U_list = U_SIRSIS7(list_rhohat21,beta2,beta1,gamma2,gamma1,s2,s1,b,mu,c);
        [Umax,imax] = max(U_list);
        rhohatSISSIR = list_rhohat21(imax);   
        rhohat=rhohatSISSIR;rhohat2d=rhohat21;
    %SIR^2 model
    elseif strcmp(type1,"SIR") && strcmp(type2,"SIR")
    
    %SISxSICAT
    elseif strcmp(type1,"SIS") && strcmp(type2,"SICAT")
    
    %SIRxSICAT
    elseif strcmp(type1,"SIR") && strcmp(type2,"SICAT")
    
    %SICAT^2 model
    elseif strcmp(type1,"SICAT") && strcmp(type2,"SICAT")
        
    else 
        error('Error. Models sorting. Cannot find any combination of models')
    end
end

