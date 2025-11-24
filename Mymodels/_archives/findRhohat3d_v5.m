function [rhohat, rhohat3d] = findRhohat3d_v5(c,type1,type2,type3,...
                                                rhohat1d,rhohat2d,...
                                                param1,param2,param3,mu,b,f)
    opt.TolFun=1e-5;
    % rhohat1 : rhohat for 1disease model (1xN) of cost c,
    % rhohat2 : rhohats of 2-disease models (1xN2)
    %This function finds the argmax of U in a 3-diseases model. 
    options = optimset('Display','off'); %options for minsearch  
    alpha1=param1.alpha;alpha2=param2.alpha;alpha3=param3.alpha;
    
    %SICTPxSEIISXSEIIS model
    if strcmp(type1,"SEIIS") && strcmp(type2,"SEIIS") && strcmp(type3,"SICTP")
        fun123 = @(rho) -U123_SICTPSEIIS2_v1(param3,param1,param2,mu,b,rho,c,f,'fsolve',opt);
        rhohat123 = min(max(fminsearch(fun123,0,options),0),max([alpha1,alpha2,alpha3]));
        list_rhohat123 = max([rhohat1d,rhohat2d,rhohat123],0);
        U_list = U_SICTPSEIIS2_v1(param3,param1,param2,mu,b,list_rhohat123,c,f);
        [~,imax] = max(U_list);
        rhohat = list_rhohat123(imax);
        rhohat3d=rhohat123;
        
    %SEIISxSEIISxSEIIIS model
    elseif strcmp(type1,"SEIIS") && strcmp(type2,"SEIIS") && strcmp(type3,"SEIIIS")
        fun = @(rho) -U123_SEIIS2SEIIIS_v4(param1,param2,param3,mu,b,rho,c,f,'fsolve',opt);
        rhohat123 = min(max(fminsearch(fun,0,options),0),max([alpha1,alpha2,alpha3])); %argmaxU12
        list_rhohat123 = max([rhohat1d,rhohat2d,rhohat123],0);
        U_list = U_SEIIS2SEIIIS_v4(param1,param2,param3,mu,b,list_rhohat123,c,f);
        [~,imax] = max(U_list);
        rhohat = list_rhohat123(imax);
        rhohat3d=rhohat123;
    
    %SICTPxSEIIISxSEIIS model
    elseif strcmp(type1,"SEIIS") && strcmp(type2,"SICTP") && strcmp(type3,"SEIIIS")  
        fun = @(rho) -U123_SICTPSEIIISSEIIS_v1(param2,param3,param1,...
                                     mu,b,rho,c,f,'fsolve',opt);
        rhohat123 = min(max(fminsearch(fun,0,options),0),max([alpha1,alpha2,alpha3])); %argmaxU12
        list_rhohat123 = max([rhohat1d,rhohat2d,rhohat123],0); %candidates for Umax
        U_list = U_SICTPSEIIISSEIIS_v1(param2,param3,param1,mu,b,list_rhohat123,c,f);
        [~,imax] = max(U_list);
        rhohat   = list_rhohat123(imax);
        rhohat3d = rhohat123; 
    else 
        error('Error. Models sorting. Cannot find any combination of models')
    end
end

