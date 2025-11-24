function [rhohat,rhohat4d] = findRhohat4d_v4(c,type1,type2,type3,type4,...
                                                rhohat1d,rhohat2d,rhohat3d,...
                                                param1,param2,param3,param4,mu,b,f)

    %This function finds the argmax of U in a 3-diseases model. 
    alpha1=param1.alpha; alpha2=param2.alpha;
    alpha3=param3.alpha; alpha4=param4.alpha;
    options = optimset('Display','off','TolFun',1e-2); %options for minsearch  
    if strcmp(type1,"SEIIS") && strcmp(type2,"SEIIS") && strcmp(type3,"SICR") && strcmp(type4,"SEIIIS")
        %disp(['avant ode system']); toc
        opt.TolFun=1e-5;
        fun1234 = @(rho) -U1234_SEIIS2SICRSEIIIS_v4(param1,param2,param3,param4,mu,b,rho,c,f,'fsolve',opt);
        
        rhohat1234 = min(max(fminsearch(fun1234,0,options),0),max([alpha1,alpha2,alpha3,alpha4]));
        %disp(['apres ode system']); toc
        list_rhohat1234 = max([rhohat1d,rhohat2d,rhohat3d,rhohat1234],0);
        U_list = U_SEIIS2SICRSEIIIS_v4(param1,param2,param3,param4,...
                                    mu,b,list_rhohat1234,c,f);

        [~,imax] = max(U_list);
        rhohat = list_rhohat1234(imax);
        rhohat4d = rhohat1234;            
    else 
        error('Error. Models sorting. Cannot find any combination of models')
    end
end

