function [rhohat, rhohat2d] = findRhohat2d_v4(c,type1,type2,rhohat1,rhohat2,param1,param2,mu,b,f)
%findRhohat2d(c,modelTypes(n),modelTypes(m),vecRhohat1d(i,n),vecRhohat1d(i,m),listParam{n},listParam{m},mu,b,alphan,alpham);
%This function finds the argmax of U in a 2-diseases model.  
alpha1 = param1.alpha; 
alpha2 = param2.alpha;
options = optimset('Display','off'); %options for minsearch
opt.TolFun=1e-8;
A = [1;-1];
bcon = [min(alpha1,alpha2);0];
%SEIIS^2 model
if strcmp(type1,"SEIIS") && strcmp(type2,"SEIIS")
    fun12 = @(rho) -U12_SEIIS2_v4(param1,param2,mu,b,rho,c,f,'fsolve',opt);
    rhohat12 = min(max(fminsearch(fun12,0,options),0),max(alpha1,alpha2));
    %rhohat12 = fmincon(fun12,0,A,bcon, 'Display', 'off');%,options);
    list_rhohat12 = max([rhohat1,rhohat2,rhohat12],0);
    U_list = U_SEIIS2_v4(param1,param2,mu,b,list_rhohat12,c,f);        
    [~,imax] = max(U_list);
    rhohatSEIIS2 = list_rhohat12(imax);
    rhohat = rhohatSEIIS2; rhohat2d=rhohat12;

%SEIISxSICR model   
elseif strcmp(type1,"SEIIS") && strcmp(type2,"SICR")     
    fun21 = @(rho) -U12_SEIISSICR_v4(param1,param2,mu,b,rho,c,f,'fsolve',opt);
    rhohat21 = min(max(fminsearch(fun21,0,options),0),max(alpha1,alpha2)); %argmaxU12     
    %rhohat21 = fmincon(fun21,0,A,bcon, 'Display', 'off');
    list_rhohat21 = max([rhohat1,rhohat2,rhohat21],0);
    U_list = U_SEIISSICR_v4(param1,param2,mu,b,list_rhohat21,c,f);
    [~,imax] = max(U_list);
    rhohatSEIISSICR = list_rhohat21(imax);   
    rhohat=rhohatSEIISSICR;rhohat2d=rhohat21;

%SEIISxSEIIIS model
elseif strcmp(type1,"SEIIS") && strcmp(type2,"SEIIIS")
    fun21 = @(rho) -U12_SEIISSEIIIS_v4(param1,param2,mu,b,rho,c,f,'fsolve',opt);
    rhohat21 = min(max(fminsearch(fun21,0,options),0),max(alpha1,alpha2)); %argmaxU12
    %rhohat21 = fmincon(fun21,0,A,bcon, 'Display', 'off');
    list_rhohat21 = max([rhohat1,rhohat2,rhohat21],0);
    U_list   = U_SEIISSEIIIS_v4(param1,param2,mu,b,list_rhohat21,c,f);
    [~,imax] = max(U_list);
    rhohat   = list_rhohat21(imax);   
    rhohat2d = rhohat21;
    
%SICRxSEIIIS
elseif strcmp(type1,"SICR") && strcmp(type2,"SEIIIS")
    fun21 = @(rho) -U12_SICRSEIIIS_v4(param1,param2,mu,b,rho,c,f,'fsolve',opt);
    rhohat21 = min(max(fminsearch(fun21,0,options),0),max(alpha1,alpha2)); %argmaxU12
    %rhohat21 = fmincon(fun21,0,A,bcon, 'Display', 'off');
    list_rhohat21 = max([rhohat1,rhohat2,rhohat21],0);
    U_list   = U_SICRSEIIIS_v4(param1,param2,mu,b,list_rhohat21,c,f);
    [~,imax] = max(U_list);
    rhohat   = list_rhohat21(imax);   
    rhohat2d = rhohat21;   
else 
    error('Error. Models sorting. Cannot find any combination of models')
end
end

