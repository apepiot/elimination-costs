function [rhohat] = findRhohat2d_v7(vecC,type1,type2,param1,param2,paramRho,mu,b,f,kit,opt)
%findRhohat2d(c,modelTypes(n),modelTypes(m),vecRhohat1d(i,n),vecRhohat1d(i,m),listParam{n},listParam{m},mu,b,alphan,alpham);
%diff with findRhohat2d_v4: SICR <- SICTP
%This function finds the argmax of U in a 2-diseases model.  
%alpha1 = param1.alpha; 
%alpha2 = param2.alpha;
options = optimset('Display','off'); %options for minsearch
solveWith='knitro-ampl';


%SEIIS^2 model
if strcmp(type1,"SEIIS") && strcmp(type2,"SEIIS")
    fun12 = @(rho) -U12_SEIIS2_v6(param1,param2,mu,b,rho,c,f,'fsolve',opt);
    %recuperer P12
    rhohat12 = fminsearch(fun12,0,options);
    rhohat12 = min(max(rhohat12,0),max(alpha1,alpha2)); %a voir, si on enelve cette ligne
    
    %rhohat12 = fmincon(fun12,0,A,bcon, 'Display', 'off');%,options);
    list_rhohat12 = max([rhohat1,rhohat2,rhohat12],0);
    U_list = U_SEIIS2_v4(param1,param2,mu,b,list_rhohat12,c,f);        
    [~,imax] = max(U_list);
    rhohatSEIIS2 = list_rhohat12(imax);
    rhohat = rhohatSEIIS2; rhohat2d=rhohat12;

%SEIISxSICTP model   
elseif strcmp(type1,"SICTP") && strcmp(type2,"SEIIS")     
    i=1;
    rhohat = zeros(length(vecC),1);
    for c=vecC
        if strcmp(solveWith,'matlab')
            fun21 = @(rho) -U12_SICTPSEIIS_v7(param1,param2,mu,b,paramRho,rho,c,f,'knitro-ampl',opt);
            [rhohat21,Cval,exitflag,output] = fminsearch(fun21,0,options);
        elseif strcmp(solveWith,'knitro-ampl')
            paramTab = {param1,param2};
            
            [rhohat21,Cval,ES,msg] = fminU_knitro_v7_bis(kit,mod,paramTab,paramRho,mu,b,c,opt);
        end
        rhohat(i)=rhohat21;
        i=i+1;        
    end
    %je veux recuperer Pval, car si Pval =0, on est a la limite et on passe
    %au modele a 1 infection
    %rhohat21 = min(max(fminsearch(fun21,0,options),0),max(alpha1,alpha2)); %argmaxU12     
    %rhohat21 = fmincon(fun21,0,A,bcon, 'Display', 'off');
%     list_rhohat21 = max([rhohat1,rhohat2,rhohat21],0);
%     U_list = U_SICTPSEIIS_v1(param2,param1,mu,b,list_rhohat21,c,f);
%     [~,imax] = max(U_list);
%     rhohatSICTPSEIIS = list_rhohat21(imax);   
%     rhohat=rhohatSICTPSEIIS;rhohat2d=rhohat21;
    

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
    
%SICTPxSEIIIS
elseif strcmp(type1,"SICTP") && strcmp(type2,"SEIIIS")
    fun21 = @(rho) -U12_SICTPSEIIIS_v1(param1,param2,mu,b,rho,c,f,'fsolve',opt);
    rhohat21 = min(max(fminsearch(fun21,0,options),0),max(alpha1,alpha2)); %argmaxU12
    %rhohat21 = fmincon(fun21,0,A,bcon, 'Display', 'off');
    list_rhohat21 = max([rhohat1,rhohat2,rhohat21],0);
    U_list   = U_SICTPSEIIIS_v1(param1,param2,mu,b,list_rhohat21,c,f);
    [~,imax] = max(U_list);
    rhohat   = list_rhohat21(imax);   
    rhohat2d = rhohat21;   
else 
    error('Error. Models sorting. Cannot find any combination of models')
end
end

