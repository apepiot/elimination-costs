function [rhohat,c0n,cnn] = findRhohat1d_v7(type,param,mu,b,vecC,f)
    %find rhohat for a 1-disease model
    options = optimset('Display','off','TolFun',1e-2); %options for minsearch 
    lC = length(vecC);
    rhohat = zeros(lC,1);
    if isequal(type,'SEIIS')
        [~,c0n] = U1_SEIISv4(param,mu,b,0,0,f); %c2 = c0
        [~,cnn] = U1_SEIISv4(param,mu,b,param.alpha,0,f); %c1=c1
        i=1;     
        for c=vecC
            fun = @(rho) -U1_SEIISv4(param,mu,b,rho,c,f);
            rhohat(i) = min(max(fminsearch(fun,0,options),0),param.alpha);
            i=i+1; %c
        end
    elseif isequal(type,'SICTP')
        [~,c0n] = U1_SICTP(param,mu,b,0,0,f); %c2 = c0
        [~,cnn] = U1_SICTP(param,mu,b,param.alpha_prep,0,f); %c1=c1
        i=1;
        for c=vecC 
            fun = @(rho) -U1_SICTP(param,mu,b,rho,c,f);
            rhohat(i) = min(max(fminsearch(fun,0,options),0),param.alpha_prev);
            i=i+1; %c
        end
    elseif isequal(type,'SEIIIS')
        [~,c0n] = U1_SEIIIS_v4(param,mu,b,0,0,f);
        [~,cnn] = U1_SEIIIS_v4(param,mu,b,param.alpha,0,f);
        i=1;
        for c=vecC
            fun = @(rho) -U1_SEIIIS_v4(param,mu,b,rho,c,f);
            rhohat(i) = min(max(fminsearch(fun,0,options),0),param.alpha);
            i=i+1; %c
        end
    end
end

