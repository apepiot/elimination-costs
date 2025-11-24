function [U12,dU,P] = U12_SEIISSEIIIS_v4(param1,param2,mu,b,rho,c,f,solveMethod,optSolver)
    %     syph STI
    % [ X1, 1, 1]
    % [ X2, 2, 1]
    % [ X3, 3, 1]
    % [ X4, 4, 1]
    % [ X5, 5, 1]
    % [ X6, 1, 2]
    % [ X7, 2, 2]
    % [ X8, 3, 2]
    % [ X9, 4, 2]
    % [X10, 5, 2]
    % [X11, 1, 3]
    % [X12, 2, 3]
    % [X13, 3, 3]
    % [X14, 4, 3]
    % [X15, 5, 3]
    % [X16, 1, 4]
    % [X17, 2, 4]
    % [X18, 3, 4]
    % [X19, 4, 4]
    % [X20, 5, 4]
    beta1=param1.beta;gamma1=param1.gamma;nu1=param1.nu;sigma1=param1.sigma;eps1=param1.eps;
    beta2=param2.beta;sigma2=param2.sigma;tau2=param2.tau;gamma12=param2.gamma1;theta2=param2.theta;gamma32=param2.gamma3;nu2=param2.nu;
    alpha1=param1.alpha;
    alpha2=param2.alpha;
    rho = min([rho,alpha1,alpha2]);
    Y0 = ones(20,1);
    if(isequal(solveMethod,'ode45'))    
        options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
        [~,ES_t] = ode45(@(t,Y)  ODE_SEIISSEIIIS_v4(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                     beta2,sigma2,tau2,nu2,gamma12,theta2,gamma32,mu,b,rho),...
                                     tspan,Y0, options);
        ES = ES_t(end,:);
    end
    if(isequal(solveMethod,'fsolve'))
        options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
        [ES] = fsolve(@(Y)  ODE_SEIISSEIIIS_v4(0,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                     beta2,sigma2,tau2,nu2,gamma12,theta2,gamma32,mu,b,rho),...
                                     Y0,options);
    end
    P = max(sum(ES(2:15))/(b/mu),0); %prevalence does not include STI symptomatics
    %P = max(sum(ES(3:5,8:20))/(b/mu),0); %prevalence without exposed ppl
    %U12 = min(max(rho,0),max([alpha1,alpha2]))*(P-c);
    %U12 = min(rho,min([alpha1,alpha2]))*(P-c) ;
    U12 = rho*(P-c); 
    dU  = 0;
end