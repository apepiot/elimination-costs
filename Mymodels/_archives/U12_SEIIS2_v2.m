function [U12,dU,P] = U12_SEIIS2_v2(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,b,mu,rho,c,alpha1,alpha2,f,solveMethod,optSolver)
    
    rho = min([rho,alpha1,alpha2]);
    Y0 = ones(16,1);
    if(isequal(solveMethod,'ode45'))
        tspan = optSolver.tspan;    
        options = odeset('RelTol',1e-3,'Stats','off');%,'OutputFcn',@odeplot);
        [~,ES_t] = ode45(@(t,Y) ODE_SEIIS2_3(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                     beta2,nu2,eps2,sigma2,gamma2,mu,b,rho),tspan,Y0,options);
        ES = ES_t(end,:);
    end
    if(isequal(solveMethod,'fsolve'))
        options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
        [ES] = fsolve(@(Y)  ODE_SEIIS2_3(1,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                     beta2,nu2,eps2,sigma2,gamma2,mu,b,rho),...
                                      Y0,options);
    end
    %P = max(sum(ES(end,2:16)),0)/(b/mu); %v2
    P = max(sum(ES([2:3,5:7,9:11])),0)/(b/mu); %v2

    %U12 = min([rho,alpha1,alpha2])*(P-c) ;
    U12 = rho*(P-c);
    dU  = 0;
end