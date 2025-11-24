function [U12,dU,P] = U12_SICRSEIIIS(betaI1,betaC1,gamma1,sigma1,theta1,...
                                     beta2,sigma2,tau2,nu2,gamma12,theta2,gamma32,...
                                     mu,b,rho,c,alpha1,alpha2,f,solveMethod,optSolver)

    % [ X1, 1, 1]
    % [ X2, 2, 1]
    % [ X3, 3, 1]
    % [ X4, 4, 1]
    % [ X5, 1, 2]
    % [ X6, 2, 2]
    % [ X7, 3, 2]
    % [ X8, 4, 2]
    % [ X9, 1, 3]
    % [X10, 2, 3]
    % [X11, 3, 3]
    % [X12, 4, 3]
    % [X13, 1, 4]
    % [X14, 2, 4]
    % [X15, 3, 4]
    % [X16, 4, 4]
    % [X17, 1, 5]
    % [X18, 2, 5]
    % [X19, 3, 5]
    % [X20, 4, 5]              
    rho = min([rho,alpha1,alpha2]);
    Y0 = ones(20,1);
    if(isequal(solveMethod,'ode45'))
        tspan = optSolver.tspan;    
        options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
        [~,ES_t] = ode45(@(t,Y) ODE_SICRSEIIIS_3(t,Y,betaI1,betaC1,gamma1,sigma1,theta1,...
                                     beta2,sigma2,tau2,nu2,gamma12,theta2,gamma32,mu,b,rho),...
                                     tspan,Y0, options);
         ES = ES_t(end,:);
    end
    if(isequal(solveMethod,'fsolve'))
        options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
        [ES] = fsolve(@(Y)  ODE_SICRSEIIIS_3(1,Y,betaI1,betaC1,gamma1,sigma1,theta1,...
                                     beta2,sigma2,tau2,nu2,gamma12,theta2,gamma32,mu,b,rho),...
                                     Y0,options);
    end
    
    %P12 = max(sum(ES(end,[2,3,5:20]))/(b/mu),0);
    PHIV = sum(ES([2,3,6,7,10,11,14,15,18,19]))/(b/mu);             
    P12 = f*PHIV + sum(ES([5,8,9,12,13,16,17,20]))/(b/mu);
    %U12 = min(max(rho,0),max([alpha1,alpha2]))*(P-c);
    P = max(P12,0); 
    %U12 = min(rho,min([alpha1,alpha2]))*(P-c) ;
    U12 = rho*(P-c);
    dU  = 0;
end