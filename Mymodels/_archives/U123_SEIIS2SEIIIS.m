function [U123,dU,P] = U123_SEIIS2SEIIIS(beta1,gamma1,nu1,sigma1,eps1,...
                                     beta2,gamma2,nu2,sigma2,eps2,...
                                     beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,...
                                     mu,b,rho,c,alpha1,alpha2,alpha3,f,solveMethod,optSolver)                                
    
    %syph, STI1, STI2                            
    % [ X1, 1, 1, 1]
    % [ X2, 2, 1, 1]
    % [ X3, 3, 1, 1]
    % [ X4, 4, 1, 1]
    % [ X5, 5, 1, 1]
    % [ X6, 1, 2, 1]
    % [ X7, 2, 2, 1]
    % [ X8, 3, 2, 1]
    % [ X9, 4, 2, 1]
    % [X10, 5, 2, 1]
    % [X11, 1, 3, 1]
    % [X12, 2, 3, 1]
    % [X13, 3, 3, 1]
    % [X14, 4, 3, 1]
    % [X15, 5, 3, 1]
    % [X16, 1, 4, 1]
    % [X17, 2, 4, 1]
    % [X18, 3, 4, 1]
    % [X19, 4, 4, 1]
    % [X20, 5, 4, 1]
    % [X21, 1, 1, 2]
    % [X22, 2, 1, 2]
    % [X23, 3, 1, 2]
    % [X24, 4, 1, 2]
    % [X25, 5, 1, 2]
    % [X26, 1, 2, 2]
    % [X27, 2, 2, 2]
    % [X28, 3, 2, 2]
    % [X29, 4, 2, 2]
    % [X30, 5, 2, 2]
    % [X31, 1, 3, 2]
    % [X32, 2, 3, 2]
    % [X33, 3, 3, 2]
    % [X34, 4, 3, 2]
    % [X35, 5, 3, 2]
    % [X36, 1, 4, 2]
    % [X37, 2, 4, 2]
    % [X38, 3, 4, 2]
    % [X39, 4, 4, 2]
    % [X40, 5, 4, 2]
    % [X41, 1, 1, 3]
    % [X42, 2, 1, 3]
    % [X43, 3, 1, 3]
    % [X44, 4, 1, 3]
    % [X45, 5, 1, 3]
    % [X46, 1, 2, 3]
    % [X47, 2, 2, 3]
    % [X48, 3, 2, 3]
    % [X49, 4, 2, 3]
    % [X50, 5, 2, 3]
    % [X51, 1, 3, 3]
    % [X52, 2, 3, 3]
    % [X53, 3, 3, 3]
    % [X54, 4, 3, 3]
    % [X55, 5, 3, 3]
    % [X56, 1, 4, 3]
    % [X57, 2, 4, 3]
    % [X58, 3, 4, 3]
    % [X59, 4, 4, 3]
    % [X60, 5, 4, 3]
    % [X61, 1, 1, 4]
    % [X62, 2, 1, 4]
    % [X63, 3, 1, 4]
    % [X64, 4, 1, 4]
    % [X65, 5, 1, 4]
    % [X66, 1, 2, 4]
    % [X67, 2, 2, 4]
    % [X68, 3, 2, 4]
    % [X69, 4, 2, 4]
    % [X70, 5, 2, 4]
    % [X71, 1, 3, 4]
    % [X72, 2, 3, 4]
    % [X73, 3, 3, 4]
    % [X74, 4, 3, 4]
    % [X75, 5, 3, 4]
    % [X76, 1, 4, 4]
    % [X77, 2, 4, 4]
    % [X78, 3, 4, 4]
    % [X79, 4, 4, 4]
    % [X80, 5, 4, 4]
    rho = min([rho,alpha1,alpha2,alpha3]);
    Y0 = ones(80,1);
    if(isequal(solveMethod,'ode45'))
        tspan = optSolver.tspan;
        options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
        [~,ES_t] = ode45(@(t,Y)  ODE_SEIIS2SEIIIS_3(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                 beta2,nu2,eps2,sigma2,gamma2,...
                                 beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,mu,b,rho),...
                                 tspan,Y0, options);
        ES = ES_t(end,:);
    end
    if(isequal(solveMethod,'fsolve'))
        options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
        [ES] = fsolve(@(Y)  ODE_SEIIS2SEIIIS_3(1,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                 beta2,nu2,eps2,sigma2,gamma2,...
                                 beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,mu,b,rho),...
                                 Y0,options);    
    end
        
    %P = max(sum(ES(end,2:80))/(b/mu),0); %v2
    P = max(sum(ES([2:15,21:35,41:55]))/(b/mu),0);
    %U123 = min(max(rho,0),max([alpha1,alpha2,alpha3]))*(P-c);
    %U123 = min(rho,min([alpha1,alpha2,alpha3]))*(P-c) ;
    U123 = rho*(P-c);
    dU  = 0;
end