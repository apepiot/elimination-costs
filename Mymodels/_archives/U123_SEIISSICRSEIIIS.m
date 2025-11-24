function [U123,dU,P] = U123_SEIISSICRSEIIIS(beta1,gamma1,nu1,sigma1,eps1,...
                                     betaI2,betaC2,gamma2,sigma2,theta2,...
                                     beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,...
                                     mu,b,rho,c,alpha1,alpha2,alpha3,f,solveMethod,optSolver)
    % HIV,Syphilis,STI                        
    % [ X1, 1, 1, 1]
    % [ X2, 2, 1, 1]
    % [ X3, 3, 1, 1]
    % [ X4, 4, 1, 1]
    % [ X5, 1, 2, 1]
    % [ X6, 2, 2, 1]
    % [ X7, 3, 2, 1]
    % [ X8, 4, 2, 1]
    % [ X9, 1, 3, 1]
    % [X10, 2, 3, 1]
    % [X11, 3, 3, 1]
    % [X12, 4, 3, 1]
    % [X13, 1, 4, 1]
    % [X14, 2, 4, 1]
    % [X15, 3, 4, 1]
    % [X16, 4, 4, 1]
    % [X17, 1, 5, 1]
    % [X18, 2, 5, 1]
    % [X19, 3, 5, 1]
    % [X20, 4, 5, 1]
    % [X21, 1, 1, 2]
    % [X22, 2, 1, 2]
    % [X23, 3, 1, 2]
    % [X24, 4, 1, 2]
    % [X25, 1, 2, 2]
    % [X26, 2, 2, 2]
    % [X27, 3, 2, 2]
    % [X28, 4, 2, 2]
    % [X29, 1, 3, 2]
    % [X30, 2, 3, 2]
    % [X31, 3, 3, 2]
    % [X32, 4, 3, 2]
    % [X33, 1, 4, 2]
    % [X34, 2, 4, 2]
    % [X35, 3, 4, 2]
    % [X36, 4, 4, 2]
    % [X37, 1, 5, 2]
    % [X38, 2, 5, 2]
    % [X39, 3, 5, 2]
    % [X40, 4, 5, 2]
    % [X41, 1, 1, 3]
    % [X42, 2, 1, 3]
    % [X43, 3, 1, 3]
    % [X44, 4, 1, 3]
    % [X45, 1, 2, 3]
    % [X46, 2, 2, 3]
    % [X47, 3, 2, 3]
    % [X48, 4, 2, 3]
    % [X49, 1, 3, 3]
    % [X50, 2, 3, 3]
    % [X51, 3, 3, 3]
    % [X52, 4, 3, 3]
    % [X53, 1, 4, 3]
    % [X54, 2, 4, 3]
    % [X55, 3, 4, 3]
    % [X56, 4, 4, 3]
    % [X57, 1, 5, 3]
    % [X58, 2, 5, 3]
    % [X59, 3, 5, 3]
    % [X60, 4, 5, 3]
    % [X61, 1, 1, 4]
    % [X62, 2, 1, 4]
    % [X63, 3, 1, 4]
    % [X64, 4, 1, 4]
    % [X65, 1, 2, 4]
    % [X66, 2, 2, 4]
    % [X67, 3, 2, 4]
    % [X68, 4, 2, 4]
    % [X69, 1, 3, 4]
    % [X70, 2, 3, 4]
    % [X71, 3, 3, 4]
    % [X72, 4, 3, 4]
    % [X73, 1, 4, 4]
    % [X74, 2, 4, 4]
    % [X75, 3, 4, 4]
    % [X76, 4, 4, 4]
    % [X77, 1, 5, 4]
    % [X78, 2, 5, 4]
    % [X79, 3, 5, 4]
    % [X80, 4, 5, 4]
    rho = min([rho,alpha1,alpha2,alpha3]);                        
    Y0 = ones(80,1);
    
    if(isequal(solveMethod,'ode45'))
        tspan = optSolver.tspan;
        options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
        [~,ES_t] = ode45(@(t,Y)  ODE_SEIISSICRSEIIIS_3(t,Y,beta1,gamma1,nu1,sigma1,eps1,...
                                         betaI2,betaC2,gamma2,sigma2,theta2,...
                                         beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,mu,b,rho),...
                                     tspan,Y0,options);
        ES = ES_t(end,:);
    end
    if(isequal(solveMethod,'fsolve'))
        options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
        [ES] = fsolve(@(Y)  ODE_SEIISSICRSEIIIS_3(1,Y,beta1,gamma1,nu1,sigma1,eps1,...
                                         betaI2,betaC2,gamma2,sigma2,theta2,...
                                         beta3,sigma3,tau3,nu3,gamma13,theta3,gamma33,mu,b,rho),...
                                      Y0,options);    
    end

    %P = max(sum(ES(end,[2:3,5:80]))/(b/mu),0);
    %PHIV = sum(ES(end,[2:4:78,3:4:79]))/(b/mu); %v2 v3
    %P123 = f*PHIV + sum(ES(end,[5:4:77,8:4:80]))/(b/mu); %v2 v3
    PHIV = sum(ES([2:4:58,3:4:59]))/(b/mu); %v3.2
    P123 = f*PHIV + sum(ES([5:4:57,8:4:60]))/(b/mu); %v3.2
    P = max(P123,0);
    %U123 = min(max(rho,0),max([alpha1,alpha2,alpha3]))*(P-c);
    %U123 = min(rho,min([alpha1,alpha2,alpha3]))*(P-c) ;
    U123 = rho*(P-c);
    dU  = 0;                                                            
 end