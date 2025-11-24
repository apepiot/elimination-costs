function [U123,dU,P] = U123_SEIIS2SICR_v4(param1,param2,param3,...
                                     mu,b,vecRho,c,f,solveMethod,optSolver)
    %     HIV,STI,STI                             
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
    % [X17, 1, 1, 2]
    % [X18, 2, 1, 2]
    % [X19, 3, 1, 2]
    % [X20, 4, 1, 2]
    % [X21, 1, 2, 2]
    % [X22, 2, 2, 2]
    % [X23, 3, 2, 2]
    % [X24, 4, 2, 2]
    % [X25, 1, 3, 2]
    % [X26, 2, 3, 2]
    % [X27, 3, 3, 2]
    % [X28, 4, 3, 2]
    % [X29, 1, 4, 2]
    % [X30, 2, 4, 2]
    % [X31, 3, 4, 2]
    % [X32, 4, 4, 2]
    % [X33, 1, 1, 3]
    % [X34, 2, 1, 3]
    % [X35, 3, 1, 3]
    % [X36, 4, 1, 3]
    % [X37, 1, 2, 3]
    % [X38, 2, 2, 3]
    % [X39, 3, 2, 3]
    % [X40, 4, 2, 3]
    % [X41, 1, 3, 3]
    % [X42, 2, 3, 3]
    % [X43, 3, 3, 3]
    % [X44, 4, 3, 3]
    % [X45, 1, 4, 3]
    % [X46, 2, 4, 3]
    % [X47, 3, 4, 3]
    % [X48, 4, 4, 3]
    % [X49, 1, 1, 4]
    % [X50, 2, 1, 4]
    % [X51, 3, 1, 4]
    % [X52, 4, 1, 4]
    % [X53, 1, 2, 4]
    % [X54, 2, 2, 4]
    % [X55, 3, 2, 4]
    % [X56, 4, 2, 4]
    % [X57, 1, 3, 4]
    % [X58, 2, 3, 4]
    % [X59, 3, 3, 4]
    % [X60, 4, 3, 4]
    % [X61, 1, 4, 4]
    % [X62, 2, 4, 4]
    % [X63, 3, 4, 4]
    % [X64, 4, 4, 4]
    beta1=param1.beta;gamma1=param1.gamma;nu1=param1.nu;sigma1=param1.sigma;eps1=param1.eps;
    beta2=param2.beta;gamma2=param2.gamma;nu2=param2.nu;sigma2=param2.sigma;eps2=param2.eps;
    betaI3=param3.betaI;betaC3=param3.betaC;gamma3=param3.gamma;sigma3=param3.sigma;theta3=param3.theta;
    alpha1=param1.alpha;  alpha2=param2.alpha; alpha3=param3.alpha;
    
    P=zeros(1,length(vecRho));
    U123=zeros(1,length(vecRho));i=1;
    for rhop=vecRho
    rho = rhop;%min([rhop,alpha1,alpha2,alpha3]);
    Y0 = ones(64,1);
    if(isequal(solveMethod,'ode45'))
        tspan = optSolver.tspan;
        options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
        [~,ES_t] = ode45(@(t,Y)  ODE_SEIIS2SICR_v4(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                     beta2,nu2,eps2,sigma2,gamma2,...
                                     betaI3,betaC3,gamma3,sigma3,theta3,mu,b,rho),...
                                     tspan,Y0,options);
        ES = ES_t(end,:);
        %P = max(sum(ES(end,[2:3,5:64]))/(b/mu),0);
        %PHIV = sum(ES(end,[2:4:62,3:4:63]))/(b/mu); %v2
        %P123 = f*PHIV+sum(ES(end,[5:4:61,8:4:64]))/(b/mu); %v2
    end
    if(isequal(solveMethod,'fsolve'))
        options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
        [ES] = fsolve(@(Y)  ODE_SEIIS2SICR_v4(1,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                      beta2,nu2,eps2,sigma2,gamma2,...
                                      betaI3,betaC3,gamma3,sigma3,theta3,mu,b,rho),...
                                      Y0,options);    
    end
    %version below : without IS for STIs
    PHIV = sum(ES([2:4:10,18:4:26,34:4:42, 3:4:11,19:4:27,35:4:43]))/(b/mu); %v2.2
    P123 = min(f*PHIV + sum(ES([5:4:9,17:4:25,33:4:41, 8:4:12,20:4:28,36:4:44]))/(b/mu),1); %v2
    
    %version below : without Exposed for STIs
    %PHIV = sum(ES([2:4:62,3:4:63]))/(b/mu); %v2.2
    %P123 = f*PHIV+sum(ES([9,13,25:4:61,12,16,28:4:64]))/(b/mu);
    
    P(i) = max(P123,0);
    %U12 = min(max(rho,0),max([alpha1,alpha2,alpha3]))*(P-c);
    %U123 = min(rho,min([alpha1,alpha2,alpha3]))*(P-c) ;
    U123(i) = rho*(P(i)-c);
    dU  = 0;
    i=i+1;
    end
end