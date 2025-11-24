function [U12,dU,P] = U12_SEIIS2_v6(param1,param2,mu,b,vecRho,c,f,solveMethod,optSolver)
beta1=param1.beta;gamma1=param1.gamma;nu1=param1.nu;sigma1=param1.sigma;eps1=param1.eps;
beta2=param2.beta;gamma2=param2.gamma;nu2=param2.nu;sigma2=param2.sigma;eps2=param2.eps;

U12 = zeros(1,length(vecRho));
dU  = zeros(1,length(vecRho));
P   = zeros(1,length(vecRho));
i=1;
for rho=vecRho
    Y0 = ones(16,1);
    if (isequal(solveMethod,'ode45'))
        tspan = optSolver.tspan;
        options = odeset('RelTol',1e-3,'Stats','off');%,'OutputFcn',@odeplot);
        [~,ES_t] = ode45(@(t,Y) ODE_SEIIS2_v6(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
            beta2,nu2,eps2,sigma2,gamma2,mu,b,rho),tspan,Y0,options);
        ES = ES_t(end,:);
        
        if min(ES)<5
            warning('probleme de convergence de ES')
        end
    elseif (isequal(solveMethod,'fsolve'))
        options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
        [ES] = fsolve(@(Y)  ODE_SEIIS2_v6(1,Y,beta1,nu1,eps1,sigma1,gamma1,...
            beta2,nu2,eps2,sigma2,gamma2,mu,b,rho),...
            Y0,options);
        ES = max(ES,0);
    elseif (isequal(solveMethod,'knitro-ampl'))
        % to do
    end
    P(i)   = P_kit(ES,{'Ct','Ng'});
    U12(i) = rho*(P(i)-c);
    dU(i)  = 0;
    i=i+1;
end
end