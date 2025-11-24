function [U12,dU,P] = U12_SEIISSICR(beta1,gamma1,nu1,sigma1,eps1,betaI2,betaC2,gamma2,sigma2,theta2,mu,b,vecRho,c,alpha1,alpha2,f,solveMethod,optSolver)
    %     HIV,STI
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
    i=1;U12=zeros(length(vecRho),1);dU=zeros(length(vecRho),1);
    for RHO=vecRho
        rho = min([RHO,alpha1,alpha2]);Y0 = ones(4*4,1);
        if(isequal(solveMethod,'ode45'))
            tspan = optSolver.tspan;
            options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
            [~,Ys] = ode45(@(t,Y) ODE_SEIISSICR_3(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                         betaI2,betaC2,gamma2,sigma2,theta2,mu,b,rho),tspan,Y0,options);
            ES=Ys(end,:);  
        end
        if(isequal(solveMethod,'fsolve'))
            options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
            [ES] = fsolve(@(Y)  ODE_SEIISSICR_3(1,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                         betaI2,betaC2,gamma2,sigma2,theta2,mu,b,rho),...
                                      Y0,options);
        end
        %f = 10;
        %disp(['attention: there is a factor ',num2str(f), ' before prevalence of HIV in U12_SEIISSICR.m']);
        %on retire les cas où un individu est symptomatic (4) pour la STI
        P12 = f*sum(ES([2:3,6:7,10:11])) + sum(ES([5,8,9,12]));
        P = max(P12,0)/(b/mu);%*(rho<max([alpha1,alpha2]));
        %P = max(sum(ES(end,[2:3,5:16]))/(b/mu),0);
        %U12 = rho*(P-c); %min(max(rho,0),max([alpha1,alpha2]))*(P-c);
        %U12 = min(rho,min([alpha1,alpha2]))*(P-c) ;
        U12(i) = rho*(P-c)*(rho>=0);
        dU(i)  = 0;
        i=i+1;
    end
end