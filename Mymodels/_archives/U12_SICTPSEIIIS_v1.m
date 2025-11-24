function [U12,dU,P] = U12_SICTPSEIIIS_v1(param1,param2,mu,b,rho,c,f,solveMethod,optSolver)

    %SICTP
    theta1 = param1.theta; betaI1 = param1.betaI;   betaC1 = param1.betaC;
    sigma1 = param1.sigma; eta1 = param1.eta; zeta1 = param1.zeta;  p1 = param1.p;
    
    %SEIIIS
    beta2=param2.beta;sigma2=param2.sigma;tau2=param2.tau;
    theta2=param2.theta;gamma32=param2.gamma3;
    
    alpha1=param1.alpha;alpha2=param2.alpha;  
    
    rho = min([rho,alpha1,alpha2]);
    Y0 = ones(7*5,1);
    if(isequal(solveMethod,'ode45'))
        tspan = optSolver.tspan;    
        options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
        [~,ES_t] = ode45(@(t,Y) ODE_SICTPSEIIIS(t,Y,betaI1,betaC1,sigma1,theta1,zeta1,eta1,p1,...
                         beta2,sigma2,gamma32,tau2,theta2, ...
                         [], mu,b,rho),...
                         tspan,Y0, options);
         ES = ES_t(end,:);
    end
    if(isequal(solveMethod,'fsolve'))
        options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
        [ES] = fsolve(@(Y)  ODE_SICTPSEIIIS(0,Y,betaI1,betaC1,sigma1,theta1,zeta1,eta1,p1,...
                         beta2,sigma2,gamma32,tau2,theta2, ...
                         [], mu,b,rho),...
                         Y0,options);
    end

    PHIV = sum(ES([2:7:30,3:7:31,5:7:33,6:7:34]))/(b/mu);             
    P12 = min(f*PHIV + sum(ES([8:7:29,11:7:32,14:7:35]))/(b/mu),1);
    P = max(P12,0);
    U12 = rho*(P-c);
    dU  = 0;
    
%     HIV1    syph1    no
%     ____    _____    __
% 
%     "S"     "S"       1
%     "I"     "S"       2 x
%     "C"     "S"       3 x
%     "P"     "S"       4
%     "Ip"    "S"       5 x
%     "Cp"    "S"       6 x
%     "T"     "S"       7
%     "S"     "E"       8 x
%     "I"     "E"       9 x
%     "C"     "E"      10 x
%     "P"     "E"      11 x
%     "Ip"    "E"      12 x
%     "Cp"    "E"      13 x
%     "T"     "E"      14 x
%     "S"     "I1"     15 x
%     "I"     "I1"     16 x
%     "C"     "I1"     17 x
%     "P"     "I1"     18 x
%     "Ip"    "I1"     19 x
%     "Cp"    "I1"     20 x
%     "T"     "I1"     21 x
%     "S"     "I2"     22 x
%     "I"     "I2"     23 x
%     "C"     "I2"     24 x
%     "P"     "I2"     25 x
%     "Ip"    "I2"     26 x
%     "Cp"    "I2"     27 x
%     "T"     "I2"     28 x
%     "S"     "I3"     29 x
%     "I"     "I3"     30 x
%     "C"     "I3"     31 x
%     "P"     "I3"     32 x
%     "Ip"    "I3"     33 x
%     "Cp"    "I3"     34 x
%     "T"     "I3"     35 x
end