function [U12,dU,P] = U123_SICTPSEIIS2_v1(param1,param2,param3,mu,b,rho,c,f,solveMethod,optSolver)

    %SICTP
    theta1 = param1.theta; betaI1 = param1.betaI;   betaC1 = param1.betaC;
    sigma1 = param1.sigma; eta1 = param1.eta; zeta1 = param1.zeta;  p1 = param1.p;
    
    %SEIIS
    beta2 = param2.beta; gamma2 = param2.gamma; nu2 = param2.nu; 
    eps2 = param2.eps; sigma2 = param2.sigma;
    
    %SEIIS
    beta3 = param3.beta; gamma3 = param3.gamma; nu3 = param3.nu; 
    eps3 = param3.eps; sigma3 = param3.sigma;

    alpha1=param1.alpha;alpha2=param2.alpha;alpha3=param3.alpha; 
    
    rho = min([rho,alpha1,alpha2,alpha3]);
    Y0 = ones(7*4*4,1);
    if(isequal(solveMethod,'ode45'))
        tspan = optSolver.tspan;    
        options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
        [~,ES_t] = ode45(@(t,Y) ODE_SICTPSEIIS2(t,Y,betaI1,betaC1,sigma1,theta1,zeta1,eta1,p1,...
                            beta2,gamma2,nu2,eps2,sigma2,...
                            beta3,gamma3,nu3,eps3,sigma3,...
                            [],mu,b,rho),...
                         tspan,Y0, options);                  
         ES = ES_t(end,:);
    end
    if(isequal(solveMethod,'fsolve'))
        options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
        [ES] = fsolve(@(Y) ODE_SICTPSEIIS2(0,Y,betaI1,betaC1,sigma1,theta1,zeta1,eta1,p1,...
                            beta2,gamma2,nu2,eps2,sigma2,...
                            beta3,gamma3,nu3,eps3,sigma3,...
                            [],mu,b,rho),...
                         Y0,options);
    end

    
    %asymptomatic prevalence:
    PHIV = sum(ES([2,3,5,6,9,10,12,13,16,17,19,20,30,31,33,34,37,38,40,41,44,45,47,48,58,...
        59,61,62,65,66,68,69,72,73,75,76]))/(b/mu);             
    P12 = min(f*PHIV + sum(ES([8,11,14,15,18,21,29,32,35,36,39,42,43,46,49,57,60,63,64,67,70,71,74,77]))/(b/mu),1);
    P = max(P12,0);
    U12 = rho*(P-c);
    dU  = 0;
    
%     HIV1    STI1    STI2    no 
%     ____    ____    ____    ___
% 
%     "S"     "S"     "S"       1 
%     "I"     "S"     "S"       2 x
%     "C"     "S"     "S"       3 x
%     "P"     "S"     "S"       4
%     "Ip"    "S"     "S"       5 x
%     "Cp"    "S"     "S"       6 x
%     "T"     "S"     "S"       7
%     "S"     "E"     "S"       8 xx
%     "I"     "E"     "S"       9 x
%     "C"     "E"     "S"      10 x
%     "P"     "E"     "S"      11 xx
%     "Ip"    "E"     "S"      12 x
%     "Cp"    "E"     "S"      13 x
%     "T"     "E"     "S"      14 xx
%     "S"     "IA"    "S"      15 xx
%     "I"     "IA"    "S"      16 x
%     "C"     "IA"    "S"      17 x
%     "P"     "IA"    "S"      18 xx
%     "Ip"    "IA"    "S"      19 x
%     "Cp"    "IA"    "S"      20 x
%     "T"     "IA"    "S"      21 xx
%     "S"     "IS"    "S"      22
%     "I"     "IS"    "S"      23 
%     "C"     "IS"    "S"      24 
%     "P"     "IS"    "S"      25
%     "Ip"    "IS"    "S"      26 
%     "Cp"    "IS"    "S"      27 
%     "T"     "IS"    "S"      28
%     "S"     "S"     "E"      29 xx
%     "I"     "S"     "E"      30 x
%     "C"     "S"     "E"      31 x
%     "P"     "S"     "E"      32 xx
%     "Ip"    "S"     "E"      33 x
%     "Cp"    "S"     "E"      34 x
%     "T"     "S"     "E"      35 xx
%     "S"     "E"     "E"      36 xx
%     "I"     "E"     "E"      37 x
%     "C"     "E"     "E"      38 x
%     "P"     "E"     "E"      39 xx
%     "Ip"    "E"     "E"      40 x
%     "Cp"    "E"     "E"      41 x
%     "T"     "E"     "E"      42 xx
%     "S"     "IA"    "E"      43 xx
%     "I"     "IA"    "E"      44 x
%     "C"     "IA"    "E"      45 x
%     "P"     "IA"    "E"      46 xx
%     "Ip"    "IA"    "E"      47 x
%     "Cp"    "IA"    "E"      48 x
%     "T"     "IA"    "E"      49 xx
%     "S"     "IS"    "E"      50 
%     "I"     "IS"    "E"      51
%     "C"     "IS"    "E"      52
%     "P"     "IS"    "E"      53
%     "Ip"    "IS"    "E"      54
%     "Cp"    "IS"    "E"      55
%     "T"     "IS"    "E"      56
%     "S"     "S"     "IA"     57 xx
%     "I"     "S"     "IA"     58 x
%     "C"     "S"     "IA"     59 x
%     "P"     "S"     "IA"     60 xx
%     "Ip"    "S"     "IA"     61 x
%     "Cp"    "S"     "IA"     62 x
%     "T"     "S"     "IA"     63 xx
%     "S"     "E"     "IA"     64 xx
%     "I"     "E"     "IA"     65 x
%     "C"     "E"     "IA"     66 x
%     "P"     "E"     "IA"     67 xx
%     "Ip"    "E"     "IA"     68 x
%     "Cp"    "E"     "IA"     69 x
%     "T"     "E"     "IA"     70 xx
%     "S"     "IA"    "IA"     71 xx
%     "I"     "IA"    "IA"     72 x
%     "C"     "IA"    "IA"     73 x
%     "P"     "IA"    "IA"     74 xx
%     "Ip"    "IA"    "IA"     75 x
%     "Cp"    "IA"    "IA"     76 x
%     "T"     "IA"    "IA"     77 xx
%     "S"     "IS"    "IA"     78 
%     "I"     "IS"    "IA"     79 
%     "C"     "IS"    "IA"     80
%     "P"     "IS"    "IA"     81
%     "Ip"    "IS"    "IA"     82
%     "Cp"    "IS"    "IA"     83
%     "T"     "IS"    "IA"     84
%     "S"     "S"     "IS"     85
%     "I"     "S"     "IS"     86
%     "C"     "S"     "IS"     87
%     "P"     "S"     "IS"     88
%     "Ip"    "S"     "IS"     89
%     "Cp"    "S"     "IS"     90
%     "T"     "S"     "IS"     91
%     "S"     "E"     "IS"     92
%     "I"     "E"     "IS"     93
%     "C"     "E"     "IS"     94
%     "P"     "E"     "IS"     95
%     "Ip"    "E"     "IS"     96
%     "Cp"    "E"     "IS"     97
%     "T"     "E"     "IS"     98
%     "S"     "IA"    "IS"     99
%     "I"     "IA"    "IS"    100
%     "C"     "IA"    "IS"    101
%     "P"     "IA"    "IS"    102
%     "Ip"    "IA"    "IS"    103
%     "Cp"    "IA"    "IS"    104
%     "T"     "IA"    "IS"    105
%     "S"     "IS"    "IS"    106
%     "I"     "IS"    "IS"    107
%     "C"     "IS"    "IS"    108
%     "P"     "IS"    "IS"    109
%     "Ip"    "IS"    "IS"    110
%     "Cp"    "IS"    "IS"    111
%     "T"     "IS"    "IS"    112
end