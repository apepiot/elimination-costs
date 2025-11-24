function [U12,dU,P] = U123_SICTPSEIIISSEIIS_v1(param1,param2,param3,mu,b,rho,c,f,solveMethod,optSolver)

    %SICTP
    theta1 = param1.theta; betaI1 = param1.betaI;   betaC1 = param1.betaC;
    sigma1 = param1.sigma; eta1 = param1.eta; zeta1 = param1.zeta;  p1 = param1.p;
    
    %SEIIIS
    beta2 = param2.beta; sigma2 = param2.sigma; tau2 = param2.tau;
    theta2 = param2.theta; gamma23 = param2.gamma3; 
    
    %SEIIS
    beta3 = param3.beta; gamma3 = param3.gamma; nu3 = param3.nu; 
    eps3 = param3.eps; sigma3 = param3.sigma;

    alpha1=param1.alpha;alpha2=param2.alpha;  
    
    rho = min([rho,alpha1,alpha2]);
    Y0 = ones(7*5*4,1);
    if(isequal(solveMethod,'ode45'))
        tspan = optSolver.tspan;    
        options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
        [~,ES_t] = ode45(@(t,Y) ODE_SICTPSEIIISSEIIS(t,Y,betaI1,betaC1,sigma1,theta1,zeta1,eta1,p1,...
                            beta2,sigma2,gamma23,tau2,theta2,...
                            beta3,gamma3,nu3,eps3,sigma3,...
                            [],mu,b,rho),...
                         tspan,Y0, options);                  
         ES = ES_t(end,:);
    end
    if(isequal(solveMethod,'fsolve'))
        options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
        [ES] = fsolve(@(Y)  ODE_SICTPSEIIISSEIIS(0,Y,betaI1,betaC1,sigma1,theta1,zeta1,eta1,p1,...
                            beta2,sigma2,gamma23,tau2,theta2,...
                            beta3,gamma3,nu3,eps3,sigma3,...
                            [],mu,b,rho),...
                         Y0,options);
    end

    %%
    %asymptomatic prevalence:
    PHIV = sum(ES([2:7:100,3:7:101,5:7:103,6:7:104]))/(b/mu);             
    P12 = min(f*PHIV + sum(ES([8:7:99,11:7:102,14:7:105]))/(b/mu),1);
    P = max(P12,0);
    U12 = rho*(P-c);
    dU  = 0;
    
%     HIV1    syph1    STI1    no 
%     ____    _____    ____    ___
% 
%     "S"     "S"      "S"       1
%     "I"     "S"      "S"       2
%     "C"     "S"      "S"       3
%     "P"     "S"      "S"       4
%     "Ip"    "S"      "S"       5
%     "Cp"    "S"      "S"       6
%     "T"     "S"      "S"       7
%     "S"     "E"      "S"       8
%     "I"     "E"      "S"       9
%     "C"     "E"      "S"      10
%     "P"     "E"      "S"      11
%     "Ip"    "E"      "S"      12
%     "Cp"    "E"      "S"      13
%     "T"     "E"      "S"      14
%     "S"     "I1"     "S"      15
%     "I"     "I1"     "S"      16
%     "C"     "I1"     "S"      17
%     "P"     "I1"     "S"      18
%     "Ip"    "I1"     "S"      19
%     "Cp"    "I1"     "S"      20
%     "T"     "I1"     "S"      21
%     "S"     "I2"     "S"      22
%     "I"     "I2"     "S"      23
%     "C"     "I2"     "S"      24
%     "P"     "I2"     "S"      25
%     "Ip"    "I2"     "S"      26
%     "Cp"    "I2"     "S"      27
%     "T"     "I2"     "S"      28
%     "S"     "I3"     "S"      29
%     "I"     "I3"     "S"      30
%     "C"     "I3"     "S"      31
%     "P"     "I3"     "S"      32
%     "Ip"    "I3"     "S"      33
%     "Cp"    "I3"     "S"      34
%     "T"     "I3"     "S"      35
%     "S"     "S"      "E"      36
%     "I"     "S"      "E"      37
%     "C"     "S"      "E"      38
%     "P"     "S"      "E"      39
%     "Ip"    "S"      "E"      40
%     "Cp"    "S"      "E"      41
%     "T"     "S"      "E"      42
%     "S"     "E"      "E"      43
%     "I"     "E"      "E"      44
%     "C"     "E"      "E"      45
%     "P"     "E"      "E"      46
%     "Ip"    "E"      "E"      47
%     "Cp"    "E"      "E"      48
%     "T"     "E"      "E"      49
%     "S"     "I1"     "E"      50
%     "I"     "I1"     "E"      51
%     "C"     "I1"     "E"      52
%     "P"     "I1"     "E"      53
%     "Ip"    "I1"     "E"      54
%     "Cp"    "I1"     "E"      55
%     "T"     "I1"     "E"      56
%     "S"     "I2"     "E"      57
%     "I"     "I2"     "E"      58
%     "C"     "I2"     "E"      59
%     "P"     "I2"     "E"      60
%     "Ip"    "I2"     "E"      61
%     "Cp"    "I2"     "E"      62
%     "T"     "I2"     "E"      63
%     "S"     "I3"     "E"      64
%     "I"     "I3"     "E"      65
%     "C"     "I3"     "E"      66
%     "P"     "I3"     "E"      67
%     "Ip"    "I3"     "E"      68
%     "Cp"    "I3"     "E"      69
%     "T"     "I3"     "E"      70
%     "S"     "S"      "IA"     71
%     "I"     "S"      "IA"     72
%     "C"     "S"      "IA"     73
%     "P"     "S"      "IA"     74
%     "Ip"    "S"      "IA"     75
%     "Cp"    "S"      "IA"     76
%     "T"     "S"      "IA"     77
%     "S"     "E"      "IA"     78
%     "I"     "E"      "IA"     79
%     "C"     "E"      "IA"     80
%     "P"     "E"      "IA"     81
%     "Ip"    "E"      "IA"     82
%     "Cp"    "E"      "IA"     83
%     "T"     "E"      "IA"     84
%     "S"     "I1"     "IA"     85
%     "I"     "I1"     "IA"     86
%     "C"     "I1"     "IA"     87
%     "P"     "I1"     "IA"     88
%     "Ip"    "I1"     "IA"     89
%     "Cp"    "I1"     "IA"     90
%     "T"     "I1"     "IA"     91
%     "S"     "I2"     "IA"     92
%     "I"     "I2"     "IA"     93
%     "C"     "I2"     "IA"     94
%     "P"     "I2"     "IA"     95
%     "Ip"    "I2"     "IA"     96
%     "Cp"    "I2"     "IA"     97
%     "T"     "I2"     "IA"     98
%     "S"     "I3"     "IA"     99
%     "I"     "I3"     "IA"    100
%     "C"     "I3"     "IA"    101
%     "P"     "I3"     "IA"    102
%     "Ip"    "I3"     "IA"    103
%     "Cp"    "I3"     "IA"    104
%     "T"     "I3"     "IA"    105
%     "S"     "S"      "IS"    106
%     "I"     "S"      "IS"    107
%     "C"     "S"      "IS"    108
%     "P"     "S"      "IS"    109
%     "Ip"    "S"      "IS"    110
%     "Cp"    "S"      "IS"    111
%     "T"     "S"      "IS"    112
%     "S"     "E"      "IS"    113
%     "I"     "E"      "IS"    114
%     "C"     "E"      "IS"    115
%     "P"     "E"      "IS"    116
%     "Ip"    "E"      "IS"    117
%     "Cp"    "E"      "IS"    118
%     "T"     "E"      "IS"    119
%     "S"     "I1"     "IS"    120
%     "I"     "I1"     "IS"    121
%     "C"     "I1"     "IS"    122
%     "P"     "I1"     "IS"    123
%     "Ip"    "I1"     "IS"    124
%     "Cp"    "I1"     "IS"    125
%     "T"     "I1"     "IS"    126
%     "S"     "I2"     "IS"    127
%     "I"     "I2"     "IS"    128
%     "C"     "I2"     "IS"    129
%     "P"     "I2"     "IS"    130
%     "Ip"    "I2"     "IS"    131
%     "Cp"    "I2"     "IS"    132
%     "T"     "I2"     "IS"    133
%     "S"     "I3"     "IS"    134
%     "I"     "I3"     "IS"    135
%     "C"     "I3"     "IS"    136
%     "P"     "I3"     "IS"    137
%     "Ip"    "I3"     "IS"    138
%     "Cp"    "I3"     "IS"    139
%     "T"     "I3"     "IS"    140
end