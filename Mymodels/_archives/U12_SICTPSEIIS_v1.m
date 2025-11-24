function [U12,dU,P] = U12_SICTPSEIIS_v1(param1,param2,mu,b,vecRho,c,f,solveMethod,optSolver)

    %SICTP
    theta1 = param1.theta; betaI1 = param1.betaI;   betaC1 = param1.betaC;
    sigma1 = param1.sigma; eta1 = param1.eta; zeta1 = param1.zeta;  p1 = param1.p;
    alpha1 = param1.alpha;
    
    %SEIIS
    beta2 = param2.beta; gamma2 = param2.gamma; nu2 = param2.nu; sigma2 = param2.sigma; eps2 = param2.eps;
    alpha2 = param2.alpha;
    
    i=1; U12=zeros(length(vecRho),1); dU=zeros(length(vecRho),1);
    for RHO=vecRho
        rho = min([RHO,alpha1,alpha2]); Y0 = ones(4*7,1);
        if(isequal(solveMethod,'ode45'))
            tspan = optSolver.tspan;
            options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
            [~,Ys] = ode45(@(t,Y) ODE_SICTPSEIIS(t,Y,betaI1,betaC1,sigma1,theta1,zeta1,eta1,p1,...
                                         beta2,gamma2,nu2,eps2,sigma2,...
                                         [],mu,b,rho),...
                                        tspan,Y0,options);                     
            ES=Ys(end,:);  
        end
        if(isequal(solveMethod,'fsolve'))
            options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
            [ES] = fsolve(@(Y)  ODE_SICTPSEIIS(0,Y,betaI1,betaC1,sigma1,theta1,zeta1,eta1,p1,...
                                         beta2,gamma2,nu2,eps2,sigma2,...
                                         [],mu,b,rho),...
                                      Y0,options);
        end
        
        %prevalence of asymptomatic and untreated:
        PHIV = f*sum(ES([2,3,5,6,9,10,12,13,16,17,19,20]));
        P12 = PHIV+sum(ES([8,11,14,15,18,21])); 
        % table2array(tabComp((tabComp.HIV1=="I" | tabComp.HIV1=="Ip" | ...
        % tabComp.HIV1=="C" | tabComp.HIV1=="Cp" | tabComp.STI1=="E" | ...
        % tabComp.STI1=="IA") & tabComp.STI1~="IS","no"))'
        P = min(max(P12,0)/(b/mu),1);
        U12(i) = rho*(P-c)*(rho>=0);
        dU(i)  = 0;
        i=i+1;
    end
    
%     HIV1    STI1    no
%     ____    ____    __
% 
%     "S"     "S"      1
%     "I"     "S"      2
%     "C"     "S"      3
%     "P"     "S"      4
%     "Ip"    "S"      5
%     "Cp"    "S"      6
%     "T"     "S"      7
%     "S"     "E"      8
%     "I"     "E"      9
%     "C"     "E"     10
%     "P"     "E"     11
%     "Ip"    "E"     12
%     "Cp"    "E"     13
%     "T"     "E"     14
%     "S"     "IA"    15
%     "I"     "IA"    16
%     "C"     "IA"    17
%     "P"     "IA"    18
%     "Ip"    "IA"    19
%     "Cp"    "IA"    20
%     "T"     "IA"    21
%     "S"     "IS"    22
%     "I"     "IS"    23
%     "C"     "IS"    24
%     "P"     "IS"    25
%     "Ip"    "IS"    26
%     "Cp"    "IS"    27
%     "T"     "IS"    28
end