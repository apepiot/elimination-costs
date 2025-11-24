function dPop=HCV(t, pop, beta, gamma,sigma,zeta,omega,b,mu,type)

    S=pop(1);
    I=pop(2);
    C=pop(3);
    A=pop(4);
    T=pop(5);

    dPop=zeros(5,1);
    
    if strcmp(type,'frequency')
        N = S+I+C+A+T;
    end
    if strcmp(type,'density')
        N = 1;
    end
    lambda = beta*(I+C)/N;
 
     dPop(1) = b -lambda*S - mu*S;%S
     dPop(2) = lambda*S - sigma*I - mu*I;%I
     dPop(3) = omega*sigma*I - gamma*C - mu*C;%C
     dPop(4) = (1-omega)*sigma*I + zeta*T - mu*A; %A
     dPop(5) = gamma*C - zeta*T - mu*T; %T
end %function
