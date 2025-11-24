function [dY] = ODE_SICR(t,Y,b, betaI, betaC, sigma, theta0, gamma0, mu,rho,type)
%strategy 2 february 2021
    
    S=Y(1);I=Y(2);C=Y(3);R=Y(4);
    
    theta = theta0+rho;
    gamma = gamma0+rho;
    
    if strcmp(type,'frequency')
        N = sum(Y);
    end
    if strcmp(type,'density')
        N = 1;
    end

    lambda = (betaI*I+betaC*C)/N;

    dS = b - lambda*S-mu*S;
    dI = lambda*S - sigma*I - gamma*I- mu*I;
    dC = sigma*I - theta*C - mu*C;
    dR = theta*C + gamma*I - mu*R;
    
    dY = [dS; dI; dC; dR];   
end

