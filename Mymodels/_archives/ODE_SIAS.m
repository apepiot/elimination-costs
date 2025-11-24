function [dY] = ODE_SIAS(t,Y,b,beta,gamma0,omega,ksi0,s,rho,mu,type)
%strategy 2 february 2021
    
    S=Y(1);I=Y(2);A=Y(3);
    
    if strcmp(type,'frequency')
        N = sum(Y);
    end
    if strcmp(type,'density')
        N = 1;
    end
    
    lambda = beta*(I+A)/N;
    gamma = gamma0+s*rho;
    ksi = ksi0+s*rho;

    dS      = b - lambda*S + gamma*I + ksi*A - mu*S; %S
    dI     = omega*lambda*S - (gamma + mu)*I; %I
    dA     = (1-omega)*lambda*S - (ksi + mu)*A; %A
    
    
    dY = [dS; dI ; dA];   
end

