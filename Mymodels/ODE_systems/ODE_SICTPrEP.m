function dY = ODE_SICTPrEP(t,Y,b,betaI,betaC,sigma,theta0,p,zeta,eta,mu,rho,type)
    %SICT model with PrEP (august 2023)
    %attention aux notations : 
    %theta0 (ici) = gamma0 (papier)
    %gamma0 (ici) = 0 (papier)    
    S=Y(1);I=Y(2);C=Y(3);P=Y(4);Ip=Y(5);Cp=Y(6);T=Y(7);   
    theta = theta0+rho;    
    if strcmp(type,'frequency')
        N = sum(Y);
    end
    if strcmp(type,'density')
        N = 1;
    end

    lambda = (betaI*(I+Ip)+betaC*(C+Cp))/N;

    dS  = (1-p)*b - lambda*S-mu*S;
    dI  = lambda*S - (sigma+rho+mu)*I;
    dC  = sigma*I - (theta+mu)*C;
    dP  = p*b - (1-zeta)*lambda*P - mu*P;
    dIp = (1-zeta)*lambda*P - (sigma+eta+mu)*Ip;
    dCp = sigma*Ip - (theta0+eta+mu)*Cp;
    dT  = rho*I+theta*C + eta*Ip + (theta0+eta)*Cp - mu*T;
        
    dY = [dS; dI; dC; dP; dIp; dCp; dT];
end

