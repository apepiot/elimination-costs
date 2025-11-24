function [M] = ODESEIIS(S,I,J,lambda,epsilon, gamma,nu,sigma,b,rho,mu)
%     ODE.dS  = b - lambda*S + (1-epsilon)*(rho+nu)*I + (gamma+rho+nu)*J - mu*S;
%     ODE.dI = lambda*S - ((1-epsilon)*(rho+nu) + epsilon*sigma)*I - mu*I;
%     ODE.dJ = epsilon*sigma*I - (gamma+nu)*J - mu*J;
    
    M = [-lambda-mu, (1-epsilon)*(rho+nu), gamma+rho+nu;...
        lambda, -(epsilon*sigma+(1-epsilon)*(nu+rho))-mu, 0 ;...
        0, epsilon*sigma, -(gamma+rho+nu+mu)];
end