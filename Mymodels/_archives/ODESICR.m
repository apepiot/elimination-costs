function [M] = ODESICR(S,I,J,lambda, theta, gamma,sigma,b,rho,mu)
%       dS = b - lambda*S-mu*S;
%       dI = lambda*S - (sigma+gamma)*I - mu*I;
%       dC = sigma*I - theta*C - mu*C;
%       dR = theta*C + gamma*I - mu*R;
    gamma = rho;
    thetap = theta+rho;
    M = [-lambda-mu, 0, 0, 0;...
        lambda, - (sigma+gamma+mu), 0,0 ;...
        0, sigma, - thetap-mu,0;...
        0,gamma,thetap,-mu];
end