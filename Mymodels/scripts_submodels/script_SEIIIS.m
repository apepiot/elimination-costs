%SEI1I2I3S (syphilis)
%% numerically
clear all
[beta,sigma,gamma10,gamma30,theta,nu,b,mu,rho] = random_parameters(true, true);
tau=1;
Y0 = [1,1,1,1,1]; tspan = 0:1:1000;
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) SEIIIS(t,Y,b,beta, sigma, gamma10, gamma30, tau, theta, nu, mu,rho),tspan,Y0, options);
T = Ys(end,:)
Nnum = sum(T) %b/mu

%% analytically
clear all
syms b lambda beta sigma gamma10 gamma30 tau theta nu mu S E I1 I2 I3 rho

N = b/mu;
lambda = beta*(I1+I2)/N;

dS  = b - lambda*S - mu*S + rho*E +(gamma10+rho)*I1+rho*I2+(nu+gamma30+rho)*I3;
dE  = lambda*S - (sigma+rho+mu)*E;
dI1 = sigma*E - (gamma10+rho+tau+mu)*I1;
dI2 = tau*I1 - (theta+rho+mu)*I2;
dI3 = theta*I2 - (nu+gamma30+rho+mu)*I3;

eqn1 = dS==0;
eqn2 = dE==0;
eqn3 = dI1==0;
eqn4 = dI2==0;
eqn5 = dI3==0;

sol = solve([eqn1,eqn2,eqn3,eqn4,eqn5], [S,E,I1,I2,I3])
%S
Smat = (b*(gamma10*mu^2 + gamma10*rho^2 + 3*mu*rho^2 + 3*mu^2*rho + mu^2*sigma + mu^2*tau + mu^2*theta + rho^2*sigma + rho^2*tau + rho^2*theta + mu^3 + rho^3 + 2*gamma10*mu*rho + gamma10*mu*sigma + gamma10*mu*theta + gamma10*rho*sigma + gamma10*rho*theta + gamma10*sigma*theta + 2*mu*rho*sigma + 2*mu*rho*tau + 2*mu*rho*theta + mu*sigma*tau + mu*sigma*theta + mu*tau*theta + rho*sigma*tau + rho*sigma*theta + rho*tau*theta + sigma*tau*theta))/(mu*sigma*(beta*mu + beta*rho + beta*tau + beta*theta))
%E
Emat = -((mu + rho + theta)*(gamma30 + mu + nu + rho)*(gamma10 + mu + rho + tau)*(b*mu^3 + b*rho^3 + 3*b*mu*rho^2 + 3*b*mu^2*rho + b*mu^2*sigma + b*mu^2*tau + b*mu^2*theta + b*rho^2*sigma + b*rho^2*tau + b*rho^2*theta + b*gamma10*mu^2 + b*gamma10*rho^2 - b*beta*mu*sigma + 2*b*gamma10*mu*rho - b*beta*rho*sigma + b*gamma10*mu*sigma + b*gamma10*mu*theta - b*beta*sigma*tau - b*beta*sigma*theta + b*gamma10*rho*sigma + b*gamma10*rho*theta + b*gamma10*sigma*theta + 2*b*mu*rho*sigma + 2*b*mu*rho*tau + 2*b*mu*rho*theta + b*mu*sigma*tau + b*mu*sigma*theta + b*mu*tau*theta + b*rho*sigma*tau + b*rho*sigma*theta + b*rho*tau*theta + b*sigma*tau*theta))/(sigma*(beta*mu^5 + beta*mu*rho^4 + 4*beta*mu^4*rho + beta*mu^4*sigma + 2*beta*mu^4*tau + 2*beta*mu^4*theta + 4*beta*mu^2*rho^3 + 6*beta*mu^3*rho^2 + beta*mu^3*tau^2 + beta*mu^3*theta^2 + beta*gamma10*mu^4 + beta*gamma30*mu^4 + beta*mu^4*nu + beta*gamma10*gamma30*mu^3 + beta*gamma10*mu^3*nu + beta*gamma10*mu*rho^3 + 3*beta*gamma10*mu^3*rho + beta*gamma30*mu*rho^3 + 3*beta*gamma30*mu^3*rho + beta*gamma30*mu^3*sigma + beta*gamma10*mu^3*tau + 2*beta*gamma30*mu^3*tau + 2*beta*gamma10*mu^3*theta + 2*beta*gamma30*mu^3*theta + beta*mu*nu*rho^3 + 3*beta*mu^3*nu*rho + beta*mu^3*nu*sigma + 2*beta*mu^3*nu*tau + 2*beta*mu^3*nu*theta + beta*mu*rho^3*sigma + 3*beta*mu^3*rho*sigma + 2*beta*mu*rho^3*tau + 6*beta*mu^3*rho*tau + 2*beta*mu*rho^3*theta + 6*beta*mu^3*rho*theta + 2*beta*mu^3*sigma*tau + 2*beta*mu^3*sigma*theta + 3*beta*mu^3*tau*theta + 3*beta*gamma10*mu^2*rho^2 + 3*beta*gamma30*mu^2*rho^2 + beta*gamma30*mu^2*tau^2 + beta*gamma10*mu^2*theta^2 + beta*gamma30*mu^2*theta^2 + 3*beta*mu^2*nu*rho^2 + beta*mu^2*nu*tau^2 + beta*mu^2*nu*theta^2 + 3*beta*mu^2*rho^2*sigma + beta*mu*rho^2*tau^2 + 2*beta*mu^2*rho*tau^2 + 6*beta*mu^2*rho^2*tau + beta*mu*rho^2*theta^2 + 2*beta*mu^2*rho*theta^2 + 6*beta*mu^2*rho^2*theta + beta*mu^2*sigma*tau^2 + beta*mu^2*sigma*theta^2 + beta*mu^2*tau*theta^2 + beta*mu^2*tau^2*theta + beta*gamma10*gamma30*mu*rho^2 + 2*beta*gamma10*gamma30*mu^2*rho + beta*gamma10*gamma30*mu^2*tau + beta*gamma10*gamma30*mu*theta^2 + 2*beta*gamma10*gamma30*mu^2*theta + beta*gamma10*mu*nu*rho^2 + 2*beta*gamma10*mu^2*nu*rho + beta*gamma10*mu^2*nu*tau + beta*gamma10*mu*nu*theta^2 + 2*beta*gamma10*mu^2*nu*theta + beta*gamma30*mu*rho^2*sigma + 2*beta*gamma30*mu^2*rho*sigma + beta*gamma10*mu*rho^2*tau + 2*beta*gamma10*mu^2*rho*tau + beta*gamma30*mu*rho*tau^2 + 2*beta*gamma30*mu*rho^2*tau + 4*beta*gamma30*mu^2*rho*tau + beta*gamma10*mu*rho*theta^2 + 2*beta*gamma10*mu*rho^2*theta + 4*beta*gamma10*mu^2*rho*theta + beta*gamma30*mu*rho*theta^2 + 2*beta*gamma30*mu*rho^2*theta + 4*beta*gamma30*mu^2*rho*theta + beta*gamma30*mu*sigma*tau^2 + 2*beta*gamma30*mu^2*sigma*tau + beta*gamma30*mu*sigma*theta^2 + 2*beta*gamma30*mu^2*sigma*theta + beta*gamma10*mu^2*tau*theta + beta*gamma30*mu*tau*theta^2 + beta*gamma30*mu*tau^2*theta + 3*beta*gamma30*mu^2*tau*theta + beta*mu*nu*rho^2*sigma + 2*beta*mu^2*nu*rho*sigma + beta*mu*nu*rho*tau^2 + 2*beta*mu*nu*rho^2*tau + 4*beta*mu^2*nu*rho*tau + beta*mu*nu*rho*theta^2 + 2*beta*mu*nu*rho^2*theta + 4*beta*mu^2*nu*rho*theta + beta*mu*nu*sigma*tau^2 + 2*beta*mu^2*nu*sigma*tau + beta*mu*nu*sigma*theta^2 + 2*beta*mu^2*nu*sigma*theta + beta*mu*nu*tau*theta^2 + beta*mu*nu*tau^2*theta + 3*beta*mu^2*nu*tau*theta + beta*mu*rho*sigma*tau^2 + 2*beta*mu*rho^2*sigma*tau + 4*beta*mu^2*rho*sigma*tau + beta*mu*rho*sigma*theta^2 + 2*beta*mu*rho^2*sigma*theta + 4*beta*mu^2*rho*sigma*theta + beta*mu*rho*tau*theta^2 + beta*mu*rho*tau^2*theta + 3*beta*mu*rho^2*tau*theta + 6*beta*mu^2*rho*tau*theta + beta*mu*sigma*tau*theta^2 + beta*mu*sigma*tau^2*theta + 3*beta*mu^2*sigma*tau*theta + beta*gamma10*gamma30*mu*rho*tau + 2*beta*gamma10*gamma30*mu*rho*theta + beta*gamma10*gamma30*mu*tau*theta + beta*gamma10*mu*nu*rho*tau + 2*beta*gamma10*mu*nu*rho*theta + beta*gamma10*mu*nu*tau*theta + 2*beta*gamma30*mu*rho*sigma*tau + 2*beta*gamma30*mu*rho*sigma*theta + beta*gamma10*mu*rho*tau*theta + 3*beta*gamma30*mu*rho*tau*theta + 2*beta*gamma30*mu*sigma*tau*theta + 2*beta*mu*nu*rho*sigma*tau + 2*beta*mu*nu*rho*sigma*theta + 3*beta*mu*nu*rho*tau*theta + 2*beta*mu*nu*sigma*tau*theta + 3*beta*mu*rho*sigma*tau*theta))
%I1
I1mat= -((mu + rho + theta)*(gamma30 + mu + nu + rho)*(b*mu^3 + b*rho^3 + 3*b*mu*rho^2 + 3*b*mu^2*rho + b*mu^2*sigma + b*mu^2*tau + b*mu^2*theta + b*rho^2*sigma + b*rho^2*tau + b*rho^2*theta + b*gamma10*mu^2 + b*gamma10*rho^2 - b*beta*mu*sigma + 2*b*gamma10*mu*rho - b*beta*rho*sigma + b*gamma10*mu*sigma + b*gamma10*mu*theta - b*beta*sigma*tau - b*beta*sigma*theta + b*gamma10*rho*sigma + b*gamma10*rho*theta + b*gamma10*sigma*theta + 2*b*mu*rho*sigma + 2*b*mu*rho*tau + 2*b*mu*rho*theta + b*mu*sigma*tau + b*mu*sigma*theta + b*mu*tau*theta + b*rho*sigma*tau + b*rho*sigma*theta + b*rho*tau*theta + b*sigma*tau*theta))/(beta*mu^5 + beta*mu*rho^4 + 4*beta*mu^4*rho + beta*mu^4*sigma + 2*beta*mu^4*tau + 2*beta*mu^4*theta + 4*beta*mu^2*rho^3 + 6*beta*mu^3*rho^2 + beta*mu^3*tau^2 + beta*mu^3*theta^2 + beta*gamma10*mu^4 + beta*gamma30*mu^4 + beta*mu^4*nu + beta*gamma10*gamma30*mu^3 + beta*gamma10*mu^3*nu + beta*gamma10*mu*rho^3 + 3*beta*gamma10*mu^3*rho + beta*gamma30*mu*rho^3 + 3*beta*gamma30*mu^3*rho + beta*gamma30*mu^3*sigma + beta*gamma10*mu^3*tau + 2*beta*gamma30*mu^3*tau + 2*beta*gamma10*mu^3*theta + 2*beta*gamma30*mu^3*theta + beta*mu*nu*rho^3 + 3*beta*mu^3*nu*rho + beta*mu^3*nu*sigma + 2*beta*mu^3*nu*tau + 2*beta*mu^3*nu*theta + beta*mu*rho^3*sigma + 3*beta*mu^3*rho*sigma + 2*beta*mu*rho^3*tau + 6*beta*mu^3*rho*tau + 2*beta*mu*rho^3*theta + 6*beta*mu^3*rho*theta + 2*beta*mu^3*sigma*tau + 2*beta*mu^3*sigma*theta + 3*beta*mu^3*tau*theta + 3*beta*gamma10*mu^2*rho^2 + 3*beta*gamma30*mu^2*rho^2 + beta*gamma30*mu^2*tau^2 + beta*gamma10*mu^2*theta^2 + beta*gamma30*mu^2*theta^2 + 3*beta*mu^2*nu*rho^2 + beta*mu^2*nu*tau^2 + beta*mu^2*nu*theta^2 + 3*beta*mu^2*rho^2*sigma + beta*mu*rho^2*tau^2 + 2*beta*mu^2*rho*tau^2 + 6*beta*mu^2*rho^2*tau + beta*mu*rho^2*theta^2 + 2*beta*mu^2*rho*theta^2 + 6*beta*mu^2*rho^2*theta + beta*mu^2*sigma*tau^2 + beta*mu^2*sigma*theta^2 + beta*mu^2*tau*theta^2 + beta*mu^2*tau^2*theta + beta*gamma10*gamma30*mu*rho^2 + 2*beta*gamma10*gamma30*mu^2*rho + beta*gamma10*gamma30*mu^2*tau + beta*gamma10*gamma30*mu*theta^2 + 2*beta*gamma10*gamma30*mu^2*theta + beta*gamma10*mu*nu*rho^2 + 2*beta*gamma10*mu^2*nu*rho + beta*gamma10*mu^2*nu*tau + beta*gamma10*mu*nu*theta^2 + 2*beta*gamma10*mu^2*nu*theta + beta*gamma30*mu*rho^2*sigma + 2*beta*gamma30*mu^2*rho*sigma + beta*gamma10*mu*rho^2*tau + 2*beta*gamma10*mu^2*rho*tau + beta*gamma30*mu*rho*tau^2 + 2*beta*gamma30*mu*rho^2*tau + 4*beta*gamma30*mu^2*rho*tau + beta*gamma10*mu*rho*theta^2 + 2*beta*gamma10*mu*rho^2*theta + 4*beta*gamma10*mu^2*rho*theta + beta*gamma30*mu*rho*theta^2 + 2*beta*gamma30*mu*rho^2*theta + 4*beta*gamma30*mu^2*rho*theta + beta*gamma30*mu*sigma*tau^2 + 2*beta*gamma30*mu^2*sigma*tau + beta*gamma30*mu*sigma*theta^2 + 2*beta*gamma30*mu^2*sigma*theta + beta*gamma10*mu^2*tau*theta + beta*gamma30*mu*tau*theta^2 + beta*gamma30*mu*tau^2*theta + 3*beta*gamma30*mu^2*tau*theta + beta*mu*nu*rho^2*sigma + 2*beta*mu^2*nu*rho*sigma + beta*mu*nu*rho*tau^2 + 2*beta*mu*nu*rho^2*tau + 4*beta*mu^2*nu*rho*tau + beta*mu*nu*rho*theta^2 + 2*beta*mu*nu*rho^2*theta + 4*beta*mu^2*nu*rho*theta + beta*mu*nu*sigma*tau^2 + 2*beta*mu^2*nu*sigma*tau + beta*mu*nu*sigma*theta^2 + 2*beta*mu^2*nu*sigma*theta + beta*mu*nu*tau*theta^2 + beta*mu*nu*tau^2*theta + 3*beta*mu^2*nu*tau*theta + beta*mu*rho*sigma*tau^2 + 2*beta*mu*rho^2*sigma*tau + 4*beta*mu^2*rho*sigma*tau + beta*mu*rho*sigma*theta^2 + 2*beta*mu*rho^2*sigma*theta + 4*beta*mu^2*rho*sigma*theta + beta*mu*rho*tau*theta^2 + beta*mu*rho*tau^2*theta + 3*beta*mu*rho^2*tau*theta + 6*beta*mu^2*rho*tau*theta + beta*mu*sigma*tau*theta^2 + beta*mu*sigma*tau^2*theta + 3*beta*mu^2*sigma*tau*theta + beta*gamma10*gamma30*mu*rho*tau + 2*beta*gamma10*gamma30*mu*rho*theta + beta*gamma10*gamma30*mu*tau*theta + beta*gamma10*mu*nu*rho*tau + 2*beta*gamma10*mu*nu*rho*theta + beta*gamma10*mu*nu*tau*theta + 2*beta*gamma30*mu*rho*sigma*tau + 2*beta*gamma30*mu*rho*sigma*theta + beta*gamma10*mu*rho*tau*theta + 3*beta*gamma30*mu*rho*tau*theta + 2*beta*gamma30*mu*sigma*tau*theta + 2*beta*mu*nu*rho*sigma*tau + 2*beta*mu*nu*rho*sigma*theta + 3*beta*mu*nu*rho*tau*theta + 2*beta*mu*nu*sigma*tau*theta + 3*beta*mu*rho*sigma*tau*theta)
%I2
I2mat = -((gamma30 + mu + nu + rho)*(b*mu^3*tau + b*rho^3*tau + b*mu^2*tau^2 + b*rho^2*tau^2 + b*gamma10*mu^2*tau - b*beta*sigma*tau^2 + b*gamma10*rho^2*tau + 2*b*mu*rho*tau^2 + 3*b*mu*rho^2*tau + 3*b*mu^2*rho*tau + b*mu*sigma*tau^2 + b*mu^2*sigma*tau + b*mu*tau^2*theta + b*mu^2*tau*theta + b*rho*sigma*tau^2 + b*rho^2*sigma*tau + b*rho*tau^2*theta + b*rho^2*tau*theta + b*sigma*tau^2*theta - b*beta*mu*sigma*tau + 2*b*gamma10*mu*rho*tau - b*beta*rho*sigma*tau + b*gamma10*mu*sigma*tau + b*gamma10*mu*tau*theta - b*beta*sigma*tau*theta + b*gamma10*rho*sigma*tau + b*gamma10*rho*tau*theta + b*gamma10*sigma*tau*theta + 2*b*mu*rho*sigma*tau + 2*b*mu*rho*tau*theta + b*mu*sigma*tau*theta + b*rho*sigma*tau*theta))/(beta*mu^5 + beta*mu*rho^4 + 4*beta*mu^4*rho + beta*mu^4*sigma + 2*beta*mu^4*tau + 2*beta*mu^4*theta + 4*beta*mu^2*rho^3 + 6*beta*mu^3*rho^2 + beta*mu^3*tau^2 + beta*mu^3*theta^2 + beta*gamma10*mu^4 + beta*gamma30*mu^4 + beta*mu^4*nu + beta*gamma10*gamma30*mu^3 + beta*gamma10*mu^3*nu + beta*gamma10*mu*rho^3 + 3*beta*gamma10*mu^3*rho + beta*gamma30*mu*rho^3 + 3*beta*gamma30*mu^3*rho + beta*gamma30*mu^3*sigma + beta*gamma10*mu^3*tau + 2*beta*gamma30*mu^3*tau + 2*beta*gamma10*mu^3*theta + 2*beta*gamma30*mu^3*theta + beta*mu*nu*rho^3 + 3*beta*mu^3*nu*rho + beta*mu^3*nu*sigma + 2*beta*mu^3*nu*tau + 2*beta*mu^3*nu*theta + beta*mu*rho^3*sigma + 3*beta*mu^3*rho*sigma + 2*beta*mu*rho^3*tau + 6*beta*mu^3*rho*tau + 2*beta*mu*rho^3*theta + 6*beta*mu^3*rho*theta + 2*beta*mu^3*sigma*tau + 2*beta*mu^3*sigma*theta + 3*beta*mu^3*tau*theta + 3*beta*gamma10*mu^2*rho^2 + 3*beta*gamma30*mu^2*rho^2 + beta*gamma30*mu^2*tau^2 + beta*gamma10*mu^2*theta^2 + beta*gamma30*mu^2*theta^2 + 3*beta*mu^2*nu*rho^2 + beta*mu^2*nu*tau^2 + beta*mu^2*nu*theta^2 + 3*beta*mu^2*rho^2*sigma + beta*mu*rho^2*tau^2 + 2*beta*mu^2*rho*tau^2 + 6*beta*mu^2*rho^2*tau + beta*mu*rho^2*theta^2 + 2*beta*mu^2*rho*theta^2 + 6*beta*mu^2*rho^2*theta + beta*mu^2*sigma*tau^2 + beta*mu^2*sigma*theta^2 + beta*mu^2*tau*theta^2 + beta*mu^2*tau^2*theta + beta*gamma10*gamma30*mu*rho^2 + 2*beta*gamma10*gamma30*mu^2*rho + beta*gamma10*gamma30*mu^2*tau + beta*gamma10*gamma30*mu*theta^2 + 2*beta*gamma10*gamma30*mu^2*theta + beta*gamma10*mu*nu*rho^2 + 2*beta*gamma10*mu^2*nu*rho + beta*gamma10*mu^2*nu*tau + beta*gamma10*mu*nu*theta^2 + 2*beta*gamma10*mu^2*nu*theta + beta*gamma30*mu*rho^2*sigma + 2*beta*gamma30*mu^2*rho*sigma + beta*gamma10*mu*rho^2*tau + 2*beta*gamma10*mu^2*rho*tau + beta*gamma30*mu*rho*tau^2 + 2*beta*gamma30*mu*rho^2*tau + 4*beta*gamma30*mu^2*rho*tau + beta*gamma10*mu*rho*theta^2 + 2*beta*gamma10*mu*rho^2*theta + 4*beta*gamma10*mu^2*rho*theta + beta*gamma30*mu*rho*theta^2 + 2*beta*gamma30*mu*rho^2*theta + 4*beta*gamma30*mu^2*rho*theta + beta*gamma30*mu*sigma*tau^2 + 2*beta*gamma30*mu^2*sigma*tau + beta*gamma30*mu*sigma*theta^2 + 2*beta*gamma30*mu^2*sigma*theta + beta*gamma10*mu^2*tau*theta + beta*gamma30*mu*tau*theta^2 + beta*gamma30*mu*tau^2*theta + 3*beta*gamma30*mu^2*tau*theta + beta*mu*nu*rho^2*sigma + 2*beta*mu^2*nu*rho*sigma + beta*mu*nu*rho*tau^2 + 2*beta*mu*nu*rho^2*tau + 4*beta*mu^2*nu*rho*tau + beta*mu*nu*rho*theta^2 + 2*beta*mu*nu*rho^2*theta + 4*beta*mu^2*nu*rho*theta + beta*mu*nu*sigma*tau^2 + 2*beta*mu^2*nu*sigma*tau + beta*mu*nu*sigma*theta^2 + 2*beta*mu^2*nu*sigma*theta + beta*mu*nu*tau*theta^2 + beta*mu*nu*tau^2*theta + 3*beta*mu^2*nu*tau*theta + beta*mu*rho*sigma*tau^2 + 2*beta*mu*rho^2*sigma*tau + 4*beta*mu^2*rho*sigma*tau + beta*mu*rho*sigma*theta^2 + 2*beta*mu*rho^2*sigma*theta + 4*beta*mu^2*rho*sigma*theta + beta*mu*rho*tau*theta^2 + beta*mu*rho*tau^2*theta + 3*beta*mu*rho^2*tau*theta + 6*beta*mu^2*rho*tau*theta + beta*mu*sigma*tau*theta^2 + beta*mu*sigma*tau^2*theta + 3*beta*mu^2*sigma*tau*theta + beta*gamma10*gamma30*mu*rho*tau + 2*beta*gamma10*gamma30*mu*rho*theta + beta*gamma10*gamma30*mu*tau*theta + beta*gamma10*mu*nu*rho*tau + 2*beta*gamma10*mu*nu*rho*theta + beta*gamma10*mu*nu*tau*theta + 2*beta*gamma30*mu*rho*sigma*tau + 2*beta*gamma30*mu*rho*sigma*theta + beta*gamma10*mu*rho*tau*theta + 3*beta*gamma30*mu*rho*tau*theta + 2*beta*gamma30*mu*sigma*tau*theta + 2*beta*mu*nu*rho*sigma*tau + 2*beta*mu*nu*rho*sigma*theta + 3*beta*mu*nu*rho*tau*theta + 2*beta*mu*nu*sigma*tau*theta + 3*beta*mu*rho*sigma*tau*theta)
%I3
I3mat = -(b*mu^3*tau*theta + b*rho^3*tau*theta + b*mu*tau^2*theta^2 + b*mu^2*tau*theta^2 + b*mu^2*tau^2*theta + b*rho*tau^2*theta^2 + b*rho^2*tau*theta^2 + b*rho^2*tau^2*theta + b*sigma*tau^2*theta^2 + b*gamma10*mu*tau*theta^2 + b*gamma10*mu^2*tau*theta - b*beta*sigma*tau*theta^2 - b*beta*sigma*tau^2*theta + b*gamma10*rho*tau*theta^2 + b*gamma10*rho^2*tau*theta + b*gamma10*sigma*tau*theta^2 + 2*b*mu*rho*tau*theta^2 + 2*b*mu*rho*tau^2*theta + 3*b*mu*rho^2*tau*theta + 3*b*mu^2*rho*tau*theta + b*mu*sigma*tau*theta^2 + b*mu*sigma*tau^2*theta + b*mu^2*sigma*tau*theta + b*rho*sigma*tau*theta^2 + b*rho*sigma*tau^2*theta + b*rho^2*sigma*tau*theta - b*beta*mu*sigma*tau*theta + 2*b*gamma10*mu*rho*tau*theta - b*beta*rho*sigma*tau*theta + b*gamma10*mu*sigma*tau*theta + b*gamma10*rho*sigma*tau*theta + 2*b*mu*rho*sigma*tau*theta)/(beta*mu^5 + beta*mu*rho^4 + 4*beta*mu^4*rho + beta*mu^4*sigma + 2*beta*mu^4*tau + 2*beta*mu^4*theta + 4*beta*mu^2*rho^3 + 6*beta*mu^3*rho^2 + beta*mu^3*tau^2 + beta*mu^3*theta^2 + beta*gamma10*mu^4 + beta*gamma30*mu^4 + beta*mu^4*nu + beta*gamma10*gamma30*mu^3 + beta*gamma10*mu^3*nu + beta*gamma10*mu*rho^3 + 3*beta*gamma10*mu^3*rho + beta*gamma30*mu*rho^3 + 3*beta*gamma30*mu^3*rho + beta*gamma30*mu^3*sigma + beta*gamma10*mu^3*tau + 2*beta*gamma30*mu^3*tau + 2*beta*gamma10*mu^3*theta + 2*beta*gamma30*mu^3*theta + beta*mu*nu*rho^3 + 3*beta*mu^3*nu*rho + beta*mu^3*nu*sigma + 2*beta*mu^3*nu*tau + 2*beta*mu^3*nu*theta + beta*mu*rho^3*sigma + 3*beta*mu^3*rho*sigma + 2*beta*mu*rho^3*tau + 6*beta*mu^3*rho*tau + 2*beta*mu*rho^3*theta + 6*beta*mu^3*rho*theta + 2*beta*mu^3*sigma*tau + 2*beta*mu^3*sigma*theta + 3*beta*mu^3*tau*theta + 3*beta*gamma10*mu^2*rho^2 + 3*beta*gamma30*mu^2*rho^2 + beta*gamma30*mu^2*tau^2 + beta*gamma10*mu^2*theta^2 + beta*gamma30*mu^2*theta^2 + 3*beta*mu^2*nu*rho^2 + beta*mu^2*nu*tau^2 + beta*mu^2*nu*theta^2 + 3*beta*mu^2*rho^2*sigma + beta*mu*rho^2*tau^2 + 2*beta*mu^2*rho*tau^2 + 6*beta*mu^2*rho^2*tau + beta*mu*rho^2*theta^2 + 2*beta*mu^2*rho*theta^2 + 6*beta*mu^2*rho^2*theta + beta*mu^2*sigma*tau^2 + beta*mu^2*sigma*theta^2 + beta*mu^2*tau*theta^2 + beta*mu^2*tau^2*theta + beta*gamma10*gamma30*mu*rho^2 + 2*beta*gamma10*gamma30*mu^2*rho + beta*gamma10*gamma30*mu^2*tau + beta*gamma10*gamma30*mu*theta^2 + 2*beta*gamma10*gamma30*mu^2*theta + beta*gamma10*mu*nu*rho^2 + 2*beta*gamma10*mu^2*nu*rho + beta*gamma10*mu^2*nu*tau + beta*gamma10*mu*nu*theta^2 + 2*beta*gamma10*mu^2*nu*theta + beta*gamma30*mu*rho^2*sigma + 2*beta*gamma30*mu^2*rho*sigma + beta*gamma10*mu*rho^2*tau + 2*beta*gamma10*mu^2*rho*tau + beta*gamma30*mu*rho*tau^2 + 2*beta*gamma30*mu*rho^2*tau + 4*beta*gamma30*mu^2*rho*tau + beta*gamma10*mu*rho*theta^2 + 2*beta*gamma10*mu*rho^2*theta + 4*beta*gamma10*mu^2*rho*theta + beta*gamma30*mu*rho*theta^2 + 2*beta*gamma30*mu*rho^2*theta + 4*beta*gamma30*mu^2*rho*theta + beta*gamma30*mu*sigma*tau^2 + 2*beta*gamma30*mu^2*sigma*tau + beta*gamma30*mu*sigma*theta^2 + 2*beta*gamma30*mu^2*sigma*theta + beta*gamma10*mu^2*tau*theta + beta*gamma30*mu*tau*theta^2 + beta*gamma30*mu*tau^2*theta + 3*beta*gamma30*mu^2*tau*theta + beta*mu*nu*rho^2*sigma + 2*beta*mu^2*nu*rho*sigma + beta*mu*nu*rho*tau^2 + 2*beta*mu*nu*rho^2*tau + 4*beta*mu^2*nu*rho*tau + beta*mu*nu*rho*theta^2 + 2*beta*mu*nu*rho^2*theta + 4*beta*mu^2*nu*rho*theta + beta*mu*nu*sigma*tau^2 + 2*beta*mu^2*nu*sigma*tau + beta*mu*nu*sigma*theta^2 + 2*beta*mu^2*nu*sigma*theta + beta*mu*nu*tau*theta^2 + beta*mu*nu*tau^2*theta + 3*beta*mu^2*nu*tau*theta + beta*mu*rho*sigma*tau^2 + 2*beta*mu*rho^2*sigma*tau + 4*beta*mu^2*rho*sigma*tau + beta*mu*rho*sigma*theta^2 + 2*beta*mu*rho^2*sigma*theta + 4*beta*mu^2*rho*sigma*theta + beta*mu*rho*tau*theta^2 + beta*mu*rho*tau^2*theta + 3*beta*mu*rho^2*tau*theta + 6*beta*mu^2*rho*tau*theta + beta*mu*sigma*tau*theta^2 + beta*mu*sigma*tau^2*theta + 3*beta*mu^2*sigma*tau*theta + beta*gamma10*gamma30*mu*rho*tau + 2*beta*gamma10*gamma30*mu*rho*theta + beta*gamma10*gamma30*mu*tau*theta + beta*gamma10*mu*nu*rho*tau + 2*beta*gamma10*mu*nu*rho*theta + beta*gamma10*mu*nu*tau*theta + 2*beta*gamma30*mu*rho*sigma*tau + 2*beta*gamma30*mu*rho*sigma*theta + beta*gamma10*mu*rho*tau*theta + 3*beta*gamma30*mu*rho*tau*theta + 2*beta*gamma30*mu*sigma*tau*theta + 2*beta*mu*nu*rho*sigma*tau + 2*beta*mu*nu*rho*sigma*theta + 3*beta*mu*nu*rho*tau*theta + 2*beta*mu*nu*sigma*tau*theta + 3*beta*mu*rho*sigma*tau*theta)

lambda = beta*(I1mat+I2mat)/N;
simplify(lambda)
amodif = gamma10*mu^2 + gamma10*rho^2 + 3*mu*rho^2 + 3*mu^2*rho + mu^2*sigma + mu^2*tau + mu^2*theta + rho^2*sigma + rho^2*tau + rho^2*theta + mu^3 + rho^3 - beta*mu*sigma + 2*gamma10*mu*rho - beta*rho*sigma + gamma10*mu*sigma + gamma10*mu*theta - beta*sigma*tau - beta*sigma*theta + gamma10*rho*sigma + gamma10*rho*theta + gamma10*sigma*theta + 2*mu*rho*sigma + 2*mu*rho*tau + 2*mu*rho*theta + mu*sigma*tau + mu*sigma*theta + mu*tau*theta + rho*sigma*tau + rho*sigma*theta + rho*tau*theta + sigma*tau*theta;
numerateur = -(gamma30+mu+nu+rho)*((mu+rho+sigma)*(mu+rho+theta)*(gamma10+mu+rho+tau)-beta*sigma*(rho+tau+theta+mu));
denominateur = gamma10*mu^2 + gamma30*mu^2 + gamma10*rho^2 + gamma30*rho^2 + mu^2*nu + 3*mu*rho^2 + 3*mu^2*rho + nu*rho^2 + mu^2*sigma + mu^2*tau + mu^2*theta + rho^2*sigma + rho^2*tau + rho^2*theta + mu^3 + rho^3 + gamma10*gamma30*mu + gamma10*gamma30*rho + gamma10*gamma30*theta + gamma10*mu*nu + 2*gamma10*mu*rho + 2*gamma30*mu*rho + gamma10*nu*rho + gamma30*mu*sigma + gamma30*mu*tau + gamma10*mu*theta + gamma30*mu*theta + gamma10*nu*theta + gamma30*rho*sigma + gamma30*rho*tau + gamma10*rho*theta + gamma30*rho*theta + 2*mu*nu*rho + gamma30*sigma*tau + gamma30*sigma*theta + mu*nu*sigma + gamma30*tau*theta + mu*nu*tau + mu*nu*theta + 2*mu*rho*sigma + 2*mu*rho*tau + 2*mu*rho*theta + nu*rho*sigma + nu*rho*tau + mu*sigma*tau + nu*rho*theta + mu*sigma*theta + mu*tau*theta + nu*sigma*tau + nu*sigma*theta + nu*tau*theta + rho*sigma*tau + rho*sigma*theta + rho*tau*theta + sigma*tau*theta;
lambdamat2 = numerateur/denominateur;
lambdamat = -((gamma30 + mu + nu + rho)*(amodif))/(gamma10*mu^2 + gamma30*mu^2 + gamma10*rho^2 + gamma30*rho^2 + mu^2*nu + 3*mu*rho^2 + 3*mu^2*rho + nu*rho^2 + mu^2*sigma + mu^2*tau + mu^2*theta + rho^2*sigma + rho^2*tau + rho^2*theta + mu^3 + rho^3 + gamma10*gamma30*mu + gamma10*gamma30*rho + gamma10*gamma30*theta + gamma10*mu*nu + 2*gamma10*mu*rho + 2*gamma30*mu*rho + gamma10*nu*rho + gamma30*mu*sigma + gamma30*mu*tau + gamma10*mu*theta + gamma30*mu*theta + gamma10*nu*theta + gamma30*rho*sigma + gamma30*rho*tau + gamma10*rho*theta + gamma30*rho*theta + 2*mu*nu*rho + gamma30*sigma*tau + gamma30*sigma*theta + mu*nu*sigma + gamma30*tau*theta + mu*nu*tau + mu*nu*theta + 2*mu*rho*sigma + 2*mu*rho*tau + 2*mu*rho*theta + nu*rho*sigma + nu*rho*tau + mu*sigma*tau + nu*rho*theta + mu*sigma*theta + mu*tau*theta + nu*sigma*tau + nu*sigma*theta + nu*tau*theta + rho*sigma*tau + rho*sigma*theta + rho*tau*theta + sigma*tau*theta);

Rp = sigma*beta*(tau+theta+rho+mu)/((theta+rho+mu)*(gamma10+rho+tau+mu)*(sigma+rho+mu));
S  = b/(mu*Rp); %ok

%% utility
clear all
syms b lambda beta sigma gamma10 gamma30 tau theta nu mu S E I1 I2 I3 rho

Rp = sigma*beta*(tau+theta+rho+mu)/((theta+rho+mu)*(gamma10+rho+tau+mu)*(sigma+rho+mu));
P = 1-1/Rp;
syms c
U = rho*(P-c);
dU = diff(U,rho);
solve(dU==0,rho)

%%
clear all;
beta=4;sigma=1/12;gamma10=6/7;gamma30=1;mu=1/35;nu=1/5;b=2;theta=10; tau=1; 
param.beta=beta;param.sigma=sigma;param.gamma1=gamma10;param.gamma3=gamma30;param.theta=theta;param.tau=tau;
syms rho
%let's find rho'
%R0 = sigma*beta*(tau+theta+mu)/((theta+mu)*(gamma10+tau+mu)*(sigma+mu));
Rpfun  = @(rho) ((sigma*beta*(tau+theta+rho+mu)./((theta+rho+mu).*(gamma10+rho+tau+mu).*(sigma+rho+mu)))-1);
alpha = fzero(Rpfun, 0); 
param.alpha=alpha;
%[U0,c0] = U_SEIIIS(b, beta, sigma, gamma10, gamma30, tau, theta, nu, mu, 0, 0) ;
%[U0,c1] = U_SEIIIS(b, beta, sigma, gamma10, gamma30, tau, theta, nu, mu, alpha, 0) ;
[U0,c0] = U_SEIIIS_v3(param,mu,0,0,1);
[U0,c1] = U_SEIIIS_v3(param,mu,param.alpha,0,1);

options = optimset('Display','off'); %options for minsearch
i=1; 
%vecC =(c1-(c0-c1)/2):(c0-c1)/200:(c0+(c0-c1)/2);
%vecC = [-0.45,-0.439,-0.38,-0.37];
vecC=-0.38
vecRhomax1 = [];
for c=vecC    
    %fun = @(rho) -U_SEIIIS_v3(param,mu,rho,c,1) ;
    %U_SEIIIS(b, beta, sigma, gamma10, gamma30, tau, theta, nu, mu, rho, c);
    %vecRhomax1(i) = min(max(fminsearch(fun,0,options),0),alpha)
    %vecRhomax1(i) = fminsearch(fun,0,options);
    %below: we solve dU/drho==0
    fun = @(rho) (1 - rho*(((mu + rho + sigma)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) + ((mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) + ((mu + rho + sigma)*(mu + rho + theta))/(beta*sigma*(mu + rho + tau + theta)) - ((mu + rho + sigma)*(mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)^2)) - ((mu + rho + sigma)*(mu + rho + theta)*(gamma10 + mu + rho + tau))/(beta*sigma*(mu + rho + tau + theta)) - c);
    vecRhomax1(i) = min(max(fsolve(fun,0,options),0),alpha);
    i=i+1; %c
end
%plot(vecC,vecRhomax1)

plot_paper1_procedure2(vecRhomax1, c1, c0, alpha, vecC,'')


%%
vecRho=0:0.0005:0.1; 
vecC=-[0.45,0.43,0.42,0.4,0.38,0.36];
for c=vecC
U=zeros(1,length(vecRho));
i=1;
for rho=vecRho  
    U(i) = U_SEIIIS_v3(param,mu,rho,c,1) ;
    Rp = sigma*beta*(tau+theta+rho+mu)/((theta+rho+mu)*(gamma10+rho+tau+mu)*(sigma+rho+mu));
    P = 1-1/Rp;
    %U(i) = rho.*(P-c);
    i=i+1; %c
end
figure(2)
plot(vecRho,U, 'DisplayName',num2str(c))
hold on
end
legend

plot([param.alpha,param.alpha],[0.04,0])

%%
rho = vecRhomax1;
Rp = sigma*beta*(tau+theta+rho+mu)./((theta+rho+mu).*(gamma10+rho+tau+mu).*(sigma+rho+mu));
numerateur = -(gamma30+mu+nu+rho).*((mu+rho+sigma).*(mu+rho+theta).*(gamma10+mu+rho+tau)-beta*sigma*(rho+tau+theta+mu));
denominateur = gamma10*mu^2 + gamma30*mu^2 + gamma10*rho.^2 + gamma30*rho.^2 + mu^2*nu + 3*mu*rho.^2 + 3*mu^2*rho + nu*rho.^2 + mu^2*sigma + mu^2*tau + mu^2*theta + rho.^2*sigma + rho.^2*tau + rho.^2*theta + mu^3 + rho.^3 + gamma10*gamma30*mu + gamma10*gamma30*rho + gamma10*gamma30*theta + gamma10*mu*nu + 2*gamma10*mu*rho + 2*gamma30*mu*rho + gamma10*nu*rho + gamma30*mu*sigma + gamma30*mu*tau + gamma10*mu*theta + gamma30*mu*theta + gamma10*nu*theta + gamma30*rho*sigma + gamma30*rho*tau + gamma10*rho*theta + gamma30*rho*theta + 2*mu*nu*rho + gamma30*sigma*tau + gamma30*sigma*theta + mu*nu*sigma + gamma30*tau*theta + mu*nu*tau + mu*nu*theta + 2*mu*rho*sigma + 2*mu*rho*tau + 2*mu*rho*theta + nu*rho*sigma + nu*rho*tau + mu*sigma*tau + nu*rho*theta + mu*sigma*theta + mu*tau*theta + nu*sigma*tau + nu*sigma*theta + nu*tau*theta + rho*sigma*tau + rho*sigma*theta + rho*tau*theta + sigma*tau*theta;
lambda2 = numerateur./denominateur
ratio = vecRhomax1.*Rp./lambda2.*(vecC>=c1); %rhoN/LambdaS
plot_paper1_procedure2(ratio, c1, c0, max(ratio), vecC)
