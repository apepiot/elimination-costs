%Verification des résultats de l'article 1

%% SIR
clear all; %close all;

syms beta gamma b mu s rho S I R c

R0 = beta/(gamma+mu);
gammap = gamma+s*rho;
Rp  = beta/(gammap+mu);
N = b/mu;
alpha = beta/s*(1-1/R0);

eqn1 = b-beta*S*I/N - mu*S==0;
eqn2 = beta*S*I/N - gammap*I -mu*I==0;
eqn3 = gammap*I-mu*R==0;

sol = solve([eqn1,eqn2,eqn3],[S,I,R])

DFS = [b/mu,0,0];
ES  = [(b*gamma + b*mu + b*rho*s)/(beta*mu),...
    -(b*gamma - b*beta + b*mu + b*rho*s)/(beta*gamma + beta*mu + beta*rho*s),...
    -(b*gamma^2 + b*rho^2*s^2 - b*beta*gamma + b*gamma*mu - b*beta*rho*s + 2*b*gamma*rho*s + b*mu*rho*s)/(beta*mu^2 + beta*gamma*mu + beta*mu*rho*s)];

%verif de la prevalence
simplify(ES(2)/N - mu/beta*(Rp-1)) %OK

P = mu/beta*(Rp-1);
U = rho*(P-c);

%dU = diff(U,rho);
dU = (mu*(beta/(gamma + mu + rho*s) - 1))/beta - c - (mu*rho*s)/(gamma + mu + rho*s)^2;

%verif de rhohat
rhohat = solve(dU==0,rho);
simplify(rhohat - beta/(R0*s)*(sqrt(R0*mu/(mu+beta*c))-1)) %OK

% [beta,zeta,gamma,sigma,s,omega,b,mu,rho] = random_parameters('false', 'false');
% c = 100;
% - ((gamma + mu)*(((beta*mu)/((gamma + mu)*(mu + beta*c)))^(1/2) - 1))/s - (gamma*mu + (beta*mu*(gamma + mu)*(mu + beta*c))^(1/2) + mu^2 + beta*c*gamma + beta*c*mu)/(s*(mu + beta*c))
% - (gamma*mu - (beta*mu*(gamma + mu)*(mu + beta*c))^(1/2) + mu^2 + beta*c*gamma + beta*c*mu)/(s*(mu + beta*c)) - ((gamma + mu)*(((beta*mu)/((gamma + mu)*(mu + beta*c)))^(1/2) - 1))/s

%verif de c1 et c2
c1 = -mu*(1-1/R0)/beta; %OK
c2 = mu/beta*(R0-1); %OK
c = c1;
rhohat = beta/(R0*s)*(sqrt(R0*mu/(mu+beta*c))-1);

% c1
simplify(rhohat - alpha)
%(gamma - beta + mu)/s + ((gamma + mu)*((beta^2/(gamma + mu)^2)^(1/2) - 1))/s
simplify((gamma - beta + mu)/s + ((gamma + mu)*(R0 - 1))/s) %0

%c2
c = c2;
rhohat = beta/(R0*s)*(sqrt(R0*mu/(mu+beta*c))-1);
simplify(rhohat - 0)%0

%R(rhohat)
rhohat = beta/(R0*s)*(sqrt(R0*mu/(mu+beta*c))-1);
simplify(beta/(gamma+s*rhohat+mu) - sqrt(R0*(mu+beta*c)/mu)) %OK

%% SICAT 
clear all; %close all;

syms beta gamma b mu s rho S I C A T c omega zeta sigma

R0  = beta*(omega*sigma+gamma+mu)/((gamma+mu)*(sigma+gamma+mu));
gammap = gamma+s*rho;
Rp = beta*(omega*sigma+gammap+mu)/((gammap+mu)*(sigma+gammap+mu));
N = b/mu;

% verif de alpha
alpha = ((beta-sigma)/2-(gamma+mu)+sqrt((beta-sigma)^2/4+beta*omega*sigma))/s;
%rho = alpha; 
%gammap = gamma+s*rho;
%Rp = simplify(beta*(omega*sigma+gammap+mu)/((gammap+mu)*(sigma+gammap+mu)))%OK

% Function R(rho)

% [beta,zeta,gamma,sigma,s,omega,b,mu,rho] = random_parameters('true', 'true');
% rho = -5:0.1:5;
% gammap = gamma+s*rho;
% Rp = beta*(omega*sigma+gammap+mu)./((gammap+mu).*(sigma+gammap+mu));
% R0  = beta*(omega*sigma+gamma+mu)/((gamma+mu)*(sigma+gamma+mu))
% plot(rho,Rp)


% Solution du systeme et prevalence
lambda = beta*(I+C)/N;
eqn1 = b - lambda*S - mu*S==0;
eqn2 = lambda*S - (sigma+gammap)*I -mu*I==0;
eqn3 = omega*sigma*I-(mu+gammap)*C==0;
eqn4 = (1-omega)*sigma*I + zeta*T - mu*A;
eqn5 = gammap*(C+I) - (mu+zeta)*T;

sol = solve([eqn1,eqn2,eqn3,eqn4,eqn5],[S,I,C,A,T]);
ES = [(b*gamma^2 + 2*b*gamma*mu + 2*b*gamma*rho*s + b*sigma*gamma + b*mu^2 + 2*b*mu*rho*s + b*sigma*mu + b*rho^2*s^2 + b*sigma*rho*s)/(beta*mu*(gamma + mu + omega*sigma + rho*s)),...
    -(b*gamma^2 + b*mu^2 + b*rho^2*s^2 - b*beta*gamma - b*beta*mu + 2*b*gamma*mu + b*gamma*sigma + b*mu*sigma - b*beta*omega*sigma - b*beta*rho*s + 2*b*gamma*rho*s + 2*b*mu*rho*s + b*rho*s*sigma)/((beta*gamma + beta*mu + beta*sigma + beta*rho*s)*(gamma + mu + omega*sigma + rho*s)),...
    -(omega*sigma*(b*gamma^2 + b*mu^2 + b*rho^2*s^2 - b*beta*gamma - b*beta*mu + 2*b*gamma*mu + b*gamma*sigma + b*mu*sigma - b*beta*omega*sigma - b*beta*rho*s + 2*b*gamma*rho*s + 2*b*mu*rho*s + b*rho*s*sigma))/((gamma + mu + omega*sigma + rho*s)*(beta*gamma^2 + 2*beta*gamma*mu + 2*beta*gamma*rho*s + beta*sigma*gamma + beta*mu^2 + 2*beta*mu*rho*s + beta*sigma*mu + beta*rho^2*s^2 + beta*sigma*rho*s)),...
    -((b*gamma^2 + b*mu^2 + b*rho^2*s^2 - b*beta*gamma - b*beta*mu + 2*b*gamma*mu + b*gamma*sigma + b*mu*sigma - b*beta*omega*sigma - b*beta*rho*s + 2*b*gamma*rho*s + 2*b*mu*rho*s + b*rho*s*sigma)*(mu^2*sigma + gamma^2*zeta - mu^2*omega*sigma + rho^2*s^2*zeta + gamma*mu*sigma + gamma*mu*zeta + gamma*sigma*zeta + mu*sigma*zeta - gamma*mu*omega*sigma + mu*rho*s*sigma + 2*gamma*rho*s*zeta - mu*omega*sigma*zeta + mu*rho*s*zeta + rho*s*sigma*zeta - mu*omega*rho*s*sigma))/(mu*(gamma + mu + omega*sigma + rho*s)*(beta*mu^3 + beta*mu^2*sigma + beta*gamma^2*zeta + beta*mu^2*zeta + 2*beta*gamma*mu^2 + beta*gamma^2*mu + beta*gamma*mu*sigma + 2*beta*gamma*mu*zeta + beta*gamma*sigma*zeta + beta*mu*sigma*zeta + 2*beta*mu^2*rho*s + beta*mu*rho^2*s^2 + beta*rho^2*s^2*zeta + 2*beta*gamma*mu*rho*s + beta*mu*rho*s*sigma + 2*beta*gamma*rho*s*zeta + 2*beta*mu*rho*s*zeta + beta*rho*s*sigma*zeta)),...
    -(b*gamma^3 + b*rho^3*s^3 - b*beta*gamma^2 + b*gamma*mu^2 + 2*b*gamma^2*mu + b*gamma^2*sigma + b*gamma*mu*sigma + 3*b*gamma^2*rho*s + b*mu^2*rho*s - b*beta*rho^2*s^2 + 3*b*gamma*rho^2*s^2 + 2*b*mu*rho^2*s^2 + b*rho^2*s^2*sigma - b*beta*gamma*mu - b*beta*gamma*omega*sigma - 2*b*beta*gamma*rho*s - b*beta*mu*rho*s + 4*b*gamma*mu*rho*s + 2*b*gamma*rho*s*sigma + b*mu*rho*s*sigma - b*beta*omega*rho*s*sigma)/(beta*mu^3 + beta*mu^2*sigma + beta*gamma^2*zeta + beta*mu^2*zeta + 2*beta*gamma*mu^2 + beta*gamma^2*mu + beta*gamma*mu*sigma + 2*beta*gamma*mu*zeta + beta*gamma*sigma*zeta + beta*mu*sigma*zeta + 2*beta*mu^2*rho*s + beta*mu*rho^2*s^2 + beta*rho^2*s^2*zeta + 2*beta*gamma*mu*rho*s + beta*mu*rho*s*sigma + 2*beta*gamma*rho*s*zeta + 2*beta*mu*rho*s*zeta + beta*rho*s*sigma*zeta)];

Sth = b/(mu*Rp);
Ith = b/(sigma+mu+gammap)*(1-1/Rp);
Cth = b/beta*(Rp-1)*omega*sigma/(omega*sigma+mu+gammap);
Ath = b/mu*((1-omega)*sigma/(sigma+mu+gammap)*(1-1/Rp) + zeta*gammap*(Rp-1)/(mu+zeta)/beta);
Tth = b*gammap/beta/(mu+zeta)*(Rp-1);

[simplify(Sth -ES(1)),...
    simplify(Ith -ES(2)),...
    simplify(Cth -ES(3)),...
    simplify(Ath -ES(4)),...
    simplify(Tth -ES(5))] %OK

P = mu/beta*(Rp-1);
simplify(P - (ES(2)+ES(3))/N) %OK

U = rho*(P-c);

dU = diff(U,rho);
%dU = (mu*((beta*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)) - 1))/beta - c - (mu*rho*((beta*s*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)^2) - (beta*s)/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)) + (beta*s*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)^2*(gamma + mu + sigma + rho*s))))/beta

% verif de c1 et c2
alpha = ((beta-sigma)/2-(gamma+mu)+sqrt((beta-sigma)^2/4+beta*omega*sigma))/s;
c1 = -mu/beta*(2*s*alpha*sqrt((beta-sigma)^2/4+beta*omega*sigma)/(beta*(gamma+s*alpha+mu+omega*sigma))); %OK
c2 = mu/beta*(R0-1); %OK
%c1
rho = alpha;
dU = (mu*((beta*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)) - 1))/beta - c - (mu*rho*((beta*s*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)^2) - (beta*s)/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)) + (beta*s*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)^2*(gamma + mu + sigma + rho*s))))/beta
solve(dU==0,c)
simplify(c1 - ((mu*((beta*(beta/2 - sigma/2 + omega*sigma + ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2))) - 1))/beta + (mu*((beta*s*(beta/2 - sigma/2 + omega*sigma + ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2))^2) - (beta*s)/((beta/2 - sigma/2 + ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2))) + (beta*s*(beta/2 - sigma/2 + omega*sigma + ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2))^2*(beta/2 + sigma/2 + ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2))))*(gamma - beta/2 + mu + sigma/2 - ((beta - sigma)^2/4 + beta*omega*sigma)^(1/2)))/(beta*s)))
(mu*(beta^2 - 2*beta*sigma + sigma^2 + 4*beta*omega*sigma)^(1/2)*(2*gamma - beta + 2*mu + sigma - (beta^2 - 2*beta*sigma + sigma^2 + 4*beta*omega*sigma)^(1/2)))/(beta^2*(beta - sigma + (beta^2 - 2*beta*sigma + sigma^2 + 4*beta*omega*sigma)^(1/2) + 2*omega*sigma)) - (2*mu*(gamma - beta/2 + mu + sigma/2 - (beta^2 - 2*beta*sigma + sigma^2 + 4*beta*omega*sigma)^(1/2)/2)*(beta^2 - 2*beta*sigma - sigma*(beta^2 - 2*beta*sigma + sigma^2 + 4*beta*omega*sigma)^(1/2) + sigma^2 + beta*(beta^2 - 2*beta*sigma + sigma^2 + 4*beta*omega*sigma)^(1/2) + 2*omega*sigma*(beta^2 - 2*beta*sigma + sigma^2 + 4*beta*omega*sigma)^(1/2) + 4*beta*omega*sigma))/(beta*((beta^2 - 2*beta*sigma + sigma^2 + 4*beta*omega*sigma)^(3/2) + beta^2*(beta^2 - 2*beta*sigma + sigma^2 + 4*beta*omega*sigma)^(1/2) - sigma^2*(beta^2 - 2*beta*sigma + sigma^2 + 4*beta*omega*sigma)^(1/2) + 2*beta*sigma^2 - 4*beta^2*sigma + 2*beta^3 - 4*beta*omega*sigma^2 + 8*beta^2*omega*sigma + 4*beta*omega^2*sigma^2))
%OK

%c2
rho = 0;
dU = (mu*((beta*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)) - 1))/beta - c - (mu*rho*((beta*s*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)^2) - (beta*s)/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)) + (beta*s*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)^2*(gamma + mu + sigma + rho*s))))/beta
solve(dU==0,c)
(mu*((beta*(gamma + mu + omega*sigma))/((gamma + mu)*(gamma + mu + sigma)) - 1))/beta -c2 %OK
 

%% SICAT : c1 
clear all;
cond = 1;
t=1
while(cond & t<1000000)

[beta,zeta,gamma,sigma,s,omega,b,mu,rho] = random_parameters('true', 'true');

R0 = beta.*(omega*sigma+mu+gamma)./((sigma+gamma+mu).*(mu+gamma));
alpha = ((beta-sigma)/2 - (gamma+mu) + sqrt(((beta-sigma)/2)^2 + beta*omega*sigma))/s;

c11 = -mu/beta*(2*s*alpha*sqrt((beta-sigma)^2/4+beta*omega*sigma)/(beta*(gamma+s*alpha+mu+omega*sigma)));
c12 = -mu/beta;
c2 = mu/beta*(R0-1);

cond = c12<=c11;
t=t+1
end
%%
c1 = max(c11,c12);

i=1; vecC =(c1-(c2-c1)/2):(c2-c1)/1000:(c2+(c2-c1)/2);
vecRhomax1 = []; vecRhomax2 = [];vecRhomax3 = [];
for c=vecC    
    fun = @(rho) U_HCV2(rho, beta, gamma, s, sigma, omega, mu, c);
    vecRhomax1(i) = max(fmincon(fun,0),0);   
    i=i+1; %c+1
end


[c11,c12,c1,c2]
plot(vecC,vecRhomax1)
hold on
plot(vecC,alpha*ones(length(vecC)))
hold on
plot(c11*ones(1,100), 0:(alpha*1.1/99):alpha*1.1)
hold on
plot(c12*ones(1,100), 0:(alpha*1.1/99):alpha*1.1)
hold on
plot(c2*ones(1,100), 0:(alpha*1.1/99):alpha*1.1)
legend('$\hat\rho$','$\rho\prime$','$c_{11}$', '$c_{12}$','$c_2$','Interpreter','latex')


