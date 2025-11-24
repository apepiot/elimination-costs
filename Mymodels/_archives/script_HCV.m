%% script HCV model (SEIR extend)
close all; 
clear all;

[beta,zeta,gamma,sigma,s,omega,b,mu,rho] = random_parameters('false', 'false');
gammap=gamma+s*rho;

R0 = beta*(omega*sigma+mu+gamma)/((sigma+mu)*(mu+gamma)) % a verifier theoriquement

pop0=[0.8,0.01,0.01,0,0];
MaxTime = 100;
[t, pop] = ode45(@(t,y) HCV(t, y, beta, gammap,sigma,zeta,omega,b,mu,'frequency'),[0 MaxTime],pop0);
%b/mu
%sum(pop(end,:))
pop(end,:)

% Prevalence at 0 when rho=alpha
%alpha = (-(gamma+mu) - beta*omega*sigma/(beta-sigma-mu))/s
alpha = (R0-1)*(sigma+mu)*(gamma+mu)/(sigma+mu-beta)/s;
%% plot
plot(t,pop)
legend('S','I','C','A','T')

%% Endemic equilibrium (numerically)
%numerically
S=pop(end,1); I=pop(end,2); C=pop(end,3);A=pop(end,4);T=pop(end,5);
N = sum(pop(end,:));
lambda = beta*(C+I)/N;
% r = mu*(omega*sigma+mu+gamma)/((mu+gamma)*(sigma+mu));
% b - (beta*r+mu)*S + mu* beta/b*r*S^2
% b/(mu*beta*r)*mu

%% Endemic equilibrium (theoretically)
%theoretically : S I C A T
R0 = beta*(omega*sigma+mu+gammap)/((sigma+mu)*(mu+gammap)) % a verifier theoriquement
S = b/(mu*R0) %S
I = b/(sigma+mu)*(1-1/R0)%I
C = omega*sigma/(omega*sigma+gammap+mu)*b*(R0-1)/beta %C
T = b/beta*(R0-1)*gammap*omega*sigma/((omega*sigma+gammap+mu)*(mu+zeta)) %T
b/(sigma+mu)*(1-1/R0)*(1+omega*sigma/(mu+gammap)) %I+C
A = (1-omega)*sigma/mu*I +zeta/mu*T %A

%% plot Rp in function of rho
rho = 0:alpha/100:(alpha*1.2);
gammap = gamma+s*rho;
Rp = beta.*(omega.*sigma+mu+gammap)./((sigma+mu).*(mu+gammap));
plot(rho,Rp)
alpha

%% plot utility
c=0.05; %mu*(1/beta-1/(sigma+mu))

rho=0:alpha/100:(alpha*1.2);
gammap = gamma+s*rho;
Rp = beta.*(omega.*sigma+mu+gammap)./((sigma+mu).*(mu+gammap));

P = mu/(sigma+mu).*(1-1./Rp).*(1+(omega*sigma)./(mu+gammap));
U=rho.*(P+c);

figure(1)
plot(rho,P)
figure(2)
plot(rho,max(U,0))

%% maxU resolution symbolic
clear all
syms beta gamma s rho sigma omega b mu c
%s=1;
gammap = gamma+s*rho;
Rp = beta.*(omega.*sigma+mu+gammap)./((sigma+mu).*(mu+gammap));
P = mu/(sigma+mu).*(1-1./Rp).*(1+(omega*sigma)./(mu+gammap));
U=rho.*(P+c);

dU=diff(U,rho);
(mu*omega*rho*s*sigma*(((mu + sigma)*(gamma + mu + rho*s))/(beta*(gamma + mu + omega*sigma + rho*s)) - 1))/((mu + sigma)*(gamma + mu + rho*s)^2) - (mu*rho*((omega*sigma)/(gamma + mu + rho*s) + 1)*((s*(mu + sigma))/(beta*(gamma + mu + omega*sigma + rho*s)) - (s*(mu + sigma)*(gamma + mu + rho*s))/(beta*(gamma + mu + omega*sigma + rho*s)^2)))/(mu + sigma) - (mu*(((mu + sigma)*(gamma + mu + rho*s))/(beta*(gamma + mu + omega*sigma + rho*s)) - 1)*((omega*sigma)/(gamma + mu + rho*s) + 1))/(mu + sigma)
c - rho*((mu*((omega*sigma)/(gamma + mu + rho*s) + 1)*((s*(mu + sigma))/(beta*(gamma + mu + omega*sigma + rho*s)) - (s*(mu + sigma)*(gamma + mu + rho*s))/(beta*(gamma + mu + omega*sigma + rho*s)^2)))/(mu + sigma) - (mu*omega*s*sigma*(((mu + sigma)*(gamma + mu + rho*s))/(beta*(gamma + mu + omega*sigma + rho*s)) - 1))/((mu + sigma)*(gamma + mu + rho*s)^2)) - (mu*(((mu + sigma)*(gamma + mu + rho*s))/(beta*(gamma + mu + omega*sigma + rho*s)) - 1)*((omega*sigma)/(gamma + mu + rho*s) + 1))/(mu + sigma)
solve(dU==0,rho)

%
% A = 1-1/Rp;
% B = (1+omega*sigma/(gammap+mu));
% ubis = rho*mu/(sigma+mu)*A*B;
% tau = omega*sigma/(gammap+mu);
% du = mu/(sigma+mu)*((1+tau-(sigma+mu))-rho*(sigma+mu)*tau^2/(beta*sigma*omega*(tau+1)) - rho*(1-(sigma+mu)/(tau+1))*tau^2/(omega*sigma))
% solve(du==0,rho)

%% solving hatrho = rho'(=alpha)
rhohat = -(gamma+mu)/s + sqrt(beta*mu*sigma*omega*(gamma+mu)/(mu*(mu+sigma-beta)-beta*c*(mu+sigma)))/s
R0 = beta*(omega*sigma+mu+gamma)/((sigma+mu)*(mu+gamma))
alpha = (R0-1)*(sigma+mu)*(gamma+mu)/(sigma+mu-beta)/s;

solve(rhohat-alpha==0,c)
-(mu*(mu - beta + sigma)*(gamma*mu - beta*mu - beta*gamma + gamma*sigma + mu*sigma + mu^2 - beta*omega*sigma))/(beta^2*omega*sigma*(mu + sigma))


%% solution of maxU
% SOLUTION POUR c=0
%%-(gamma*mu - beta*mu - beta*gamma + gamma*sigma + (beta*omega*sigma*(gamma + mu)*(mu - beta + sigma))^(1/2) + mu*sigma + mu^2)/(mu*s - beta*s + s*sigma)
-(gamma*mu - beta*mu - beta*gamma + gamma*sigma - (beta*omega*sigma*(gamma + mu)*(mu - beta + sigma))^(1/2) + mu*sigma + mu^2)/(mu*s - beta*s + s*sigma)
%forme simplifie
(sqrt(beta*omega*sigma*(gamma+mu)/(mu+sigma-beta)) - gamma-mu)/s 
%(sqrt(R0-beta/(mu+sigma-beta))-1)*(gamma+mu)/s

% SOLUTION POUR c neq 0
-(beta*mu^2 - gamma*mu^2 - mu^2*sigma + (-beta*mu*omega*sigma*(gamma + mu)*(beta*mu - mu*sigma - mu^2 + beta*c*mu + beta*c*sigma))^(1/2) - mu^3 + beta*gamma*mu - gamma*mu*sigma + beta*c*mu^2 + beta*c*mu*sigma + beta*c*gamma*mu + beta*c*gamma*sigma)/(beta*mu*s - mu^2*s - mu*s*sigma + beta*c*mu*s + beta*c*s*sigma)
%%sol2 = -(beta*mu^2 - gamma*mu^2 - mu^2*sigma - (-beta*mu*omega*sigma*(gamma + mu)*(beta*mu - mu*sigma - mu^2 + beta*c*mu + beta*c*sigma))^(1/2) - mu^3 + beta*gamma*mu - gamma*mu*sigma + beta*c*mu^2 + beta*c*mu*sigma + beta*c*gamma*mu + beta*c*gamma*sigma)/(beta*mu*s - mu^2*s - mu*s*sigma + beta*c*mu*s + beta*c*s*sigma)

%forme generale simplifiee (hatrho>0)
-(gamma+mu)/s + sqrt(beta*mu*sigma*omega*(gamma+mu)/(mu*(mu+sigma-beta)-beta*c*(mu+sigma)))/s

%DETAILS des calculs
% souslaracine = -beta*mu*omega*sigma*(gamma + mu)*(beta*mu - mu*sigma - mu^2 + beta*c*mu + beta*c*sigma);
% num = mu*(gamma+mu)*(sigma+mu-beta)-beta*c*(gamma+mu)*(sigma+mu) + sqrt(souslaracine)
% denom = mu*(beta-mu-sigma)+beta*c*(mu+sigma);%ok
% res = num/denom/s
%-(gamma+mu)/s + sqrt(souslaracine)/denom/s
c=-mu*(sigma+mu-beta)*((mu+sigma-beta)*(gamma+mu)-beta*omega*sigma)/(beta^2*omega*sigma*(mu+sigma))
C1 = mu*(sigma+mu-beta)/beta/(mu+sigma);
c = C1*s*alpha/(beta*omega*sigma)*(sigma+mu-beta)

%% Plot de \hat\rho en fonction de c
%clear all; close all;
%[beta,zeta,gamma,sigma,s,omega,b,mu,rho] = random_parameters('false', 'false');
%clear rho;

%bounds of c
upperbound = mu*(1/beta-1/(sigma+mu))
lowerbound = mu*(1/beta-1/(sigma+mu) - omega*sigma/((mu+sigma)*(gamma+mu)))
step = (upperbound - lowerbound)/100;
c = lowerbound*1.2:step:upperbound*1.2;

R0 = beta*(omega*sigma+mu+gamma)/((sigma+mu)*(mu+gamma)) %to check
alpha = (R0-1)*(sigma+mu)*(gamma+mu)/(sigma+mu-beta)/s

rhohat0 = -(gamma+mu)/s + sqrt(beta*mu*sigma*omega*(gamma+mu)./(mu*(mu+sigma-beta)-beta*c*(mu+sigma)))/s;
rhohat = rhohat0.*(c>lowerbound & c<upperbound) + zeros(1,length(c)) + alpha.*(c>=upperbound)

gammap = gamma+s*rhohat;
Rp = beta.*(omega.*sigma+mu+gammap)./((sigma+mu).*(mu+gammap));

figure(3)
plot(c,min(rhohat,alpha),'LineWidth',2,'Color','k')
title('HCV: $\hat\rho$ in function of $c$' , 'Interpreter','latex')
hold on
plot(upperbound,0,'+','Linewidth',1,'Color','k')
text(upperbound*0.85,-0.6,'$\mu(\frac{1}{\beta}-\frac{1}{\sigma+\mu})$','Interpreter','latex')
plot(lowerbound,0,'+','Linewidth',1,'Color','k')
text(lowerbound*1.2,-0.6,'$\mu \left(\frac{1}{\beta}-\frac{1}{\sigma+\mu} - \frac{\omega \sigma}{(\mu+\sigma)(\gamma+\mu)}\right)$','Interpreter','latex')
plot(lowerbound*1.15,alpha,'+','Linewidth',1,'Color','k')
text(lowerbound*1.15-0.01,alpha,'$\rho\prime$','Interpreter','latex')
plot(0,0,'+','Linewidth',1,'Color','k')
text(-0.0025,-0.6,'$0$','Interpreter','latex')

xlabel('c','Interpreter','latex')
ylabel('$\hat \rho$','Interpreter','latex')
set(gca,'YTickLabel',[],'XTickLabel',[]);
ylim([0,alpha*1.2])
xlim([lowerbound*1.15,upperbound*1.2])


figure(4)
plot(c,max(Rp,1),'LineWidth',2,'Color','k')
title('HCV: $\mathtt R(\hat\rho)$ in function of $c$' , 'Interpreter','latex')
hold on
plot(upperbound,0.95,'+','Linewidth',1,'Color','k')
text(upperbound*0.85,0.95-0.01,'$\mu(\frac{1}{\beta}-\frac{1}{\sigma+\mu})$','Interpreter','latex')
plot(lowerbound,0.95,'+','Linewidth',1,'Color','k')
text(lowerbound*1.2,0.95-0.01,'$\mu \left(\frac{1}{\beta}-\frac{1}{\sigma+\mu} - \frac{\omega \sigma}{(\mu+\sigma)(\gamma+\mu)}\right)$','Interpreter','latex')
plot(lowerbound*1.15,R0,'+','Linewidth',1,'Color','k')
text(lowerbound*1.15-0.02,R0,'$\mathtt R(0)$','Interpreter','latex')
plot(0,0.95,'+','Linewidth',1,'Color','k')
text(-0.0025,0.95-0.02,'$0$','Interpreter','latex')

xlabel('c','Interpreter','latex')
ylabel('$\mathtt R(\hat\rho)$','Interpreter','latex')
set(gca,'YTickLabel',[],'XTickLabel',[]);
ylim([0.95,R0*1.05])
xlim([lowerbound*1.15,upperbound*1.2])
