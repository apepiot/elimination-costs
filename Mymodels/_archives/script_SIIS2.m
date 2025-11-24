%script for the SIAS model
[beta1,beta2,gamma10,gamma20,s,s2,b,mu] = random_parameters('true', 'true');
omega = 0.4;
rho = 0.5;
eps1 = 0.4;eps2=0.6;
sigma1 = 365/14;
sigma2 = 365/21;
nu1 = 1;
nu2 = 365/450;

% Parametres du systeme d'ODE 
tspan = 0:1:100;
Y0 = [99; 1; 1 ; 1 ; 1 ; 1 ; 1 ;1 ; 1];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) SIIS2(t,Y,b,beta1,beta2,gamma10,gamma20,rho,eps1,eps2,sigma1,sigma2,nu1,nu2,mu,'frequency'),tspan,Y0, options);
sum(Ys(end,:))
b/mu

Yend = Ys(end,:)

%% solving the ODE system
clear all

syms b beta1 beta2 gamma10 gamma20 rho eps1 eps2 sigma1 sigma2 nu1 nu2 mu
syms S IA0 IS0 I0A IAA ISA I0S IAS ISS

%rho=0;nu1=0;nu2=0;
gamma1 = gamma10+rho; gamma2=gamma20+rho;
eps1=0;

R1 = beta1*(eps1*sigma1 + gamma1+mu+nu1)/((gamma1+mu+nu1)*(eps1*sigma1+mu+(1-eps1)*(nu1+rho)));
R2 = beta2*(eps2*sigma2 + gamma2+mu+nu2)/((gamma2+mu+nu2)*(eps2*sigma2+mu+(1-eps2)*(nu2+rho)));

%Lambda1 = beta1*(N-S-I0A-I0S)/N;
%Lambda2 = beta2*(N-S-IA0-IS0)/N;
Lambda1 = beta1*(1-1/R1);
Lambda2 = beta2*(1-1/R2);
N = b/mu;

dS   = b - (Lambda1 + Lambda2)*S + (1-eps1)*nu1*IA0 + (1-eps2)*nu2*I0A +...
    (nu1+gamma10)*IS0 + (nu2+gamma20)*I0S + rho*(N-S) - eps1*rho*IA0 - eps2*rho*I0A - mu*S;
dIA0 = Lambda1*S - Lambda2*IA0 - ((1-eps1)*(nu1+rho) + eps1*sigma1 + mu)*IA0 + (1-eps2)*nu2*IAA + (nu2+gamma20)*IAS;
dI0A = Lambda2*S - Lambda1*I0A - ((1-eps2)*(nu2+rho) + eps2*sigma2 + mu)*I0A + (nu1+gamma10)*ISA +(1-eps1)*nu1*IAA;
dIS0 = eps1*sigma1*IA0 - Lambda2*IS0 + (1-eps2)*nu2*ISA - (nu1+gamma10+rho + mu)*IS0 + (nu2+gamma20)*ISS;
dI0S = eps2*sigma2*I0A - Lambda1*I0S + (nu1+gamma10)*ISS + (1-eps1)*nu1*IAS - (nu2+gamma20+rho + mu)*I0S;
dISA = eps1*sigma1*IAA + Lambda2*IS0 - ((1-eps2)*nu2 + nu1+gamma10 + eps2*sigma2 + rho + mu)*ISA;
dIAS = Lambda1*I0S + eps2*sigma2*IAA - ((1-eps1)*nu1 + nu2+gamma20 + eps1*sigma1 + rho + mu)*IAS;
dIAA = Lambda1*I0A + Lambda2*IA0 - ((1-eps1)*nu1 + eps1*sigma1 + (1-eps2)*nu2 + eps2*sigma2 + rho +mu)*IAA;
dISS = eps1*sigma1*IAS + eps2*sigma2*ISA - (nu1+gamma10+rho + nu2+gamma20 + mu)*ISS;

sol = solve([dS==0,dIA0==0,dI0A==0,dIS0==0,dI0S==0,dISA==0,dIAS==0,dIAA==0,dISS==0],[S ,IA0 ,IS0, I0A ,IAA ,ISA ,I0S ,IAS ,ISS]);

%% utility function
[beta1,beta2,gamma10,gamma20,s,s2,b,mu] = random_parameters('true', 'true');
omega = 0.4;
eps1 = 0.4;eps2=0.6;
sigma1 = 365/14;
sigma2 = 365/21;
nu1 = 1;
nu2 = 365/450;

% Parametres du systeme d'ODE 
tspan = 0:1:100;
Y0 = [99; 1; 1 ; 1 ; 1 ; 1 ; 1 ;1 ; 1];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) SIIS2(t,Y,b,beta1,beta2,gamma10,gamma20,rho,eps1,eps2,sigma1,sigma2,nu1,nu2,mu,'frequency'),tspan,Y0, options);

%% argmaxU
clear all; close all;
% [beta1,beta2,gamma10,gamma20,s,s2,b,mu] = random_parameters('true', 'true');
% omega = 0.4;
% eps1 = 0.4;eps2=0.6;
% sigma1 = 365/14;
% sigma2 = 365/21;
% nu1 = 1;
% nu2 = 365/450;

P1=0.074;eps1=11/100;sigma1=1/(11/365);gamma10=1/((30-11)/365);nu1=1/(497/365); mu=1/35;
beta1 = -((mu + eps1*sigma1 - nu1*(eps1 - 1))*(gamma10 + mu + nu1))/((P1 - 1)*(gamma10 + mu + nu1 + eps1*sigma1));
R10 = beta1*(eps1*sigma1 + gamma10+mu+nu1)/((gamma10+mu+nu1)*(eps1*sigma1+mu+(1-eps1)*nu1));

P2=0.07;eps2=70/100;sigma2=1/(5/365);gamma20=1/(5/365);nu2=1/(5/12);mu=1/35;
beta2 = -((mu + eps2*sigma2 - nu2*(eps2 - 1))*(gamma20 + mu + nu2))/((P2 - 1)*(gamma20 + mu + nu2 + eps2*sigma2));
R20 = beta2*(eps2*sigma2 + gamma20+mu+nu2)/((gamma20+mu+nu2)*(eps2*sigma2+mu+(1-eps2)*nu2));
b=2;

% %SIIS
alpha1 = - nu1 - (beta1 - gamma10 - 2*mu + gamma10*eps1 + mu*eps1 - eps1*sigma1 +...
    (gamma10^2*eps1^2 - 2*gamma10^2*eps1 + gamma10^2 + 2*gamma10*mu*eps1^2 - 2*gamma10*mu*eps1 +...
    2*gamma10*eps1^2*sigma1- 2*gamma10*eps1*sigma1 - 2*beta1*gamma10*eps1 + 2*beta1*gamma10 +...
    mu^2*eps1^2 + 2*mu*eps1^2*sigma1 - 2*beta1*mu*eps1 + eps1^2*sigma1^2 -...
    4*beta1*eps1^2*sigma1 + 2*beta1*eps1*sigma1 + beta1^2)^(1/2))/(2*(eps1 - 1));
alpha2 = - nu2 - (beta2 - gamma20 - 2*mu + gamma20*eps2 + mu*eps2 - eps2*sigma2 + (gamma20^2*eps2^2 - 2*gamma20^2*eps2 + gamma20^2 + 2*gamma20*mu*eps2^2 - 2*gamma20*mu*eps2 +...
    2*gamma20*eps2^2*sigma2- 2*gamma20*eps2*sigma2 - 2*beta2*gamma20*eps2 + 2*beta2*gamma20 + mu^2*eps2^2 + 2*mu*eps2^2*sigma2 - 2*beta2*mu*eps2 + eps2^2*sigma2^2 -...
    4*beta2*eps2^2*sigma2 + 2*beta2*eps2*sigma2 + beta2^2)^(1/2))/(2*(eps2 - 1));

% [U0,c0] = U_SIIS(p, beta, sigma, gamma0, mu, nu, 0, 0);
% [U1,c1] = U_SIIS(p, beta, sigma, gamma0, mu, nu, alpha, 0);
options = optimset('Display','off','MaxFunEvals',25); %options for minsearch
i=1; %vecC =(c1-(c2-c1)/2):(c2-c1)/1000:(c2+(c2-c1)/2);
vecC = [-10:0.5:-5.5,-5.4:0.1:9.5,9.6:0.5:12];
vecRhomax12=zeros(1,length(vecC));vecRhomax1=vecRhomax12;
vecRhomax2=vecRhomax12; vecRhomax=vecRhomax12;
for c=vecC    
    fun12 = @(rho) -U12_SIIS2(b,beta1,beta2,gamma10,gamma20,eps1,eps2,sigma1,sigma2,nu1,nu2,mu,rho,c);
    rho12 = min(max(fminsearch(fun12,0,options),0),max(alpha1,alpha2));
    vecRhomax12(i) = rho12;
    fun1 = @(rho) -U_SIIS(eps1, beta1, sigma1, gamma10, mu, nu1, rho, c);
    rho1 = min(max(fminsearch(fun1,0,options),0),alpha1);
    vecRhomax1(i) = rho1;
    fun2 = @(rho) -U_SIIS(eps2, beta2, sigma2, gamma20, mu, nu2, rho, c);
    rho2 = min(max(fminsearch(fun2,0,options),0),alpha2);
    vecRhomax2(i) = rho2;
    
    rho_cand = [rho12,rho1,rho2];
    [U_cand] = U12_SIIS2(b,beta1,beta2,gamma10,gamma20,eps1,eps2,sigma1,sigma2,nu1,nu2,mu,rho_cand,c);
    [Umax, imax] = max(U_cand);
    vecRhomax(i) = rho_cand(imax);
    
    i=i+1 %c
end

plot(vecC,vecRhomax, vecC, vecRhomax1, vecC, vecRhomax2)
legend('$\hat\rho$','$\hat\rho_1$', '$\hat\rho_2$','Interpreter','latex')
xlabel('$c$','Interpreter','latex')
ylabel('$\hat\rho$','Interpreter','latex')

title([{['$\beta_1$=',num2str(round(beta1,2)), ' $\beta_2$=',num2str(round(beta2,2)),' $\gamma_1(0)$=',num2str(round(gamma10,2)),' $\gamma_2(0)$=',...
    num2str(round(gamma20,2)), ' $\epsilon_1$=', num2str(eps1),' $\epsilon_2$=', num2str(eps2), ' $\mu$=', num2str(round(mu,2))]},...
    {['$\nu_1$=',num2str(round(nu1,2)), ' $\nu_2$=',num2str(round(nu2,2)),' $\sigma_1$=',num2str(round(sigma1,2)),' $\sigma_2$=',...
    num2str(round(sigma2,2)),' $\mathtt R_1(0)$=', num2str(round(R10,2)), ' $\mathtt R_2(0)=$', num2str(round(R20,2)),...
    ' $\rho_1\prime=$',  num2str(round(alpha1,2)),' $ \rho_2\prime=$',  num2str(round(alpha2,2))  ]}],...
    'Interpreter','latex')

% i=1; U=[]; vecRho=-1:0.1:10;
% for rho=vecRho  
%     U(i) = U_SIIS2(b,beta1,beta2,gamma10,gamma20,eps1,eps2,sigma1,sigma2,nu1,nu2,mu,rho,c);
%     i=i+1 %c
% end

%plot(vecRho,U)
