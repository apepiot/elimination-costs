%script to test the speed of the four diseases model
clear all;
%parameters
mu=1/35;
R1     = 1.08;
sigma1 = 365/11;
nu1    = 1/1.36;
gamma1 = 365/(31.5-11); 
eps1   = 0.11;
beta1  = R1*(sigma1+mu)*(gamma1+nu1+mu)*(nu1+mu)/(sigma1*(gamma1*(1-eps1)+mu+nu1));
souslaracine = @(gamma0,beta,nu,p,sigma) (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2);
alpha1 = max((beta1*eps1*sigma1 + sqrt(souslaracine(gamma1,beta1,nu1,eps1,sigma1)))/(2*(gamma1+mu+nu1)) - (2*mu+nu1+sigma1)/2,0);

%SEIIS2 %gono
R2     = 1.08;
sigma2 = 365/5;
nu2    = 12/6;
eps2   = 0.5;
gamma2 = 365/5; 
beta2  = R2*(sigma2+mu)*(gamma2+nu2+mu)*(nu2+mu)/(sigma2*(gamma2*(1-eps2)+mu+nu2));
alpha2 = (beta2*eps2*sigma2 + sqrt(souslaracine(gamma2,beta2,nu2,eps2,sigma2)))/(2*(gamma2+mu+nu2)) - (2*mu+nu2+sigma2)/2;

%SICR 
PHIV = 0.161;
RSICR     = 1/(1-PHIV);%2.2; %P=1-1/R so R=1/(1-P);
sigma3 = 365/(8.2*7); %(1-0.39)*365/(8.2*7)%0.6*(); %8.2 semaines
theta3 = 1/7;%1/9.8;
gamma3 = 0;%0.39*365/(8.2*7); %0.4*(); %proportion de diag au stade precoce : 40
ratioBeta = 9;
betaC3 = RSICR*(sigma3+gamma3+mu)*(theta3+mu)/(ratioBeta*(theta3+mu)+1);
betaI3 = ratioBeta*betaC3;
alpha3 = betaI3/2 - gamma3/2 - mu - sigma3/2 - theta3/2 + (betaI3^2 - 2*betaI3*gamma3 - 2*betaI3*sigma3 + 2*betaI3*theta3 +...
    gamma3^2 + 2*gamma3*sigma3 - 2*gamma3*theta3+ sigma3^2 - 2*sigma3*theta3 + 4*betaC3*sigma3 + theta3^2)^(1/2)/2;

%SEIIIS/S(syphilis)
PS = 0.077;
RS      = 1/(1-PS);
sigma4  = 365/25;
tau4    = 365/45;%(1-0.55*0.31)*365/45; %0.6*
theta4  = 12/(3.6);
gamma14 = 0;%(0.55*0.31)*365/45; %0.2*
gamma34 = 1/20;%1/5;
nu4     = 0;
beta4 = RS*(theta4+mu)*(gamma14+tau4+mu)*(sigma4+mu)/(sigma4*(tau4+theta4+mu));
Rpfun  = @(rho) (sigma4*beta4*(tau4+theta4+rho+mu)./((theta4+rho+mu).*(gamma14+rho+tau4+mu).*(sigma4+rho+mu))-1);
alpha4 = fzero(Rpfun, 0); 
b=2;
%%
options = optimset('Display','off'); %options for minsearch  
vecC=linspace(-0.3,1,10); tic
vecTimes = zeros(length(vecC),1);
vecRhohat = zeros(length(vecC),1);
i=1; 
for c=vecC
    fun1234 = @(rho) -U1234_SEIIS2SICRSEIIIS(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,...
            betaI3,betaC3,gamma3,sigma3,theta3,beta4,sigma4,tau4,gamma14,theta4,gamma34,nu4,mu,b,rho,c,alpha1,alpha2,alpha3,alpha4,1);
       
    vecRhohat(i) = min(max(fminsearch(fun1234,0,options),0),max([alpha1,alpha2,alpha3,alpha4]));
    disp(c)
    vecTimes(i) = toc;
    i=i+1;
end

%save('ws_4dis_003.mat')
%load('ws_4dis_002.mat')
