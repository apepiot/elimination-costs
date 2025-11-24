function [Rp,Lambdap,alpha,Ptot,Pun] = Rp_SICTP(betaI,betaC,theta0,gamma0,sigma,zeta,eta,p,mu,binput,rho)

if gamma0 ~=0
    error('the function has not been addressed when gamma0 is different from 0')
end

%%To change (input parameters as well)
Rp = (1-p)*(betaI*(mu+rho+theta0)+betaC*sigma)./((mu+rho+sigma)*(mu+rho+theta0)) + p*(1-zeta)*(betaI*(eta+mu+theta0)+betaC*sigma)./((eta+mu+sigma)*(eta+mu+theta0));

a = (1-p)./(mu+rho+sigma);
b = p*(1-zeta)./(eta+mu+sigma);
c = (1-p)./((sigma+rho+mu)*(theta0+rho+mu));
d = p*(1-zeta)./((eta+mu+sigma)*(eta+mu+theta0));
Lambdap=-(mu*zeta - 2*mu + mu*(a^2*betaI^2*zeta^2 - 2*a^2*betaI^2*zeta + a^2*betaI^2 - 2*a*b*betaI^2*zeta + 2*a*b*betaI^2 + 2*a*betaC*betaI*c*sigma*zeta^2 - 4*a*betaC*betaI*c*sigma*zeta + 2*a*betaC*betaI*c*sigma - 2*a*betaC*betaI*d*sigma*zeta + 2*a*betaC*betaI*d*sigma - 2*a*betaI*zeta^2 + 2*a*betaI*zeta + b^2*betaI^2 - 2*b*betaC*betaI*c*sigma*zeta + 2*b*betaC*betaI*c*sigma + 2*b*betaC*betaI*d*sigma - 2*b*betaI*zeta + betaC^2*c^2*sigma^2*zeta^2 - 2*betaC^2*c^2*sigma^2*zeta + betaC^2*c^2*sigma^2 - 2*betaC^2*c*d*sigma^2*zeta + 2*betaC^2*c*d*sigma^2 + betaC^2*d^2*sigma^2 - 2*betaC*c*sigma*zeta^2 + 2*betaC*c*sigma*zeta - 2*betaC*d*sigma*zeta + zeta^2).^(1/2) + a*betaI*mu + b*betaI*mu + betaC*c*mu*sigma + betaC*d*mu*sigma - a*betaI*mu*zeta - betaC*c*mu*sigma*zeta)./(2*(zeta - 1));
Lambdap=max(Lambdap,0);

%Prevalences
I  = binput*(1-p)*Lambdap./((Lambdap+mu)*(sigma+rho+mu));
C  = binput*Lambdap*sigma*(1-p)./((Lambdap+mu)*(sigma+rho+mu)*(theta0+rho+mu));
Ip = binput*Lambdap*p*(1-zeta)./((Lambdap+mu-Lambdap*zeta)*(sigma+eta+mu));
Cp = binput*Lambdap*p*sigma*(1-zeta)./((eta+theta0+mu)*(Lambdap+mu-Lambdap*zeta)*(eta+sigma+mu));
S  = binput*(1-p)./(Lambdap+mu);
P  = binput*p./(Lambdap+mu-Lambdap*zeta);

N=binput/mu;
%disp((I+C+Ip+Cp)/(N-(S+P)));
%disp((N-(S+I+C+Ip+Cp+P))/N);
Ptot = 1 - mu/binput*(S+P);
Pun  = (I+C+Ip+Cp)*mu/binput;

%alpha_temp solution de Pun==0
PrevThresh=0;
alpha = alphaSICTP(betaI,betaC,theta0,gamma0,sigma,zeta,eta,p,mu,binput,PrevThresh);
end