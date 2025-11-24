function [Rp,Lambdap,alpha] = Rp_SICR_v4(betaI,betaC,theta,sigma,gamma_sicr,mu,b,rho)
%attention a la definition de theta et gamma_sicr. 
%dans la plupart des cas, ici theta<-gamma et gamma_sicr=0;

%%Case lambda = (betaI*I+betaC*C)/N;
Rp = (betaI.*(theta+rho+mu)+betaC*sigma)./((sigma+gamma_sicr+rho+mu).*(theta+rho+mu)); %OK (voir page 53 du recap)
%IHIV = b*(1-1./Rp)./(sigma+gamma_sicr+rho+mu);
%CHIV = b*sigma*(1-1./Rp)./((theta+rho+mu).*(sigma+gamma_sicr+rho+mu));
%Lambdap = max((betaI*IHIV+betaC*CHIV)/(b/mu),0);
Lambdap = mu*(Rp-1);
alpha1 = betaI/2 - gamma_sicr/2 - mu - sigma/2 - theta/2 - (betaI^2 - 2*betaI*gamma_sicr - 2*betaI*sigma + 2*betaI*theta + gamma_sicr^2 + 2*gamma_sicr*sigma - 2*gamma_sicr*theta + sigma^2 - 2*sigma*theta + 4*betaC*sigma + theta^2)^(1/2)/2;
alpha2 = betaI/2 - gamma_sicr/2 - mu - sigma/2 - theta/2 + (betaI^2 - 2*betaI*gamma_sicr - 2*betaI*sigma + 2*betaI*theta + gamma_sicr^2 + 2*gamma_sicr*sigma - 2*gamma_sicr*theta + sigma^2 - 2*sigma*theta + 4*betaC*sigma + theta^2)^(1/2)/2;
alpha = max(alpha1,alpha2);
end

