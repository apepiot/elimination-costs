function [dY] = SIIS2(t,Y,b,beta1,beta2,gamma10,gamma20,rho,eps1,eps2,sigma1,sigma2,nu1,nu2,mu,type)
% SIISxSIIS
% what's changed ? there is an additional voluntary testing rate, from I12
% to S et from I1 to S and from I2 to S
    
    S=Y(1);IA0=Y(2);IS0=Y(3);
    I0A=Y(4);IAA=Y(5);ISA=Y(6);
    I0S=Y(7);IAS=Y(8);ISS=Y(9); 
    
    if strcmp(type,'frequency')
        N = sum(Y);
    end
    if strcmp(type,'density')
        N = 1;
    end

    Lambda1 = beta1*(N-S-I0A-I0S)/N;
    Lambda2 = beta2*(N-S-IA0-IS0)/N;
    
    dS   = b - (Lambda1 + Lambda2)*S + (1-eps1)*nu1*IA0 + (1-eps2)*nu2*I0A +...
        (nu1+gamma10)*IS0 + (nu2+gamma20)*I0S + rho*(N-S) - eps1*rho*IA0 - eps2*rho*I0A - mu*S;
    dIA0 = Lambda1*S - Lambda2*IA0 - ((1-eps1)*(nu1+rho) + eps1*sigma1 + mu)*IA0 + (1-eps2)*nu2*IAA + (nu2+gamma20)*IAS;
    dI0A = Lambda2*S - Lambda1*I0A - ((1-eps2)*(nu2+rho) + eps2*sigma2 + mu)*I0A + (nu1+gamma10)*ISA +(1-eps1)*nu1*IAA;
    dIS0 = eps1*sigma1*IA0 - Lambda2*IS0 + (1-eps2)*nu2*ISA -(nu1+gamma10+rho + mu)*IS0 + (nu2+gamma20)*ISS;
    dI0S = eps2*sigma2*I0A - Lambda1*I0S + (nu1+gamma10)*ISS + (1-eps1)*nu1*IAS - (nu2+gamma20+rho + mu)*I0S;
    dISA = eps1*sigma1*IAA + Lambda2*IS0 - ((1-eps2)*nu2 + nu1+gamma10 + eps2*sigma2 + rho + mu)*ISA;
    dIAS = Lambda1*I0S + eps2*sigma2*IAA - ((1-eps1)*nu1 + nu2+gamma20 + eps1*sigma1 + rho + mu)*IAS;
    dIAA = Lambda1*I0A + Lambda2*IA0 - ((1-eps1)*nu1 + eps1*sigma1 + (1-eps2)*nu2 + eps2*sigma2 + rho +mu)*IAA;
    dISS = eps1*sigma1*IAS + eps2*sigma2*ISA - (nu1+gamma10+rho + nu2+gamma20 + mu)*ISS;

    dY = [dS; dIA0 ; dIS0 ; dI0A ; dIAA ; dISA ; dI0S ; dIAS ; dISS];
 
end

