function [dY] = SISxSIS_V2(t,Y,b,beta1, beta2,e1,e2,gamma1,gamma2, s, delta, mu, type)
% Version 2 of the SISxSIS
% what's changed ? there is an additional voluntary testing rate, from I12
% to S et from I1 to S and from I2 to S
    
    S=Y(1);I1=Y(2);I2=Y(3);I12=Y(4);
    
    if strcmp(type,'frequency')
        N = sum(Y);
    end
    if strcmp(type,'density')
        N = 1;
    end

    lambda1 = (beta1*I1+e1*I12)/N;
    lambda2 = (beta2*I2+e2*I12)/N;
    
    dS      = b - (lambda1 + lambda2)*S + (gamma1 + s*delta)*I1 + (gamma2 + s*delta)*I2 + s*delta*I12 - mu*S; %S
    dI1     = lambda1*S - lambda2*I1 + gamma2*I12 - (gamma1 + s*delta)*I1 - mu*I1;
    dI2     = lambda2*S - lambda1*I2 + gamma1*I12 - (gamma2 + s*delta)*I2 - mu*I2;
    dI12    = lambda2*I1 + lambda1*I2 - (gamma1 + gamma2 + s*delta)*I12 - mu*I12;
    
    dY = [dS; dI1; dI2; dI12];   
    
    
 
    
end

