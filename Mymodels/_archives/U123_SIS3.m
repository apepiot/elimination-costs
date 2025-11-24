function [U,P123] = U123_SIS3(rho,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c)
    % Function U123 of testing for the SISxSISxSIS model
    gamma1p = gamma1+s1*rho;
    gamma2p = gamma2+s2*rho;
    gamma3p = gamma3+s3*rho;
    R1p = beta1./(gamma1p+mu);
    R2p = beta2./(gamma2p+mu);
    R3p = beta3./(gamma3p+mu);
    N = b/mu;
    lambda1=beta1.*(1-1./R1p);
    lambda2=beta2.*(1-1./R2p);
    lambda3=beta3.*(1-1./R3p);
    gamma12t = s1*s2*rho;
    gamma13t = s1*s3*rho;
    gamma23t = s2*s3*rho;
    gamma1t2 = gamma1p-s1*s2*rho;
    gamma1t3 = gamma1p-s1*s3*rho;
    gamma2t1 = gamma2p-s1*s2*rho;
    gamma2t3 = gamma2p-s3*s2*rho;
    gamma3t1 = gamma3p-s1*s3*rho;
    gamma3t2 = gamma3p-s3*s2*rho;
    
    %I123
    S13   = b/mu*(gamma1t3./R3p+gamma3t1./R1p+mu+gamma13t)./(beta1+beta3-mu-gamma13t); %T(1)+T(3)
    I13_2 = b/mu*(1-1./R1p-1./R3p)+S13; %I13+I123 (prevalence des coinfectes 1 et 3 I13, voir SIS^2)
    S23   = b/mu*(gamma2t3./R3p+gamma3t2./R2p+mu+gamma23t)./(beta2+beta3-mu-gamma23t); %T(1)+T(2)
    I23_1 = b/mu*(1-1./R2p-1./R3p)+S23; %I23+I123, T(6)+T(8)
    S12   = b/mu*(gamma1t2./R2p+gamma2t1./R1p+mu+gamma12t)./(beta1+beta2-mu-gamma12t); %T(1)+T(4)
    I12_3 = b/mu*(1-1./R1p-1./R2p)+S12; %I12+I123 T(5)+T(8)
    gamma123 = gamma1+gamma2+gamma3+(1-(1-s1)*(1-s2)*(1-s3)).*rho;
    
    I123th = (lambda1.*I23_1+lambda2.*I13_2+lambda3.*I12_3)./(lambda1+lambda2+lambda3+gamma123+mu);
    
    %S123
    Sth = S12+S23+S13+b/mu.*(1-1./R1p-1./R2p-1./R3p)-I123th;
    P123 = 1-Sth./N;

    %utility
    U = rho.*(P123-c);
end

