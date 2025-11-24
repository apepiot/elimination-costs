function [dY] = SIRSIS2(t,Y,b,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,rho,mu,type)
% SIS3 (a terminer)
    S=Y(1);  I1=Y(2); I2=Y(3); I3=Y(4);
    I12=Y(5);I23=Y(6);I13=Y(7);I123=Y(8);
    R1=Y(9); IR2=Y(10); IR3=Y(11); IR23=Y(12);
    if strcmp(type,'frequency')
        N = sum(Y);
    end
    if strcmp(type,'density')
        N = 1;
    end

    lambda1 = beta1*(I1+I12+I13+I123)/N;
    lambda2 = beta2*(I2+I12+I23+I123+IR2+IR23)/N;
    lambda3 = beta3*(I3+I13+I23+I123+IR3+IR23)/N;
    
    gamma1p = gamma1+s1*rho;
    gamma2p = gamma2+s2*rho;
    gamma3p = gamma3+s3*rho;
    gamma12t = s1*s2*rho;
    gamma13t = s1*s3*rho;
    gamma23t = s2*s3*rho;
    gamma123t = s1*s2*s3*rho;
    gamma12 = gamma1p+gamma2p-gamma12t;
    gamma23 = gamma2p+gamma3p-gamma23t;
    gamma13 = gamma1p+gamma3p-gamma13t;
    gamma123 = gamma1+gamma2+gamma3+(1-(1-s1)*(1-s2)*(1-s3))*rho;
    gamma1t2 = gamma1p-s1*s2*rho;
    gamma1t3 = gamma1p-s1*s3*rho;
    gamma2t1 = gamma2p-s1*s2*rho;
    gamma2t3 = gamma2p-s3*s2*rho;
    gamma3t1 = gamma3p-s1*s3*rho;
    gamma3t2 = gamma3p-s3*s2*rho;
    gamma1t  = gamma1p-gamma12t-gamma13t+gamma123t;%I123->23
    gamma2t  = gamma2p-gamma12t-gamma23t+gamma123t;%I123->13
    gamma3t  = gamma3p-gamma13t-gamma23t+gamma123t;%I123->12
    
   %from S,                                 I1,                         I2,                         I3,                         I12,                I23,                I13,                I123,                   R1,                     IR2,                    IR3,                    IR23
    MAT = [-(lambda1+lambda2+lambda3+mu),   0,                          gamma2p,                    gamma3p,                    0,                  gamma23t,           0,                  0,                      0,                      0,                      0,                      0;... 
        lambda1,                            -lambda2-lambda3-gamma1p-mu, 0,                         0,                          gamma2t1,           0,                  gamma3t1,           gamma23t-gamma123t,     0,                      0,                      0,                      0;... 
        lambda2,                            0,                          -lambda1-lambda3-gamma2p-mu,0,                          0,                  gamma3t2,           0,                  0,                      0,                      0,                      0,                      0;... 
        lambda3,                            0,                          0,                          -lambda1-lambda2-gamma3p-mu,0,                  gamma2t3,           0,                  0,                      0,                      0,                      0,                      0;... 
        0,                                  lambda2,                    lambda1,                    0,                          -lambda3-gamma12-mu,0,                  0,                  gamma3t,                0,                      0,                      0,                      0;... 
        0,                                  0,                          lambda3,                    lambda2,                    0,                  -lambda1-gamma23-mu,0,                  0,                      0,                      0,                      0,                      0;... 
        0,                                  lambda3,                    0,                          lambda1,                    0,                  0,                  -lambda2-gamma13-mu,gamma2t,                0,                      0,                      0,                      0;... 
        0,                                  0,                          0,                          0,                          lambda3,            lambda1,            lambda2,            -(gamma123+mu),         0,                      0,                      0,                      0;...      
        0,                                  gamma1p,                    0,                          0,                          gamma12t,           0,                  gamma13t,           gamma123t,              -(lambda2+lambda3+mu),  gamma2p,                gamma3p,                gamma23t;...
        0,                                  0,                          0,                          0,                          gamma1t2,           0,                  0,                  gamma13t-gamma123t,     lambda2,                -(gamma2p+lambda3+mu),  0,                      gamma3t2;...
        0,                                  0,                          0,                          0,                          0,                  0,                  gamma1t3,           gamma12t-gamma123t,     lambda3,                0,                      -(gamma3p+lambda2+mu),  gamma2t3;...
        0,                                  0,                          0,                          0,                          0,                  0,                  0,                  gamma1t,                0,                      lambda3,                lambda2,                -(gamma23+mu);...
        ];
    
    dY = [b;0;0;0;0;0;0;0;0;0;0;0] + MAT*Y;   
 
end

