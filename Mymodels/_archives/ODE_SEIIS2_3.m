function [dY] = ODE_SEIIS2_3(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                 beta2,nu2,eps2,sigma2,gamma2,mu,b,rho)

    %    N = sum(Y);
    %     totInf1 = sum(Y(indexInf1));
    %     Lambda1   = beta1*totInf1/N;   
    %     totInf2 = sum(Y(indexInf2));
    %     Lambda2   = beta2*totInf2/N;
    
    R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rho))./((mu + sigma1 + rho).*(gamma1 + mu + nu1).*(mu + nu1 + rho)); 
    Lambda1 = max(beta1*(R1-1)*(sigma1+rho+mu)/(beta1+(sigma1+rho+mu)*R1),0);
    R2 = (beta2*sigma2*(gamma2*(1-eps2) + mu + nu2 + eps2*rho))./((mu + sigma2 + rho).*(gamma2 + mu + nu2).*(mu + nu2 + rho)); 
    Lambda2 = max(beta2*(R2-1)*(sigma2+rho+mu)/(beta2+(sigma2+rho+mu)*R2),0);

    dY = [b - Y(1)*(Lambda1 + Lambda2 + mu) + Y(2)*rho + Y(5)*rho + Y(6)*rho + Y(7)*rho + Y(10)*rho + Y(11)*rho + Y(4)*(gamma1 + nu1) + Y(13)*(gamma2 + nu2) + Y(3)*(nu1 + rho) + Y(9)*(nu2 + rho);...
     Lambda1*Y(1) + Y(10)*nu2 - Y(2)*(Lambda2 + mu + rho + sigma1) + Y(14)*(gamma2 + nu2);...%E1
     Y(11)*nu2 - Y(3)*(Lambda2 + mu + nu1 + rho) + Y(15)*(gamma2 + nu2) - Y(2)*sigma1*(eps1 - 1);...%IA1
     Y(12)*nu2 - Y(4)*(Lambda2 + gamma1 + mu + nu1) + Y(16)*(gamma2 + nu2) + Y(2)*eps1*sigma1;...%IS1
     Lambda2*Y(1) + Y(7)*nu1 - Y(5)*(Lambda1 + mu + rho + sigma2) + Y(8)*(gamma1 + nu1);...%E2
     Lambda2*Y(2) + Lambda1*Y(5) - Y(6)*(mu + rho + sigma1 + sigma2);...%E1E2
     Lambda2*Y(3) - Y(7)*(mu + nu1 + rho + sigma2) - Y(6)*sigma1*(eps1 - 1);...%IA1E2
     Lambda2*Y(4) - Y(8)*(gamma1 + mu + nu1 + sigma2) + Y(6)*eps1*sigma1;...%IS1E2
     Y(11)*nu1 - Y(9)*(Lambda1 + mu + nu2 + rho) + Y(12)*(gamma1 + nu1) - Y(5)*sigma2*(eps2 - 1);...%IA2
     Lambda1*Y(9) - Y(10)*(mu + nu2 + rho + sigma1) - Y(6)*sigma2*(eps2 - 1);...%E1IA2
     - Y(11)*(mu + nu1 + nu2 + rho) - Y(7)*sigma2*(eps2 - 1) - Y(10)*sigma1*(eps1 - 1);...%IA1IA2
     Y(10)*eps1*sigma1 - Y(12)*(gamma1 + mu + nu1 + nu2) - Y(8)*sigma2*(eps2 - 1);...%IS1IA2
     Y(15)*nu1 - Y(13)*(Lambda1 + gamma2 + mu + nu2) + Y(16)*(gamma1 + nu1) + Y(5)*eps2*sigma2;...%IS2
     Lambda1*Y(13) - Y(14)*(gamma2 + mu + nu2 + sigma1) + Y(6)*eps2*sigma2;...%E1IS2
     Y(7)*eps2*sigma2 - Y(15)*(gamma2 + mu + nu1 + nu2) - Y(14)*sigma1*(eps1 - 1);...%IA1IS2
     Y(8)*eps2*sigma2 - Y(16)*(gamma1 + gamma2 + mu + nu1 + nu2) + Y(14)*eps1*sigma1];%IS1IS2
end

