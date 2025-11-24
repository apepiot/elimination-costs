function [dY] = ODE_SEIIS2_v4(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                 beta2,nu2,eps2,sigma2,gamma2,mu,b,rho)

    [~,Lambda1] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,rho);
    [~,Lambda2] = Rp_SEIIS_v4(beta2,nu2,eps2,sigma2,gamma2,mu,b,rho);

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

