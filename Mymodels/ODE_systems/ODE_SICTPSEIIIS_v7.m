function [dY] = ODE_SICTPSEIIIS_v7(t,Y,betaIh,betaCh,thetah,sigmah,zetah,ph,...
                                        betas,sigmas,gamma3s,taus,thetas,...
                                        rho_h,rho_s,rho_hs,...
                                        eta_h_prep,eta_s_prep,eta_s_art,...
                                        VTunderART,mu,b)
                                    
dY = zeros(35,1);
%syms dY ;
N = b/mu;
Lambdah  = betaIh*sum(Y([2:7:30,5:7:33]))/N + betaCh*sum(Y([3:7:31,6:7:34]))/N;
Lambdas  = betas*sum(Y([15:35]))/N;

%syms Lambdah
%syms Lambdas

dY(1) = Y(29)*(gamma3s + rho_s + rho_hs) - Y(1)*(Lambdah + Lambdas + mu) - b*(ph - 1) + Y(8)*(rho_s + rho_hs) + Y(15)*(rho_s + rho_hs) + Y(22)*(rho_s + rho_hs);
dY(2) = Lambdah*Y(1) - Y(2)*(Lambdas + mu + rho_h + rho_hs + sigmah) + Y(9)*rho_s + Y(16)*rho_s + Y(23)*rho_s + Y(30)*(gamma3s + rho_s);
dY(3) = Y(10)*rho_s + Y(17)*rho_s + Y(24)*rho_s - Y(3)*(Lambdas + mu + rho_h + rho_hs + thetah) + Y(2)*sigmah + Y(31)*(gamma3s + rho_s);
dY(4) = Y(11)*eta_s_prep + Y(18)*eta_s_prep + Y(25)*eta_s_prep + b*ph - Y(4)*(Lambdas + mu - Lambdah*(zetah - 1)) + Y(32)*(eta_s_prep + gamma3s);
dY(5) = Y(12)*eta_s_prep + Y(19)*eta_s_prep + Y(26)*eta_s_prep - Y(5)*(Lambdas + eta_h_prep + mu + sigmah) + Y(33)*(eta_s_prep + gamma3s) - Lambdah*Y(4)*(zetah - 1);
dY(6) = Y(13)*eta_s_prep + Y(20)*eta_s_prep + Y(27)*eta_s_prep + Y(5)*sigmah - Y(6)*(Lambdas + eta_h_prep + mu + thetah) + Y(34)*(eta_s_prep + gamma3s);
dY(7) = Y(14)*(eta_s_art + rho_s + VTunderART*rho_hs) + Y(21)*(eta_s_art + rho_s + VTunderART*rho_hs) + Y(28)*(eta_s_art + rho_s + VTunderART*rho_hs) + Y(5)*eta_h_prep + Y(9)*rho_hs + Y(10)*rho_hs + Y(16)*rho_hs + Y(17)*rho_hs + Y(23)*rho_hs + Y(24)*rho_hs + Y(30)*rho_hs + Y(31)*rho_hs + Y(3)*(rho_h + rho_hs + thetah) + Y(35)*(eta_s_art + gamma3s + rho_s + VTunderART*rho_hs) - Y(7)*(Lambdas + mu) + Y(6)*(eta_h_prep + thetah) + Y(2)*(rho_h + rho_hs);
dY(8) = Lambdas*Y(1) - Y(8)*(Lambdah + mu + rho_s + rho_hs + sigmas);
dY(9) = Lambdah*Y(8) + Lambdas*Y(2) - Y(9)*(mu + rho_h + rho_s + rho_hs + sigmah + sigmas);
dY(10) = Lambdas*Y(3) + Y(9)*sigmah - Y(10)*(mu + rho_h + rho_s + rho_hs + sigmas + thetah);
dY(11) = Lambdas*Y(4) - Y(11)*(eta_s_prep + mu + sigmas - Lambdah*(zetah - 1));
dY(12) = Lambdas*Y(5) - Y(12)*(eta_h_prep + eta_s_prep + mu + sigmah + sigmas) - Lambdah*Y(11)*(zetah - 1);
dY(13) = Lambdas*Y(6) + Y(12)*sigmah - Y(13)*(eta_h_prep + eta_s_prep + mu + sigmas + thetah);
dY(14) = Lambdas*Y(7) + Y(12)*eta_h_prep + Y(9)*rho_h + Y(13)*(eta_h_prep + thetah) - Y(14)*(eta_s_art + mu + rho_s + sigmas + VTunderART*rho_hs) + Y(10)*(rho_h + thetah);
dY(15) = Y(8)*sigmas - Y(15)*(Lambdah + mu + rho_s + rho_hs + taus);
dY(16) = Lambdah*Y(15) + Y(9)*sigmas - Y(16)*(mu + rho_h + rho_s + rho_hs + sigmah + taus);
dY(17) = Y(16)*sigmah + Y(10)*sigmas - Y(17)*(mu + rho_h + rho_s + rho_hs + taus + thetah);
dY(18) = Y(11)*sigmas - Y(18)*(eta_s_prep + mu + taus - Lambdah*(zetah - 1));
dY(19) = Y(12)*sigmas - Y(19)*(eta_h_prep + eta_s_prep + mu + sigmah + taus) - Lambdah*Y(18)*(zetah - 1);
dY(20) = Y(19)*sigmah + Y(13)*sigmas - Y(20)*(eta_h_prep + eta_s_prep + mu + taus + thetah);
dY(21) = Y(19)*eta_h_prep + Y(16)*rho_h + Y(14)*sigmas + Y(20)*(eta_h_prep + thetah) - Y(21)*(eta_s_art + mu + rho_s + taus + VTunderART*rho_hs) + Y(17)*(rho_h + thetah);
dY(22) = Y(15)*taus - Y(22)*(Lambdah + mu + rho_s + rho_hs + thetas);
dY(23) = Lambdah*Y(22) + Y(16)*taus - Y(23)*(mu + rho_h + rho_s + rho_hs + sigmah + thetas);
dY(24) = Y(23)*sigmah + Y(17)*taus - Y(24)*(mu + rho_h + rho_s + rho_hs + thetah + thetas);
dY(25) = Y(18)*taus - Y(25)*(eta_s_prep + mu + thetas - Lambdah*(zetah - 1));
dY(26) = Y(19)*taus - Y(26)*(eta_h_prep + eta_s_prep + mu + sigmah + thetas) - Lambdah*Y(25)*(zetah - 1);
dY(27) = Y(26)*sigmah + Y(20)*taus - Y(27)*(eta_h_prep + eta_s_prep + mu + thetah + thetas);
dY(28) = Y(26)*eta_h_prep + Y(23)*rho_h + Y(21)*taus + Y(27)*(eta_h_prep + thetah) - Y(28)*(eta_s_art + mu + rho_s + thetas + VTunderART*rho_hs) + Y(24)*(rho_h + thetah);
dY(29) = Y(22)*thetas - Y(29)*(Lambdah + gamma3s + mu + rho_s + rho_hs);
dY(30) = Lambdah*Y(29) + Y(23)*thetas - Y(30)*(gamma3s + mu + rho_h + rho_s + rho_hs + sigmah);
dY(31) = Y(30)*sigmah + Y(24)*thetas - Y(31)*(gamma3s + mu + rho_h + rho_s + rho_hs + thetah);
dY(32) = Y(25)*thetas - Y(32)*(eta_s_prep + gamma3s + mu - Lambdah*(zetah - 1));
dY(33) = Y(26)*thetas - Y(33)*(eta_h_prep + eta_s_prep + gamma3s + mu + sigmah) - Lambdah*Y(32)*(zetah - 1);
dY(34) = Y(33)*sigmah - Y(34)*(eta_h_prep + eta_s_prep + gamma3s + mu + thetah) + Y(27)*thetas;
dY(35) = Y(33)*eta_h_prep + Y(30)*rho_h + Y(28)*thetas - Y(35)*(eta_s_art + gamma3s + mu + rho_s + VTunderART*rho_hs) + Y(34)*(eta_h_prep + thetah) + Y(31)*(rho_h + thetah);

end

