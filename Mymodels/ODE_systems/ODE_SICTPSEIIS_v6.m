function [dY] = ODE_SICTPSEIIS_v7(t,Y,betaIh,betaCh,thetah,sigmah,zetah,ph,...
                                         betaCt, gammaCt, nuCt, epsCt, sigmaCt,...
                                         rho_h,rho_c,rho_hc,...
                                         eta_h_prep,eta_c_prep,eta_c_art,...
                                         VTunderART,mu,b)
dY=zeros(28,1);
N=b/mu;
LambdaCt = betaCt*sum(Y(15:28))/N;
Lambdah  = betaIh*sum(Y([2:7:23,5:7:26]))/N + betaCh*sum(Y([3:7:24,6:7:27]))/N;
rho_cg=0;

dY(1) = Y(8)*(rho_c + rho_cg + rho_hc) - Y(1)*(Lambdah + LambdaCt + mu) - b*(ph - 1) + Y(15)*(nuCt + rho_c + rho_cg + rho_hc) + Y(22)*(gammaCt + nuCt);
dY(2) = Lambdah*Y(1) - Y(2)*(LambdaCt + mu + rho_h + rho_hc + sigmah) + Y(16)*(nuCt + rho_c + rho_cg) + Y(23)*(gammaCt + nuCt) + Y(9)*(rho_c + rho_cg);
dY(3) = Y(2)*sigmah - Y(3)*(LambdaCt + mu + rho_h + rho_hc + thetah) + Y(17)*(nuCt + rho_c + rho_cg) + Y(24)*(gammaCt + nuCt) + Y(10)*(rho_c + rho_cg);
dY(4) = Y(11)*eta_c_prep + b*ph + Y(25)*(eta_c_prep + gammaCt + nuCt) - Y(4)*(LambdaCt + mu - Lambdah*(zetah - 1)) + Y(18)*(eta_c_prep + nuCt);
dY(5) = Y(12)*eta_c_prep + Y(26)*(eta_c_prep + gammaCt + nuCt) - Y(5)*(LambdaCt + eta_h_prep + mu + sigmah) + Y(19)*(eta_c_prep + nuCt) - Lambdah*Y(4)*(zetah - 1);
dY(6) = Y(13)*eta_c_prep + Y(5)*sigmah + Y(27)*(eta_c_prep + gammaCt + nuCt) - Y(6)*(LambdaCt + eta_h_prep + mu + thetah) + Y(20)*(eta_c_prep + nuCt);
dY(7) = Y(5)*eta_h_prep + Y(9)*rho_hc + Y(10)*rho_hc + Y(16)*rho_hc + Y(17)*rho_hc + Y(28)*(eta_c_art + gammaCt + nuCt) + Y(3)*(rho_h + rho_hc + thetah) + Y(14)*(eta_c_art + rho_c + rho_cg + VTunderART*rho_hc) - Y(7)*(LambdaCt + mu) + Y(6)*(eta_h_prep + thetah) + Y(21)*(eta_c_art + nuCt + rho_c + rho_cg + VTunderART*rho_hc) + Y(2)*(rho_h + rho_hc);
dY(8) = LambdaCt*Y(1) - Y(8)*(Lambdah + mu + rho_c + rho_cg + rho_hc + sigmaCt);
dY(9) = Lambdah*Y(8) + LambdaCt*Y(2) - Y(9)*(mu + rho_c + rho_h + rho_cg + rho_hc + sigmah + sigmaCt);
dY(10) = LambdaCt*Y(3) + Y(9)*sigmah - Y(10)*(mu + rho_c + rho_h + rho_cg + rho_hc + sigmaCt + thetah);
dY(11) = LambdaCt*Y(4) - Y(11)*(eta_c_prep + mu + sigmaCt - Lambdah*(zetah - 1));
dY(12) = LambdaCt*Y(5) - Y(12)*(eta_c_prep + eta_h_prep + mu + sigmah + sigmaCt) - Lambdah*Y(11)*(zetah - 1);
dY(13) = LambdaCt*Y(6) + Y(12)*sigmah - Y(13)*(eta_c_prep + eta_h_prep + mu + sigmaCt + thetah);
dY(14) = LambdaCt*Y(7) + Y(12)*eta_h_prep + Y(9)*rho_h - Y(14)*(eta_c_art + mu + rho_c + rho_cg + sigmaCt + VTunderART*rho_hc) + Y(13)*(eta_h_prep + thetah) + Y(10)*(rho_h + thetah);
dY(15) = - Y(15)*(Lambdah + mu + nuCt + rho_c + rho_cg + rho_hc) - Y(8)*sigmaCt*(epsCt - 1);
dY(16) = Lambdah*Y(15) - Y(16)*(mu + nuCt + rho_c + rho_h + rho_cg + rho_hc + sigmah) - Y(9)*sigmaCt*(epsCt - 1);
dY(17) = Y(16)*sigmah - Y(17)*(mu + nuCt + rho_c + rho_h + rho_cg + rho_hc + thetah) - Y(10)*sigmaCt*(epsCt - 1);
dY(18) = - Y(18)*(eta_c_prep + mu + nuCt - Lambdah*(zetah - 1)) - Y(11)*sigmaCt*(epsCt - 1);
dY(19) = - Y(19)*(eta_c_prep + eta_h_prep + mu + nuCt + sigmah) - Lambdah*Y(18)*(zetah - 1) - Y(12)*sigmaCt*(epsCt - 1);
dY(20) = Y(19)*sigmah - Y(20)*(eta_c_prep + eta_h_prep + mu + nuCt + thetah) - Y(13)*sigmaCt*(epsCt - 1);
dY(21) = Y(19)*eta_h_prep + Y(16)*rho_h - Y(21)*(eta_c_art + mu + nuCt + rho_c + rho_cg + VTunderART*rho_hc) + Y(20)*(eta_h_prep + thetah) + Y(17)*(rho_h + thetah) - Y(14)*sigmaCt*(epsCt - 1);
dY(22) = Y(8)*epsCt*sigmaCt - Y(22)*(Lambdah + gammaCt + mu + nuCt);
dY(23) = Lambdah*Y(22) - Y(23)*(gammaCt + mu + nuCt + sigmah) + Y(9)*epsCt*sigmaCt;
dY(24) = Y(23)*sigmah - Y(24)*(gammaCt + mu + nuCt + thetah) + Y(10)*epsCt*sigmaCt;
dY(25) = Y(11)*epsCt*sigmaCt - Y(25)*(eta_c_prep + gammaCt + mu + nuCt - Lambdah*(zetah - 1));
dY(26) = Y(12)*epsCt*sigmaCt - Y(26)*(eta_c_prep + eta_h_prep + gammaCt + mu + nuCt + sigmah) - Lambdah*Y(25)*(zetah - 1);
dY(27) = Y(26)*sigmah - Y(27)*(eta_c_prep + eta_h_prep + gammaCt + mu + nuCt + thetah) + Y(13)*epsCt*sigmaCt;
dY(28) = Y(26)*eta_h_prep + Y(24)*thetah - Y(28)*(eta_c_art + gammaCt + mu + nuCt) + Y(27)*(eta_h_prep + thetah) + Y(14)*epsCt*sigmaCt;


end

