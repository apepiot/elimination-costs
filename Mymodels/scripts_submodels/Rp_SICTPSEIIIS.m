function [eigvalues] = Rp_SICTPSEIIIS(betaIh,betaCh,thetah,sigmah,zetah,ph,...
                                      betas, sigmas, taus, thetas, gamma3s,...
                                      rho_h,rho_s,rho_hs,...
                                      eta_h_prep,eta_s_prep,eta_s_art,...
                                      VTunderART,mu,b)

syms Y [35,1]
syms F [35,1];
syms V [35,1];

N=b/mu;
Lambdas = betas*sum(Y(8:35))/N;
Lambdah = betaIh*sum(Y([2:7:30,5:7:33]))/N + betaCh*sum(Y([3:7:31,6:7:34]))/N;

F(1)= 0;
F(2) = Y1*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b);
F(3) = 0;
F(4) = 0;
F(5) = Y4*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b);
F(6) = 0;
F(7) = 0;
F(8) = (Y1*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
F(9) = Y8*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) + (Y2*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
F(10) =(Y3*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
F(11) = (Y4*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
F(12) = Y11*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) + (Y5*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
F(13) = (Y6*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
F(14) = (Y7*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
F(15) = 0;
F(16) = Y15*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b);
F(17) = 0;
F(18) = 0;
F(19) = Y18*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b);
F(20) = 0;
F(21) = 0;
F(22) = 0;
F(23) = Y22*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b);
F(24) = 0;
F(25) = 0;
F(26) = Y25*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b);
F(27) = 0;
F(28) = 0;
F(29) = 0;
F(30) = Y29*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b);
F(31) = 0;
F(32) = 0;
F(33) = Y32*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b);
F(34) = 0;
F(35) = 0;



V(1) = Y9*rho_hs - Y1*(Lambdah + Lambdas + mu) + Y10*rho_hs + Y16*rho_hs + Y17*rho_hs + Y23*rho_hs + Y24*rho_hs + Y30*rho_hs + Y31*rho_hs + Y29*(gamma3s + rho_s + rho_hs) + Y2*(rho_h + rho_hs) + Y3*(rho_h + rho_hs) + Y8*(rho_s + rho_hs) + Y15*(rho_s + rho_hs) + Y22*(rho_s + rho_hs);
V(2) = Lambdah*Y1 - Y2*(Lambdas + mu + rho_h + rho_hs + sigmah) + Y9*rho_s + Y16*rho_s + Y23*rho_s - Y1*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) + Y30*(gamma3s + rho_s);
V(3) = Y10*rho_s + Y17*rho_s + Y24*rho_s - Y3*(Lambdas + mu + rho_h + rho_hs + thetah) + Y2*sigmah + Y31*(gamma3s + rho_s);
V(4) = Y32*gamma3s - Y4*(Lambdas + mu - Lambdah*(zetah - 1));
V(5) = Y33*gamma3s - Y5*(Lambdas + eta_h_prep + mu + sigmah) - Y4*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) - Lambdah*Y4*(zetah - 1);
V(6) = Y34*gamma3s + Y5*sigmah - Y6*(Lambdas + eta_h_prep + mu + thetah);
V(7) = Y35*(gamma3s + rho_s + VTunderART*rho_hs) + Y5*eta_h_prep + Y3*thetah + Y14*(rho_s + VTunderART*rho_hs) + Y21*(rho_s + VTunderART*rho_hs) + Y28*(rho_s + VTunderART*rho_hs) - Y7*(Lambdas + mu) + Y6*(eta_h_prep + thetah);
V(8) = Lambdas*Y1 - Y8*(Lambdah + mu + rho_s + rho_hs + sigmas) + Y9*rho_h + Y10*rho_h - (Y1*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
V(9) = Lambdah*Y8 + Lambdas*Y2 - Y8*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) - Y9*(mu + rho_h + rho_s + rho_hs + sigmah + sigmas) - (Y2*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
V(10) = Lambdas*Y3 + Y9*sigmah - Y10*(mu + rho_h + rho_s + rho_hs + sigmas + thetah) - (Y3*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
V(11) = Lambdas*Y4 - Y11*(eta_s_prep + mu + sigmas - Lambdah*(zetah - 1)) - (Y4*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
V(12) = Lambdas*Y5 - Y12*(eta_h_prep + eta_s_prep + mu + sigmah + sigmas) - Y11*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) - Lambdah*Y11*(zetah - 1) - (Y5*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
V(13) = Lambdas*Y6 + Y12*sigmah - Y13*(eta_h_prep + eta_s_prep + mu + sigmas + thetah) - (Y6*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
V(14) = Lambdas*Y7 + Y12*eta_h_prep + Y10*thetah + Y13*(eta_h_prep + thetah) - Y14*(eta_s_art + mu + rho_s + sigmas + VTunderART*rho_hs) - (Y7*betas*mu*(Y8 + Y9 + Y10 + Y11 + Y12 + Y13 + Y14 + Y15 + Y16 + Y17 + Y18 + Y19 + Y20 + Y21 + Y22 + Y23 + Y24 + Y25 + Y26 + Y27 + Y28 + Y29 + Y30 + Y31 + Y32 + Y33 + Y34 + Y35))/b;
V(15) = Y16*rho_h + Y17*rho_h - Y15*(Lambdah + mu + rho_s + rho_hs + taus) + Y8*sigmas;
V(16) = Lambdah*Y15 + Y9*sigmas - Y15*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) - Y16*(mu + rho_h + rho_s + rho_hs + sigmah + taus);
V(17) = Y16*sigmah + Y10*sigmas - Y17*(mu + rho_h + rho_s + rho_hs + taus + thetah);
V(18) = Y11*sigmas - Y18*(eta_s_prep + mu + taus - Lambdah*(zetah - 1));
V(19) = Y12*sigmas - Y19*(eta_h_prep + eta_s_prep + mu + sigmah + taus) - Y18*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) - Lambdah*Y18*(zetah - 1);
V(20) = Y19*sigmah + Y13*sigmas - Y20*(eta_h_prep + eta_s_prep + mu + taus + thetah);
V(21) = Y19*eta_h_prep + Y14*sigmas + Y17*thetah + Y20*(eta_h_prep + thetah) - Y21*(eta_s_art + mu + rho_s + taus + VTunderART*rho_hs);
V(22) = Y23*rho_h + Y24*rho_h - Y22*(Lambdah + mu + rho_s + rho_hs + thetas) + Y15*taus;
V(23) = Lambdah*Y22 + Y16*taus - Y22*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) - Y23*(mu + rho_h + rho_s + rho_hs + sigmah + thetas);
V(24) = Y23*sigmah + Y17*taus - Y24*(mu + rho_h + rho_s + rho_hs + thetah + thetas);
V(25) = Y18*taus - Y25*(eta_s_prep + mu + thetas - Lambdah*(zetah - 1));
V(26) = Y19*taus - Y26*(eta_h_prep + eta_s_prep + mu + sigmah + thetas) - Y25*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) - Lambdah*Y25*(zetah - 1);
V(27) = Y26*sigmah + Y20*taus - Y27*(eta_h_prep + eta_s_prep + mu + thetah + thetas);
V(28) = Y26*eta_h_prep + Y21*taus + Y24*thetah + Y27*(eta_h_prep + thetah) - Y28*(eta_s_art + mu + rho_s + thetas + VTunderART*rho_hs);
V(29) = Y30*rho_h - Y29*(Lambdah + gamma3s + mu + rho_s + rho_hs) + Y31*rho_h + Y22*thetas;
V(30) = Lambdah*Y29 + Y23*thetas - Y29*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) - Y30*(gamma3s + mu + rho_h + rho_s + rho_hs + sigmah);
V(31) = Y30*sigmah + Y24*thetas - Y31*(gamma3s + mu + rho_h + rho_s + rho_hs + thetah);
V(32) = Y25*thetas - Y32*(eta_s_prep + gamma3s + mu - Lambdah*(zetah - 1));
V(33) = Y26*thetas - Y33*(eta_h_prep + eta_s_prep + gamma3s + mu + sigmah) - Y32*((betaIh*mu*(Y2 + Y5 + Y9 + Y12 + Y16 + Y19 + Y23 + Y26 + Y30 + Y33))/b + (betaCh*mu*(Y3 + Y6 + Y10 + Y13 + Y17 + Y20 + Y24 + Y27 + Y31 + Y34))/b) - Lambdah*Y32*(zetah - 1);
V(34) = Y33*sigmah - Y34*(eta_h_prep + eta_s_prep + gamma3s + mu + thetah) + Y27*thetas;
V(35) = Y33*eta_h_prep + Y31*thetah + Y28*thetas - Y35*(eta_s_art + gamma3s + mu + rho_s + VTunderART*rho_hs) + Y34*(eta_h_prep + thetah);

syms dF [35,35];
syms dV [35,35];
for i=1:35
    for j=1:35
        dF(i,j) = diff(F(i),Y(j));
        dV(i,j) = diff(V(i),Y(j));
    end
end



disp('calcul de V^-1')
Vmoins1 = dV^(-1);
disp('calcul des eigenvalues de FVmoins1')
eigvalues = eig(F*Vmoins1);

end

