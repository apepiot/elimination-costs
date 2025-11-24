function [F,V] = FVevaluatedInY()
syms Y [35,1]
syms mu b;
N=b/mu;
syms betaCh betaIh betas
syms LambdaCt; syms nuCt; syms epsCt; syms sigmaCt; syms gammaCt;
syms LambdaNg; syms nuNg; syms epsNg; syms sigmaNg; syms gammaNg;
syms betaIh; syms betaCh; 

syms Lambdah thetah sigmah ph eta_h_prep zetah
syms Lambdas sigmas taus thetas gamma1s gamma3s nus
% targeted testing rate:
syms rho_h rho_s rho_hs
syms eta_s_prep
syms eta_h_prep
syms eta_s_art
syms VTunderART
Lambdah = betaIh*sum(Y([2:7:30,5:7:33]))/N + betaCh*sum(Y([3:7:31,6:7:34]))/N;
Lambdas = betas*sum(Y(15:35))/N;
F = [Lambdah*Y1, 0, Lambdah*Y4, 0, 0, Lambdas*Y1, Lambdah*Y8 + Lambdas*Y2, Lambdas*Y3, Lambdas*Y4, Lambdah*Y11 + Lambdas*Y5, Lambdas*Y6, Lambdas*Y7, 0, Lambdah*Y15, 0, 0, Lambdah*Y18, 0, 0, 0, Lambdah*Y22, 0, 0, Lambdah*Y25, 0, 0, 0, Lambdah*Y29, 0, 0, Lambdah*Y32, 0, 0, 0, 0];
V = [Y2*(Lambdas + mu + rho_h + rho_hs + sigmah) - Y9*rho_s - Y16*rho_s - Y23*rho_s - Y30*(gamma3s + rho_s), Y3*(Lambdas + mu + rho_h + rho_hs + thetah) - Y17*rho_s - Y24*rho_s - Y10*rho_s - Y2*sigmah - Y31*(gamma3s + rho_s), Lambdah*Y4 - Y12*eta_s_prep - Y19*eta_s_prep - Y26*eta_s_prep + Y5*(Lambdas + eta_h_prep + mu + sigmah) - Y33*(eta_s_prep + gamma3s) + Lambdah*Y4*(zetah - 1), Y6*(Lambdas + eta_h_prep + mu + thetah) - Y20*eta_s_prep - Y27*eta_s_prep - Y5*sigmah - Y13*eta_s_prep - Y34*(eta_s_prep + gamma3s), Y7*(Lambdas + mu) - Y21*(eta_s_art + rho_s + VTunderART*rho_hs) - Y28*(eta_s_art + rho_s + VTunderART*rho_hs) - Y5*eta_h_prep - Y9*rho_hs - Y10*rho_hs - Y16*rho_hs - Y17*rho_hs - Y23*rho_hs - Y24*rho_hs - Y30*rho_hs - Y31*rho_hs - Y3*(rho_h + rho_hs + thetah) - Y35*(eta_s_art + gamma3s + rho_s + VTunderART*rho_hs) - Y14*(eta_s_art + rho_s + VTunderART*rho_hs) - Y6*(eta_h_prep + thetah) - Y2*(rho_h + rho_hs), Y8*(Lambdah + mu + rho_s + rho_hs + sigmas), Y9*(mu + rho_h + rho_s + rho_hs + sigmah + sigmas), Y10*(mu + rho_h + rho_s + rho_hs + sigmas + thetah) - Y9*sigmah, Y11*(eta_s_prep + mu + sigmas - Lambdah*(zetah - 1)), Lambdah*Y11 + Y12*(eta_h_prep + eta_s_prep + mu + sigmah + sigmas) + Lambdah*Y11*(zetah - 1), Y13*(eta_h_prep + eta_s_prep + mu + sigmas + thetah) - Y12*sigmah, Y14*(eta_s_art + mu + rho_s + sigmas + VTunderART*rho_hs) - Y9*rho_h - Y13*(eta_h_prep + thetah) - Y12*eta_h_prep - Y10*(rho_h + thetah), Y15*(Lambdah + mu + rho_s + rho_hs + taus) - Y8*sigmas, Y16*(mu + rho_h + rho_s + rho_hs + sigmah + taus) - Y9*sigmas, Y17*(mu + rho_h + rho_s + rho_hs + taus + thetah) - Y10*sigmas - Y16*sigmah, Y18*(eta_s_prep + mu + taus - Lambdah*(zetah - 1)) - Y11*sigmas, Lambdah*Y18 - Y12*sigmas + Y19*(eta_h_prep + eta_s_prep + mu + sigmah + taus) + Lambdah*Y18*(zetah - 1), Y20*(eta_h_prep + eta_s_prep + mu + taus + thetah) - Y13*sigmas - Y19*sigmah, Y21*(eta_s_art + mu + rho_s + taus + VTunderART*rho_hs) - Y16*rho_h - Y14*sigmas - Y20*(eta_h_prep + thetah) - Y19*eta_h_prep - Y17*(rho_h + thetah), Y22*(Lambdah + mu + rho_s + rho_hs + thetas) - Y15*taus, Y23*(mu + rho_h + rho_s + rho_hs + sigmah + thetas) - Y16*taus, Y24*(mu + rho_h + rho_s + rho_hs + thetah + thetas) - Y17*taus - Y23*sigmah, Y25*(eta_s_prep + mu + thetas - Lambdah*(zetah - 1)) - Y18*taus, Lambdah*Y25 - Y19*taus + Y26*(eta_h_prep + eta_s_prep + mu + sigmah + thetas) + Lambdah*Y25*(zetah - 1), Y27*(eta_h_prep + eta_s_prep + mu + thetah + thetas) - Y20*taus - Y26*sigmah, Y28*(eta_s_art + mu + rho_s + thetas + VTunderART*rho_hs) - Y23*rho_h - Y21*taus - Y27*(eta_h_prep + thetah) - Y26*eta_h_prep - Y24*(rho_h + thetah), Y29*(Lambdah + gamma3s + mu + rho_s + rho_hs) - Y22*thetas, Y30*(gamma3s + mu + rho_h + rho_s + rho_hs + sigmah) - Y23*thetas, Y31*(gamma3s + mu + rho_h + rho_s + rho_hs + thetah) - Y24*thetas - Y30*sigmah, Y32*(eta_s_prep + gamma3s + mu - Lambdah*(zetah - 1)) - Y25*thetas, Lambdah*Y32 + Y33*(eta_h_prep + eta_s_prep + gamma3s + mu + sigmah) - Y26*thetas + Lambdah*Y32*(zetah - 1), Y34*(eta_h_prep + eta_s_prep + gamma3s + mu + thetah) - Y33*sigmah - Y27*thetas, Y35*(eta_s_art + gamma3s + mu + rho_s + VTunderART*rho_hs) - Y30*rho_h - Y28*thetas - Y33*eta_h_prep - Y34*(eta_h_prep + thetah) - Y31*(rho_h + thetah), Y1*(Lambdah + Lambdas + mu) - Y29*(gamma3s + rho_s + rho_hs) - Y8*(rho_s + rho_hs) - Y15*(rho_s + rho_hs) - Y22*(rho_s + rho_hs), Y4*(Lambdas + mu - Lambdah*(zetah - 1)) - Y18*eta_s_prep - Y25*eta_s_prep - Y11*eta_s_prep - Y32*(eta_s_prep + gamma3s)];

end

