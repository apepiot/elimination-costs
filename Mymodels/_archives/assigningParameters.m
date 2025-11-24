
% Reassign data - all instances
% HIV data
betaIh = ampl.getParameter('betaIh'); betaIh.setValues([paramHIV.betaI]);
betaCh = ampl.getParameter('betaCh'); betaCh.setValues([paramHIV.betaC]);
sigmah = ampl.getParameter('sigmah'); sigmah.setValues([paramHIV.sigma]);
thetah = ampl.getParameter('thetah'); thetah.setValues([paramHIV.theta]);
zetah = ampl.getParameter('zetah'); zetah.setValues([paramHIV.zeta]);
eta_h_prep = ampl.getParameter('eta_h_prep'); eta_h_prep.setValues([paramHIV.eta]);
ph = ampl.getParameter('ph'); ph.setValues([paramHIV.p]);

% syphilis data
betas = ampl.getParameter('betas'); betas.setValues([paramS.beta]);
sigmas = ampl.getParameter('sigmas'); sigmas.setValues([paramS.sigma]);
gamma3s = ampl.getParameter('gamma3s'); gamma3s.setValues([paramS.gamma3]);
taus = ampl.getParameter('taus'); taus.setValues([paramS.tau]);
thetas = ampl.getParameter('thetas'); thetas.setValues([paramS.theta]);

% Ct data
betaX = ampl.getParameter('betaX'); betaX.setValues([paramC.beta]);
gammaX = ampl.getParameter('gammaX'); gammaX.setValues([paramC.gamma]);
nuX = ampl.getParameter('nuX'); nuX.setValues([paramC.nu]);
epsX = ampl.getParameter('epsX'); epsX.setValues([paramC.eps]);
sigmaX = ampl.getParameter('sigmaX'); sigmaX.setValues([paramC.sigma]);

%Ng
betaY = ampl.getParameter('betaY'); betaY.setValues([paramG.beta]);
gammaY = ampl.getParameter('gammaY'); gammaY.setValues([paramG.gamma]);
nuY = ampl.getParameter('nuY'); nuY.setValues([paramG.nu]);
epsY = ampl.getParameter('epsY'); epsY.setValues([paramG.eps]);
sigmaY = ampl.getParameter('sigmaY'); sigmaY.setValues([paramG.sigma]);

%General parameters
mu_apml = ampl.getParameter('mu'); mu_apml.setValues([mu]);
b_ampl = ampl.getParameter('b'); b_ampl.setValues([b]);

%Testing parameters (baseline testing rate)
rho_h = ampl.getParameter('rho_h'); rho_h.setValues([paramRho.rho_h]);
rho_s = ampl.getParameter('rho_s'); rho_s.setValues([paramRho.rho_s]);
rho_c = ampl.getParameter('rho_c'); rho_c.setValues([paramRho.rho_c]);
rho_g = ampl.getParameter('rho_g'); rho_g.setValues([paramRho.rho_g]);

%Testing rates (kit)
rho_hs   = ampl.getParameter('rho_hs'); rho_hs.setValues([paramRho.rho_hs]);
rho_hc   = ampl.getParameter('rho_hc'); rho_hc.setValues([paramRho.rho_hc]);
rho_hg   = ampl.getParameter('rho_hg'); rho_hg.setValues([paramRho.rho_hg]);
rho_sc   = ampl.getParameter('rho_sc'); rho_sc.setValues([paramRho.rho_sc]);
rho_sg   = ampl.getParameter('rho_sg'); rho_sg.setValues([paramRho.rho_sg]);
rho_cg   = ampl.getParameter('rho_cg'); rho_cg.setValues([paramRho.rho_cg]);
rho_hsc  = ampl.getParameter('rho_hsc'); rho_hsc.setValues([paramRho.rho_hsc]);
rho_hsg  = ampl.getParameter('rho_hsg'); rho_hsg.setValues([paramRho.rho_hsg]);
rho_hcg  = ampl.getParameter('rho_hcg'); rho_hcg.setValues([paramRho.rho_hcg]);
rho_scg  = ampl.getParameter('rho_scg'); rho_scg.setValues([paramRho.rho_scg]);
rho_hscg = ampl.getParameter('rho_hscg'); rho_hscg.setValues([paramRho.rho_hscg]);

eta_s_prep = ampl.getParameter('eta_s_prep'); eta_s_prep.setValues([paramRho.eta_s_prep]);
eta_c_prep = ampl.getParameter('eta_c_prep'); eta_c_prep.setValues([paramRho.eta_c_prep]);
eta_g_prep = ampl.getParameter('eta_g_prep'); eta_g_prep.setValues([paramRho.eta_g_prep]);
eta_s_art  = ampl.getParameter('eta_s_art'); eta_s_art.setValues([paramRho.eta_s_art]);
eta_c_art  = ampl.getParameter('eta_c_art'); eta_c_art.setValues([paramRho.eta_c_art]);
eta_g_art  = ampl.getParameter('eta_g_art'); eta_g_art.setValues([paramRho.eta_g_art]);
