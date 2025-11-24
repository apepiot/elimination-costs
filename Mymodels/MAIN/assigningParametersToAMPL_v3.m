function ampl = assigningParametersToAMPL_v3(paramTab,paramRho,mu,b,ampl,kit)
paramHIV = getParams('h',paramTab);
paramS = getParams('s',paramTab);
paramC = getParams('c',paramTab);
paramG = getParams('g',paramTab);

betaIh = ampl.getParameter('betaIh'); betaIh.setValues([paramHIV.betaI]);
betaCh = ampl.getParameter('betaCh'); betaCh.setValues([paramHIV.betaC]);
sigmah = ampl.getParameter('sigmah'); sigmah.setValues([paramHIV.sigma]);
thetah = ampl.getParameter('thetah'); thetah.setValues([paramHIV.theta0]);
zetah = ampl.getParameter('zetah'); zetah.setValues([paramHIV.zeta]);
eta_h_prep = ampl.getParameter('eta_h_prep'); eta_h_prep.setValues([paramRho.eta_h_prep]);
ph = ampl.getParameter('ph_mod'); ph.setValues([paramHIV.p]);

% syphilis data
betas = ampl.getParameter('betas'); betas.setValues([paramS.beta]);
sigmas = ampl.getParameter('sigmas'); sigmas.setValues([paramS.sigma]);
gamma3s = ampl.getParameter('gamma3s'); gamma3s.setValues([paramS.gamma30]);
taus = ampl.getParameter('taus'); taus.setValues([paramS.tau]);
thetas = ampl.getParameter('thetas'); thetas.setValues([paramS.theta]);
eta_s_prep = ampl.getParameter('eta_s_prep'); eta_s_prep.setValues([paramRho.eta_s_prep]);
eta_s_art = ampl.getParameter('eta_s_art'); eta_s_art.setValues([paramRho.eta_s_art]);

% Ct data
betaX = ampl.getParameter('beta_c'); betaX.setValues([paramC.beta]);
gammaX = ampl.getParameter('gamma_c'); gammaX.setValues([paramC.gamma]);
nuX = ampl.getParameter('nu_c'); nuX.setValues([paramC.nu]);
epsX = ampl.getParameter('eps_c'); epsX.setValues([paramC.eps]);
sigmaX = ampl.getParameter('sigma_c'); sigmaX.setValues([paramC.sigma]);
eta_c_prep = ampl.getParameter('eta_c_prep'); eta_c_prep.setValues([paramRho.eta_c_prep]);
eta_c_art = ampl.getParameter('eta_c_art'); eta_c_art.setValues([paramRho.eta_c_art]);

%Ng
betaY = ampl.getParameter('beta_g'); betaY.setValues([paramG.beta]);
gammaY = ampl.getParameter('gamma_g'); gammaY.setValues([paramG.gamma]);
nuY = ampl.getParameter('nu_g'); nuY.setValues([paramG.nu]);
epsY = ampl.getParameter('eps_g'); epsY.setValues([paramG.eps]);
sigmaY = ampl.getParameter('sigma_g'); sigmaY.setValues([paramG.sigma]);
eta_g_prep = ampl.getParameter('eta_g_prep'); eta_g_prep.setValues([paramRho.eta_g_prep]);
eta_g_art = ampl.getParameter('eta_g_art'); eta_g_art.setValues([paramRho.eta_g_art]);


%General parameters
mu_apml = ampl.getParameter('mu'); mu_apml.setValues([mu]);
b_ampl = ampl.getParameter('b');   b_ampl.setValues([b]);

%Routine testing under PrEP/ART
vt_under_art  = ampl.getParameter('VTunderART');  vt_under_art.setValues([paramRho.VTunderART]);

%Testing parameters (baseline testing rate)
if ~isequal(kit,{'HIV'}) && ~isequal(kit,'h')
    rho_h = ampl.getParameter('rho_h'); rho_h.setValues([paramRho.rho_h]);
end
if ~isequal(kit,{'syphilis'}) && ~isequal(kit,'s')
    rho_s = ampl.getParameter('rho_s'); rho_s.setValues([paramRho.rho_s]);
end
if ~isequal(kit,{'Ct'}) && ~isequal(kit,'c')
    rho_c = ampl.getParameter('rho_c'); rho_c.setValues([paramRho.rho_c]);
end
if ~isequal(kit,{'Ng'}) && ~isequal(kit,'g')
    rho_g = ampl.getParameter('rho_g'); rho_g.setValues([paramRho.rho_g]);
end

%Testing rates (kit)
if ~isequal(kit,{'HIV','syphilis'}) && ~isequal(kit,'hs')
    rho_hs   = ampl.getParameter('rho_hs'); rho_hs.setValues([paramRho.rho_hs]);
end
if ~isequal(kit,{'HIV','Ct'}) && ~isequal(kit,'hc')
    rho_hc   = ampl.getParameter('rho_hc'); rho_hc.setValues([paramRho.rho_hc]);
end
if ~isequal(kit,{'HIV','Ng'}) && ~isequal(kit,'hg')
    rho_hg   = ampl.getParameter('rho_hg'); rho_hg.setValues([paramRho.rho_hg]);
end
if ~isequal(kit,{'syphilis','Ct'}) && ~isequal(kit,'sc')
    rho_sc   = ampl.getParameter('rho_sc'); rho_sc.setValues([paramRho.rho_sc]);
end
if ~isequal(kit,{'syphilis','Ng'}) && ~isequal(kit,'sg')
    rho_sg   = ampl.getParameter('rho_sg'); rho_sg.setValues([paramRho.rho_sg]);
end
if ~isequal(kit,{'Ct','Ng'}) && ~isequal(kit,'cg')
    rho_cg   = ampl.getParameter('rho_cg'); rho_cg.setValues([paramRho.rho_cg]);
end
if ~isequal(kit,{'HIV','syphilis','Ct'}) && ~isequal(kit,'hsc')
    rho_hsc  = ampl.getParameter('rho_hsc'); rho_hsc.setValues([paramRho.rho_hsc]);
end
if ~isequal(kit,{'HIV','syphilis','Ng'}) && ~isequal(kit,'hsg')
    rho_hsg  = ampl.getParameter('rho_hsg'); rho_hsg.setValues([paramRho.rho_hsg]);
end
if ~isequal(kit,{'HIV','Ct','Ng'}) && ~isequal(kit,'hcg')
    rho_hcg  = ampl.getParameter('rho_hcg'); rho_hcg.setValues([paramRho.rho_hcg]);
end
if ~isequal(kit,{'syphilis','Ct','Ng'}) && ~isequal(kit,'scg')
    rho_scg  = ampl.getParameter('rho_scg'); rho_scg.setValues([paramRho.rho_scg]);
end
if ~isequal(kit,{'HIV','syphilis','Ct','Ng'}) && ~isequal(kit,'hscg')
    rho_hscg = ampl.getParameter('rho_hscg'); rho_hscg.setValues([paramRho.rho_hscg]);
end

end