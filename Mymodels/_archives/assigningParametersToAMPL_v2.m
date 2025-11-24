function ampl = assigningParametersToAMPL_v2(paramTab,paramRho,mu,b,ampl,k,mod)

if contains(mod,'h')
    paramHIV = getParams('h',paramTab);
    betaIh = ampl.getParameter('betaIh'); betaIh.setValues([paramHIV.betaI]);
    betaCh = ampl.getParameter('betaCh'); betaCh.setValues([paramHIV.betaC]);
    sigmah = ampl.getParameter('sigmah'); sigmah.setValues([paramHIV.sigma]);
    thetah = ampl.getParameter('thetah'); thetah.setValues([paramHIV.theta0]);
    zetah = ampl.getParameter('zetah');   zetah.setValues([paramHIV.zeta]);
    ph = ampl.getParameter('ph');         ph.setValues([paramHIV.p]);
    eta_h_prep = ampl.getParameter('eta_h_prep');  eta_h_prep.setValues([paramRho.eta_h_prep]);
    rho_h = ampl.getParameter('rho_h'); rho_h.setValues([paramRho.rho_h]);
end
if contains(mod,'s')
    paramS = getParams('s',paramTab);
    warning('to do')
end
if contains(mod,'c')
    paramC = getParams('c',paramTab);
    betaX = ampl.getParameter('betaX');   betaX.setValues([paramC.beta]);
    gammaX = ampl.getParameter('gammaX'); gammaX.setValues([paramC.gamma]);
    nuX = ampl.getParameter('nuX');       nuX.setValues([paramC.nu]);
    epsX = ampl.getParameter('epsX');     epsX.setValues([paramC.eps]);
    sigmaX = ampl.getParameter('sigmaX'); sigmaX.setValues([paramC.sigma]);
    eta_X_prep = ampl.getParameter('eta_X_prep');     eta_X_prep.setValues([paramRho.eta_c_prep]);
    eta_X_art  = ampl.getParameter('eta_X_art');      eta_X_art.setValues([paramRho.eta_c_art]);
    rho_X = ampl.getParameter('rho_X'); rho_X.setValues([paramRho.rho_c]);
end

if contains(mod,'g') && ~contains(k,'c')
    paramG = getParams('g',paramTab);
    betaX = ampl.getParameter('betaX');   betaX.setValues([paramG.beta]);
    gammaX = ampl.getParameter('gammaX'); gammaX.setValues([paramG.gamma]);
    nuX = ampl.getParameter('nuX');       nuX.setValues([paramG.nu]);
    epsX = ampl.getParameter('epsX');     epsX.setValues([paramG.eps]);
    sigmaX = ampl.getParameter('sigmaX'); sigmaX.setValues([paramG.sigma]);
    eta_X_prep = ampl.getParameter('eta_X_prep');     eta_X_prep.setValues([paramRho.eta_g_prep]);
    eta_X_art  = ampl.getParameter('eta_X_art');      eta_X_art.setValues([paramRho.eta_g_art]);
    rho_X = ampl.getParameter('rho_X'); rho_X.setValues([paramRho.rho_g]);
elseif contains(mod,'g') && contains(k,'c')
    paramG = getParams('g',paramTab);
    betaY = ampl.getParameter('betaY');   betaY.setValues([paramG.beta]);
    gammaY = ampl.getParameter('gammaY'); gammaY.setValues([paramG.gamma]);
    nuY = ampl.getParameter('nuY');       nuY.setValues([paramG.nu]);
    epsY = ampl.getParameter('epsY');     epsY.setValues([paramG.eps]);
    sigmaY = ampl.getParameter('sigmaY'); sigmaY.setValues([paramG.sigma]);
    eta_Y_prep = ampl.getParameter('eta_Y_prep');     eta_Y_prep.setValues([paramRho.eta_g_prep]);
    eta_Y_art  = ampl.getParameter('eta_Y_art');      eta_Y_art.setValues([paramRho.eta_g_art]);
    rho_Y = ampl.getParameter('rho_Y'); rho_Y.setValues([paramRho.rho_g]);
end

if length(mod)>=2
    if contains(mod,'h') && contains(mod,'c') && ~strcmp(k,'hc')
        rho_hX   = ampl.getParameter('rho_hX'); rho_hX.setValues([paramRho.rho_hc]);
    end
    if contains(mod,'h') && contains(mod,'g') && ~contains(mod,'c') && ~strcmp(k,'hg')
        rho_hX   = ampl.getParameter('rho_hX'); rho_hX.setValues([paramRho.rho_hg]);
    end
    if contains(mod,'h') && contains(mod,'g') && contains(mod,'c') && ~strcmp(k,'hg')
        rho_hY   = ampl.getParameter('rho_hY'); rho_hY.setValues([paramRho.rho_hg]);
    end
    if contains(mod,'h') && contains(mod,'s') && ~strcmp(k,'hs')
        rho_hs   = ampl.getParameter('rho_hs'); rho_hs.setValues([paramRho.rho_hs]);
    end
    if contains(mod,'s') && contains(mod,'c') && ~strcmp(k,'sc')
        rho_sX   = ampl.getParameter('rho_sX'); rho_sX.setValues([paramRho.rho_sc]);
    end
    if contains(mod,'s') && contains(mod,'g') && ~contains(mod,'c') && ~strcmp(k,'sg')
        rho_sX   = ampl.getParameter('rho_sX'); rho_sX.setValues([paramRho.rho_sg]);
    end
    if contains(mod,'s') && contains(mod,'g') && contains(mod,'c') && ~strcmp(k,'sg')
        rho_sY   = ampl.getParameter('rho_sY'); rho_sY.setValues([paramRho.rho_sg]);
    end
    
    if contains(mod,'c') && contains(mod,'g')
        rho_XY   = ampl.getParameter('rho_XY'); rho_XY.setValues([paramRho.rho_cg]);
    end
end

if length(k)>=3
    if (strcmp(k,'hcg') || strcmp(k,'hscg')) && ~strcmp(k,'hcg')
        rho_hXY   = ampl.getParameter('rho_hXY'); rho_hXY.setValues([paramRho.rho_hcg]);
    end
    if (strcmp(k,'hsc') || strcmp(k,'hscg')) && ~strcmp(k,'hsc')
        rho_hsX   = ampl.getParameter('rho_hsX'); rho_hsX.setValues([paramRho.rho_hsc]);
    end
    if strcmp(k,'hsg') && ~strcmp(k,'hsg')
        rho_hsX   = ampl.getParameter('rho_hsX'); rho_hsX.setValues([paramRho.rho_hsg]);
    end
    if strcmp(k,'hscg') && ~strcmp(k,'hsg')
        rho_hsY   = ampl.getParameter('rho_hsY'); rho_hsY.setValues([paramRho.rho_hsg]);
    end
    if strcmp(k,'scg') && ~strcmp(k,'scg')
        rho_sXY   = ampl.getParameter('rho_sXY'); rho_sXY.setValues([paramRho.rho_scg]);
    end
end

if length(k)>=4 && ~strcmp(k,'hscg')
    rho_hsXY   = ampl.getParameter('rho_hsXY'); rho_hsXY.setValues([paramRho.rho_hscg]);
end

%General parameters
mu_apml = ampl.getParameter('mu'); mu_apml.setValues([mu]);
b_ampl = ampl.getParameter('b');   b_ampl.setValues([b]);

%Routine testing under PrEP/ART
vt_under_art  = ampl.getParameter('VTunderART');  vt_under_art.setValues([paramRho.VTunderART]);

end