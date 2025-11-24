function [paramTab,mu,vecAlphas] = sampleParameters_v3_extent(sampleCt,sampleNg,sampleHIV,sampleS,b,pHIV)

[paramTab,mu,vecAlphas] = sampleParameters_v3(sampleCt,sampleNg,sampleHIV,sampleS,b);   %Ct,Ng,HIV,syph  
paramTab{3}.modelType='SICTP';
paramTab{3}.eta  = 4;
paramTab{3}.zeta = round(randPERT(46,60,71,1)/100,5);
paramTab{3}.alpha_p0 = paramTab{3}.alpha_prev;
paramTab{3}.p = pHIV;

[paramTab{3}.R_prep_base,~,~,paramTab{3}.Ptot_prep_base,paramTab{3}.Pun_prep_base] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
    paramTab{3}.theta0,paramTab{3}.gamma0,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
    paramTab{3}.p,paramTab{3}.mu,b,paramTab{3}.rhob);

[paramTab{3}.R_prep_0,~,paramTab{3}.alpha_prep,paramTab{3}.Ptot_prep_0,paramTab{3}.Pun_prep_0] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
    paramTab{3}.theta0,paramTab{3}.gamma0,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
    paramTab{3}.p,paramTab{3}.mu,b,0);

vecAlphas(3) = paramTab{3}.alpha_prep;
end

