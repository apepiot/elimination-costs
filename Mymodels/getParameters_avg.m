addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\INITIALISATION\parametersSets')
% Lecture des paramètres
ID_ech = 4000000;
mu=1/30.6;b=100000;

%% Ct
P1VT = 0.053; paramTab{1}.P_base=P1VT;
R1VT    = 1/(1-P1VT);         
paramTab{1}.R_base = R1VT;
paramTab{1}.disease='Ct';
eps1 = 0.11; paramTab{1}.eps = eps1;
gamma1 = 365/14; paramTab{1}.gamma = gamma1;
paramTab{1}.mini_d = 'c';
paramTab{1}.modelType = 'SEIIS';
paramTab{1}.mu = mu;
nu1 = 12/12; paramTab{1}.nu = nu1;
paramTab{1}.pi = b;
rhob1 = 0.12; paramTab{1}.rhob = rhob1;
sigma1 = 365/11; paramTab{1}.sigma = sigma1;
paramTab{1}.IDech_id = ID_ech;
beta1   = R1VT/((sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rhob1))./((mu + sigma1 + rhob1).*(gamma1 + mu + nu1).*(mu + nu1 + rhob1)));
paramTab{1}.beta = beta1;
[R10,~,alpha1] = Rp_SEIIS_v4(beta1,nu1,eps1,sigma1,gamma1,mu,b,0);
paramTab{1}.alpha = alpha1;
paramTab{1}.R0 = R10;
P10     = 1-1/R10;
paramTab{1}.P0 = P10;

%% Ng
P2VT = 0.047; paramTab{2}.P_base=P2VT;
R2VT    = 1/(1-P2VT);         
paramTab{2}.R_base = R2VT;
paramTab{2}.disease='Ng';
eps2 = 0.8; paramTab{2}.eps = eps2;
gamma2 = 365/14; paramTab{2}.gamma = gamma2;
paramTab{2}.mini_d = 'g';
paramTab{2}.modelType = 'SEIIS';
paramTab{2}.mu = mu;
nu2 = 12/6; paramTab{2}.nu = nu2;
paramTab{2}.pi = b;
rhob2=0.12; paramTab{2}.rhob = rhob2;
sigma2 = 365/5; paramTab{2}.sigma = sigma2;
paramTab{2}.IDech_id = ID_ech;
beta2   = R2VT/((sigma2*(gamma2*(1-eps2) + mu + nu2 + eps2*rhob2))./((mu + sigma2 + rhob2).*(gamma2 + mu + nu2).*(mu + nu2 + rhob2)));
paramTab{2}.beta = beta2;
[R20,~,alpha2] = Rp_SEIIS_v4(beta2,nu2,eps2,sigma2,gamma2,mu,b,0);
paramTab{2}.alpha = alpha2;
paramTab{2}.R0 = R20;
P20     = 1-1/R20;
paramTab{2}.P0 = P20;

%% syphilis
PS_VT     = 6.6/100;
RS_VT     = 1/(1-PS_VT);
rhobs     = 1/2.48;
sigmaS    = 365/25;
tauS      = 365/46;
thetaS    = 12/3.6;
gamma1S0  = 0;
gamma1S   = gamma1S0+rhobs;
gamma3S0  = 1/20; 
gamma3S   = gamma3S0+rhobs;
nuS       = 0;
betaS = RS_VT/((sigmaS.*((gamma3S+mu+nuS).*(mu+rhobs+tauS+thetaS)+tauS*thetaS))./((mu+rhobs+sigmaS).*(mu+rhobs+thetaS).*(gamma3S+mu+nuS).*(gamma1S+mu+tauS)));
[R0S,~,alphaS] = Rp_SEIIIS_v4(betaS,sigmaS,tauS,nuS,gamma1S0,thetaS,gamma3S0,mu,b,0);
P0S = 1-1/R0S;

paramTab{4}.beta = betaS; paramTab{4}.sigma=sigmaS; paramTab{4}.tau=tauS; 
paramTab{4}.gamma1_base=gamma1S; paramTab{4}.gamma10=gamma1S0;
paramTab{4}.gamma3_base=gamma3S; paramTab{4}.gamma30=gamma3S0;
paramTab{4}.theta=thetaS;
paramTab{4}.nu=nuS; paramTab{4}.alpha=alphaS;
paramTab{4}.rhob = rhobs;
paramTab{4}.R_base = RS_VT; paramTab{4}.R0 = R0S;
paramTab{4}.disease = 'syphilis'; paramTab{4}.modelType = 'SEIIIS'; paramTab{4}.mini_d = 's';
paramTab{4}.P_base = PS_VT; paramTab{4}.P0 = P0S; 
paramTab{4}.mu = mu; paramTab{4}.pi=b;
paramTab{4}.IDech_id = ID_ech;

%% HIV 
%prevalence and not undiagnosed prevalence, with voluntary testing
PHIV_prevagay_VT = 16.1/100;
R_sict_VT   = 1/(1-PHIV_prevagay_VT);           %P=1-1/R so R=1/(1-P);          %check
sigma_sict  = 52/8.2; 
rhobh       = 1/1.9;                                        %check
theta_sict0 = 1/4.4;           %1/9.8;                         %check, theta_sict0: basé sur les symptomes
theta_sict  = theta_sict0+rhobh;
gamma_sict0 = 0;                                %I->T                           %check
gamma_sict  = rhobh+gamma_sict0;
ratioBeta   = 9.1;           %betaI/betaC                    %check
betaC_sict  = R_sict_VT*(sigma_sict+gamma_sict+mu)*(theta_sict+mu)/(ratioBeta*(theta_sict+mu)+sigma_sict); %check
betaI_sict  = ratioBeta*betaC_sict; %check
[R0_sict,~,alpha_sict] = Rp_SICR_v4(betaI_sict,betaC_sict,theta_sict0,sigma_sict,gamma_sict0,mu,b,0); %check
P0_sict = 1-1/R0_sict;

%undiagnosed prevalence :
Pund_base = mu*(theta_sict+sigma_sict+mu)/(betaI_sict*(theta_sict+mu)+betaC_sict*sigma_sict)*(R_sict_VT-1);
Pund_0    = mu*(theta_sict0+sigma_sict+mu)/(betaI_sict*(theta_sict0+mu)+betaC_sict*sigma_sict)*(R0_sict-1);

paramTab{3}.betaI = betaI_sict; paramTab{3}.betaC = betaC_sict; 
paramTab{3}.gamma0 = gamma_sict0; paramTab{3}.theta0 = theta_sict0;
paramTab{3}.gamma_base = gamma_sict; paramTab{3}.theta_base = theta_sict;
paramTab{3}.sigma=sigma_sict; paramTab{3}.alpha_prev = alpha_sict;
paramTab{3}.rhob = rhobh;
paramTab{3}.R_prev_base = R_sict_VT; paramTab{3}.R_prev_0 = R0_sict; 
paramTab{3}.Pund_prev_base = Pund_base; paramTab{3}.Pund_prev_0 = Pund_0;
paramTab{3}.Ptot_prev_base = PHIV_prevagay_VT; paramTab{3}.Ptot_prev_0 = P0_sict;
paramTab{3}.disease = 'HIV'; paramTab{3}.modelType = 'SICT'; paramTab{3}.mini_d = 'h';
paramTab{3}.mu = mu; paramTab{3}.pi=b;
paramTab{3}.IDech_id = ID_ech;

% adding
paramTab{3}.modelType='SICTP';
paramTab{3}.eta  = 4;
paramTab{3}.zeta = 60/100;
paramTab{3}.alpha_p0 = paramTab{3}.alpha_prev;
paramTab{3}.p = pHIV;

[paramTab{3}.R_prep_base,~,~,paramTab{3}.Ptot_prep_base,paramTab{3}.Pun_prep_base] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
    paramTab{3}.theta0,paramTab{3}.gamma0,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
    paramTab{3}.p,paramTab{3}.mu,b,paramTab{3}.rhob);

[paramTab{3}.R_prep_0,~,paramTab{3}.alpha_prep,paramTab{3}.Ptot_prep_0,paramTab{3}.Pun_prep_0] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
    paramTab{3}.theta0,paramTab{3}.gamma0,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
    paramTab{3}.p,paramTab{3}.mu,b,0);




%%
verbose=0;
paramSolver.tolP0 = 0.5e-4;
paramSolver.maxBndAlpha=20;
paramSolver.nbRelanceMax=5;
paramSolver.timeLimit = 20; %seconds
paramSolver.iterMaxDicho = 30;
paramSolver.tolAlpha = 1e-4;
paramSolver.method_alpha = 'dicho';
paramSolver.timeSolver = 20;
%opt.TolP0=0.5e-4;

pathRes = ['.\ParameterAnalysis\results_',num2str(paramNo),'\_round_',num2str(roundNo),'\'];
paramRho = table2struct(readtable([pathRes ,'paramRho.txt']));

paramRho.eta_s_prep=4;
%paramRho.eta_c_prep=0;

%paramRho.VTunderART=0;
