%% Test du systeme d'ODE (avec ode45)
%clear all;
close all;

%% 0. Parameters
%clear all;
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\')
b=100;

newSet=0; defaultSet=0;
if newSet
    pHIV = 0.3;
    [paramTab,mu,~] = sampleParameters_v3_extent(true,true,true,true,b,pHIV);   %Ct,Ng,HIV,syph
elseif defaultSet && ~newSet
elseif ~defaultSet && ~newSet
    mu=1/35;
    paramTab{1}.beta = 1.3759;
    paramTab{1}.gamma = 24.8922;
    paramTab{1}.nu = 0.5246;
    paramTab{1}.eps = 0.0486;
    paramTab{1}.sigma=40.2979;
    [~,~,paramTab{1}.alpha] = Rp_SEIIS_v4(paramTab{1}.beta,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,paramTab{1}.gamma,mu,b,0);
    
    paramTab{2}.beta=27.9007;
    paramTab{2}.gamma=34.9018;
    paramTab{2}.nu=1.2615;
    paramTab{2}.eps=0.9444;
    paramTab{2}.sigma=35.6416;
    [~,~,paramTab{2}.alpha] = Rp_SEIIS_v4(paramTab{2}.beta,paramTab{2}.nu,paramTab{2}.eps,paramTab{2}.sigma,paramTab{2}.gamma,mu,b,0);
    
    paramTab{3}.betaI = 13.194799683714759;
    paramTab{3}.betaC = 1.406761238295781;
    paramTab{3}.sigma = 6.890538634865351;
    paramTab{3}.theta0 = 0.204544709146265; %%%%% verifier que c'est bien ça et pas gamma
    paramTab{3}.gamma0 = 0;
    paramTab{3}.zeta = 0.631739291427463;
    paramTab{3}.eta = 4;
    paramTab{3}.p = 0.2;
    [paramTab{3}.R_prep_0,~,paramTab{3}.alpha] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
        paramTab{3}.theta0,paramTab{3}.gamma0,paramTab{3}.sigma,paramTab{3}.zeta,...
        paramTab{3}.eta,paramTab{3}.p,mu,b,0);
    
    paramTab{4}.beta  = 0.9130;
    paramTab{4}.sigma = 9.4303;
    paramTab{4}.gamma10=0.;
    paramTab{4}.gamma30=0.0388;
    paramTab{4}.tau = 18.5060;
    paramTab{4}.theta=1.7839;
    paramTab{4}.nu = 0;
end

paramRho.eta_c_prep = 1;
paramRho.eta_s_prep = 1;
paramRho.eta_h_prep = paramTab{3}.eta;
paramRho.eta_g_prep = 4;
paramRho.eta_s_art  = 2;
paramRho.eta_c_art  = 0.1;
paramRho.eta_g_art  = 2;
paramRho.VTunderART  = 1;

paramRho.rho_h      = 0.0;
paramRho.rho_s      = 0.;
paramRho.rho_c      = 0.;
paramRho.rho_g      = 0.;
paramRho.rho_hs     = 0.0;
paramRho.rho_hc     = 0.0;
paramRho.rho_hg     = 0.0;
paramRho.rho_sc     = 0.0;
paramRho.rho_sg     = 0.0;
paramRho.rho_cg     = 0.0;
paramRho.rho_hsc    = 0.;
paramRho.rho_hsg    = 0.0;
paramRho.rho_hcg    = 0.0;
paramRho.rho_scg    = 0.0;
paramRho.rho_hscg   = 0.0;

N = b/mu;

% Input parameters (all in one)
allParameters4d = {paramTab{3}.betaI,paramTab{3}.betaC,paramTab{3}.sigma,paramTab{3}.theta0,paramTab{3}.zeta,paramTab{3}.eta,paramTab{3}.p,...
    paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.gamma30,paramTab{4}.tau,paramTab{4}.theta,...
    paramTab{1}.beta,paramTab{1}.gamma,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,...
    paramTab{2}.beta,paramTab{2}.gamma,paramTab{2}.nu,paramTab{2}.eps,paramTab{2}.sigma,...
    [],mu,b,...
    paramRho.rho_h, paramRho.rho_s, paramRho.rho_c, paramRho.rho_g,...
    paramRho.rho_hs, paramRho.rho_hc,paramRho.rho_hg,...
    paramRho.rho_sc, paramRho.rho_sg,...
    paramRho.rho_cg,...
    paramRho.rho_hsc, paramRho.rho_hsg, paramRho.rho_hcg, paramRho.rho_scg,...
    paramRho.rho_hscg,...
    paramRho.eta_s_prep, paramRho.eta_c_prep, paramRho.eta_g_prep,...
    paramRho.eta_s_art, paramRho.eta_c_art, paramRho.eta_g_art, paramRho.VTunderART};
%--------------------------------------------------------------------------%

%% Recreating tabComp
createTabComp;
tabComp_4dis = tabComp;

dis = ["HIV","syph"];
boxesSICTP  = ["S","I","C","P","Ip","Cp","T"]; %should be the same order than in the ODE_SICTPrEP.m
boxesSEIIIS = ["S","E","I1","I2","I3"];        %should be the same order than in the ODESEIIIS.m
boxesSEIIS  = ["S","E","IA","IS"];             %should be the same order than in the ODESEIIS.m
%Create the table of compartments
tabComp_hs = createTableComp(0,1,1,boxesSEIIS,boxesSICTP,boxesSEIIIS,dis);
tabComp_hs.no = [1:35]';

dis = ["HIV","Ct"];
boxesSICTP  = ["S","I","C","P","Ip","Cp","T"]; %should be the same order than in the ODE_SICTPrEP.m
boxesSEIIIS = ["S","E","I1","I2","I3"];        %should be the same order than in the ODESEIIIS.m
boxesSEIIS  = ["S","E","IA","IS"];             %should be the same order than in the ODESEIIS.m
%Create the table of compartments
tabComp_hc = createTableComp(1,1,0,boxesSEIIS,boxesSICTP,boxesSEIIIS,dis);
tabComp_hc.no = [1:28]';
%--------------------------------------------------------------------------%

% syms betaIh betaCh sigmah thetah zetah eta_h_prep ph,...
% syms betas sigmas gamma3s taus thetas,...
% syms betaCt gammaCt nuCt epsCt sigmaCt,...
% syms betaNg gammaNg nuNg epsNg sigmaNg,...
% syms tabComb mu b,...
% syms rho_h rho_s rho_c rho_g,...
% syms rho_hs rho_hc rho_hg,...
% syms rho_sc rho_sg,...
% syms rho_cg,...
% syms rho_hsc rho_hsg  rho_hcg  rho_scg,...
% syms rho_hscg,...
% syms eta_s_prep eta_c_prep eta_g_prep,...
% syms eta_s_art eta_c_art eta_g_art VTunderART

if (1)
    %% ODE system with the 4-disease model
    addpath 'C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ODEsystems\'
    Y0 = ones(560,1); tspan=[1,500];
    [res4dis] = ode45(@(t,Y) ODE_SICTPSEIIISSEIIS2_v7(t,Y,allParameters4d{:}),tspan, Y0);
    
    %Total of the population
    disp(['Pop. total_num - Pop tot theoric.:'])
    disp(sum(res4dis.y(:,end)) - b/mu)
    
    %% Comparaison avec les modèles à 1 infection
    [ES_sictp,ES_seiiis,ES_ct,ES_ng] = singleDisES(paramTab,paramRho,mu,b);
    % SICTP
    % Obtenu avec le 4-d-mod
    disp(['HIV pop., with the model of 4-dis.:'])
    disp([sum(res4dis.y(tabComp(tabComp.HIV=="S",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.HIV=="I",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.HIV=="C",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.HIV=="P",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.HIV=="Ip",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.HIV=="Cp",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.HIV=="T",:).no,end))]'/(b/mu))
    disp(['HIV pop., with the model SICTP:'])
    disp(ES_sictp/(b/mu));
    
    
    % SEIIIS
    disp(['Syphilis pop., with the model of 4-dis.:'])
    disp([sum(res4dis.y(tabComp(tabComp.syph=="S",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.syph=="E",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.syph=="I1",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.syph=="I2",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.syph=="I3",:).no,end))]'/(b/mu))
    
    disp(['Syphilis pop., with the model SEIIIS.:'])
    disp(ES_seiiis/(b/mu));
    
    % SEIIS (Ct)
    disp(['Ct pop., with the model of 4-dis.:'])
    disp([sum(res4dis.y(tabComp(tabComp.Ct=="S",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.Ct=="E",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.Ct=="IA",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.Ct=="IS",:).no,end))]'/(b/mu))
    
    disp(['Ct pop., with the model SEIIS:'])
    disp(ES_ct/(b/mu));
    
    % SEIIS (Ng)
    disp(['Ng pop., with the model of 4-dis.:'])
    disp([sum(res4dis.y(tabComp(tabComp.Ng=="S",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.Ng=="E",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.Ng=="IA",:).no,end));...
        sum(res4dis.y(tabComp(tabComp.Ng=="IS",:).no,end))]'/(b/mu))
    
    disp(['Ng pop., with the model SEIIS:'])
    disp(ES_ng/(b/mu));
    
    ES_hscg = res4dis.y(:,end);
    
end

%% Comparaison avec les modèles à 2 infections

addpath 'C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ODEsystems\'
Y0 = ones(35,1); tspan=[1,600];

allParameters_hs = {paramTab{3}.betaI,paramTab{3}.betaC,paramTab{3}.theta0,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.p,...
    paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.gamma30,paramTab{4}.tau,paramTab{4}.theta,...
    paramRho.rho_h, paramRho.rho_s,paramRho.rho_hs,...
    paramRho.eta_h_prep,paramRho.eta_s_prep,paramRho.eta_s_art,...
    paramRho.VTunderART,mu,b};

[ttot,res_hs] = ode45(@(t,Y) ODE_SICTPSEIIIS_v7(t,Y,allParameters_hs{:}),tspan,Y0);
%plot(ttot,res_hs)
ES_hs = res_hs(end,:);
disp('HIV x syphilis')

ES_hs_hscg=[];
i=1;
for stateS = ["S","E","I1","I2","I3"]
    for stateHIV=["S","I","C","P","Ip","Cp","T"]
        nos_stateHIV_S = tabComp(tabComp.HIV==stateHIV & tabComp.syph==stateS,:).no;
        %eqns = M_hscg*X;
        ES_hs_hscg(i) = sum(ES_hscg(nos_stateHIV_S));
        i=i+1;
    end
end

[tabComp_hs, array2table([ES_hs_hscg;ES_hs;ES_hs_hscg-ES_hs]'/N)]
max(abs(ES_hs_hscg-ES_hs)/N)

% rho_hsg=0;
% rho_hsc=0;
% rho_hscg=0;
% rho_sg = 0;
% rho_sc=0;
% rho_scg=0;
% rho_hg=0;
% rho_hc=0;

%(Y2 + Y57+ Y92 + Y127 + Y162 + Y197 + Y232 + Y267 + Y302 + Y337 + Y372 + Y407 + Y442 + Y477 + Y512 + Y547)*thetas - Lambdah*(Y64 + Y99 + Y134 + Y169 + Y204 + Y239 + Y274 + Y309 + Y344 + Y379 + Y414 + Y449 + Y484 + Y519 + Y554 + Y29) - (Y29 + Y64 + Y99 + Y134 + Y169 + Y204 + Y239 + Y274 + Y309 + Y344 + Y379 + Y414 + Y449 + Y484 + Y519 + Y554)*gamma3s - (Y29 + Y64 + Y99 + Y134 + Y169 + Y204 + Y239 + Y274 + Y309 + Y344 + Y379 + Y414 + Y449 + Y484 + Y519 + Y554)*mu - (Y29 + Y64 + Y99 + Y134 + Y169 + Y239)*rho_s - Y29*rho_sc - Y29*rho_sg - Y29*rho_hs - Y64*rho_sc - Y64*rho_sg - Y64*rho_hs - Y99*rho_sc - Y99*rho_sg - Y99*rho_hs - Y204*rho_s - Y29*rho_scg - Y29*rho_hsc - Y29*rho_hsg - Y134*rho_sg - Y134*rho_hs - Y64*rho_scg - Y64*rho_hsc - Y169*rho_sc - Y64*rho_hsg - Y169*rho_sg - Y169*rho_hs - Y274*rho_s - Y99*rho_scg - Y99*rho_hsc - Y204*rho_sc - Y99*rho_hsg - Y204*rho_sg - Y204*rho_hs - Y309*rho_s - Y29*rho_hscg - Y239*rho_sc - Y134*rho_hsg - Y239*rho_sg - Y239*rho_hs - Y344*rho_s - Y64*rho_hscg - Y169*rho_scg - Y169*rho_hsc - Y169*rho_hsg - Y274*rho_sg - Y274*rho_hs - Y379*rho_s - Y99*rho_hscg - Y204*rho_scg - Y204*rho_hsc - Y309*rho_sc - Y204*rho_hsg - Y309*rho_sg - Y309*rho_hs - Y414*rho_s - Y239*rho_scg - Y239*rho_hsc - Y344*rho_sc - Y239*rho_hsg - Y344*rho_sg - Y344*rho_hs - Y449*rho_s - Y169*rho_hscg - Y379*rho_sc - Y274*rho_hsg - Y379*rho_sg - Y379*rho_hs - Y484*rho_s - Y204*rho_hscg - Y309*rho_scg - Y309*rho_hsc - Y309*rho_hsg - Y414*rho_sg - Y414*rho_hs - Y519*rho_s - Y239*rho_hscg - Y344*rho_scg - Y344*rho_hsc - Y449*rho_sc - Y344*rho_hsg - Y449*rho_hs - Y554*rho_s - Y379*rho_scg - Y379*rho_hsc - Y484*rho_sc - Y379*rho_hsg - Y484*rho_hs - Y309*rho_hscg - Y519*rho_sc - Y414*rho_hsg - Y519*rho_hs - Y344*rho_hscg - Y449*rho_hsc - Y554*rho_hs - Y379*rho_hscg - Y484*rho_hsc - Y519*rho_hsc
% thetas*(Y2 + Y57 + Y92 + Y127 + Y162 + Y197 + Y232 + Y267 + Y302 + Y337 + Y372 + Y407 + Y442 + Y477 + Y512 + Y547) ...
% - (Y64 + Y99 + Y134 + Y169 + Y204 + Y239 + Y274 + Y309 + Y344 + Y379 + Y414 + Y29 + Y449+ Y484 + Y519 + Y554)*rho_hs ...
% - (Y204 + Y274 + Y309 + Y344 + Y379 + Y414 + Y449 + Y484 + Y519 + Y554 + Y29 + Y64 + Y99 + Y134 + Y169 + Y239) ...
% - Lambdah*(Y29 + Y64 + Y99 + Y134 + Y169 + Y204 + Y239 + Y274 + Y309 + Y344 + Y379 + Y414 + Y449 + Y484 + Y519 + Y554) ...
% - gamma3s*(Y29 + Y64 + Y99 + Y134 + Y169 + Y204 + Y239 + Y274 + Y309 + Y344 + Y379 + Y414 + Y449 + Y484 + Y519 + Y554) ...
% - mu*(Y29 + Y64 + Y99 + Y134 + Y169 + Y204 + Y239 + Y274 + Y309 + Y344 + Y379 + Y414 + Y449 + Y484 + Y519 + Y554);
%
% Y(12)*eta_s_prep + Y(19)*eta_s_prep + Y(26)*eta_s_prep - Y(5)*(Lambdas + eta_h_prep + mu + sigmah) + Y(33)*(eta_s_prep + gamma3s) - Lambdah*Y(4)*(zetah - 1);
%Y22*thetas - Lambdah*Y64 - Lambdah*Y99 - Lambdah*Y134 - Lambdah*Y169 - Lambdah*Y204 - Lambdah*Y239 - Lambdah*Y274 - Lambdah*Y309 - Lambdah*Y344 - Lambdah*Y379 - Lambdah*Y414 - Lambdah*Y449 - Lambdah*Y484 - Lambdah*Y519 - Lambdah*Y554 - Y29*gamma3s - Y64*gamma3s - Y99*gamma3s - Y134*gamma3s - Y169*gamma3s - Y204*gamma3s - Y239*gamma3s - Y274*gamma3s - Y309*gamma3s - Y344*gamma3s - Y379*gamma3s - Y414*gamma3s - Y449*gamma3s - Y484*gamma3s - Y519*gamma3s - Y554*gamma3s - Y29*mu - Y64*mu - Y99*mu - Y134*mu - Y169*mu - Y204*mu - Y239*mu - Y274*mu - Y309*mu - Y344*mu - Y379*mu - Y414*mu - Y449*mu - Y484*mu - Y519*mu - Y554*mu - Y29*rho_s - Y64*rho_s - Y99*rho_s - Y29*rho_sc - Y29*rho_sg - Y29*rho_hs - Y134*rho_s - Y64*rho_sc - Y64*rho_sg - Y64*rho_hs - Y169*rho_s - Y99*rho_sc - Y99*rho_sg - Y99*rho_hs - Y204*rho_s - Y29*rho_scg - Y29*rho_hsc - Y29*rho_hsg - Y134*rho_sg - Y134*rho_hs - Y239*rho_s - Y64*rho_scg - Y64*rho_hsc - Y169*rho_sc - Y64*rho_hsg - Y169*rho_sg - Y169*rho_hs - Y274*rho_s - Y99*rho_scg - Y99*rho_hsc - Y204*rho_sc - Y99*rho_hsg - Y204*rho_sg - Y204*rho_hs - Y309*rho_s - Y29*rho_hscg - Y239*rho_sc - Y134*rho_hsg - Y239*rho_sg - Y239*rho_hs - Y344*rho_s - Y64*rho_hscg - Y169*rho_scg - Y169*rho_hsc - Y169*rho_hsg - Y274*rho_sg - Y274*rho_hs - Y379*rho_s - Y99*rho_hscg - Y204*rho_scg - Y204*rho_hsc - Y309*rho_sc - Y204*rho_hsg - Y309*rho_sg - Y309*rho_hs - Y414*rho_s - Y239*rho_scg - Y239*rho_hsc - Y344*rho_sc - Y239*rho_hsg - Y344*rho_sg - Y344*rho_hs - Y449*rho_s - Y169*rho_hscg - Y379*rho_sc - Y274*rho_hsg - Y379*rho_sg - Y379*rho_hs - Y484*rho_s - Y204*rho_hscg - Y309*rho_scg - Y309*rho_hsc - Y309*rho_hsg - Y414*rho_sg - Y414*rho_hs - Y519*rho_s - Y239*rho_hscg - Y344*rho_scg - Y344*rho_hsc - Y449*rho_sc - Y344*rho_hsg - Y449*rho_hs - Y554*rho_s - Y379*rho_scg - Y379*rho_hsc - Y484*rho_sc - Y379*rho_hsg - Y484*rho_hs - Y309*rho_hscg - Y519*rho_sc - Y414*rho_hsg - Y519*rho_hs - Y344*rho_hscg - Y449*rho_hsc - Y554*rho_hs - Y379*rho_hscg - Y484*rho_hsc - Y519*rho_hsc - Lambdah*Y29 + Y57*thetas + Y92*thetas + Y127*thetas + Y162*thetas + Y197*thetas + Y232*thetas + Y267*thetas + Y302*thetas + Y337*thetas + Y372*thetas + Y407*thetas + Y442*thetas + Y477*thetas + Y512*thetas + Y547*thetas


Y0 = ones(28,1); tspan=[1,300];
allParameters_hc = {paramTab{3}.betaI,paramTab{3}.betaC,paramTab{3}.theta0,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.p,...
    paramTab{1}.beta,paramTab{1}.gamma,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,...
    paramRho.rho_h, paramRho.rho_c,paramRho.rho_hc,...
    paramRho.eta_h_prep,paramRho.eta_c_prep,paramRho.eta_c_art,...
    paramRho.VTunderART,mu,b};

[ttot,res_hc] = ode45(@(t,Y) ODE_SICTPSEIIS_v7(t,Y,allParameters_hc{:}),tspan,Y0);
ES_hc = res_hc(end,:);
disp('HIV x Ct')

i=1; ES_hc_hscg=[];
for stateC = ["S","E","IA","IS"]
    for stateHIV=["S","I","C","P","Ip","Cp","T"]
        nos_stateHIV_C = tabComp(tabComp.HIV==stateHIV & tabComp.Ct==stateC,:).no;
        %eqns = M_hscg*X;
        ES_hc_hscg(i) = sum(ES_hscg(nos_stateHIV_C));
        i=i+1;
    end
end

[tabComp_hc, array2table([ES_hc_hscg;ES_hc;ES_hc_hscg-ES_hc]'/N)]
max(abs(ES_hc_hscg-ES_hc)/N)


if (0)
    syms Y [560 1]
    syms X [35 1]
    syms betaIh betaCh thetah sigmah zetah ph
    syms betas sigmas gamma3s taus thetas
    syms rho_h rho_s rho_hs
    syms eta_h_prep eta_s_prep eta_s_art
    syms VTunderART mu b
    syms Lambdah
    syms Lambdas
    syms rho_hc
    syms rho_sc
    syms rho_sg
    syms rho_scg
    syms rho_hsc
    syms rho_hsg
    syms rho_hscg
    t = 0;
    dY_hs = ODE_SICTPSEIIIS_v7(t,Y,betaIh,betaCh,thetah,sigmah,zetah,ph,...
        betas,sigmas,gamma3s,taus,thetas,...
        rho_h,rho_s,rho_hs,eta_h_prep,eta_s_prep,eta_s_art,...
        VTunderART,mu,b);
    
    
    %S for HIV and S for syphilis
    
    dY_hs(1)
    SS_hs_hs = X29*(gamma3s + rho_s + rho_hs) - X1*(Lambdah + Lambdas + mu) - b*(ph - 1) + X8*(rho_s + rho_hs) + X15*(rho_s + rho_hs) + X22*(rho_s + rho_hs)
    
    %SS_4 = expand(simplify(sum(dY(nos_stateHIV_S))));
    SS_hs_hscg = b - Lambdah*Y1 - Lambdas*Y1 - Lambdah*Y36 - Lambdas*Y36 - Lambdah*Y71 - Lambdas*Y71 - Lambdah*Y106 - Lambdas*Y106 - Lambdah*Y141 - Lambdas*Y141 - Lambdah*Y176 - Lambdas*Y176 - Lambdah*Y211 - Lambdas*Y211 - Lambdah*Y246 - Lambdas*Y246 - Lambdah*Y281 - Lambdas*Y281 - Lambdah*Y316 - Lambdas*Y316 - Lambdah*Y351 - Lambdas*Y351 - Lambdah*Y386 - Lambdas*Y386 - Lambdah*Y421 - Lambdas*Y421 - Lambdah*Y456 - Lambdas*Y456 - Lambdah*Y491 - Lambdas*Y491 - Lambdah*Y526 - Lambdas*Y526 + Y29*gamma3s + Y64*gamma3s + Y99*gamma3s + Y134*gamma3s + Y169*gamma3s + Y204*gamma3s + Y239*gamma3s + Y274*gamma3s + Y309*gamma3s + Y344*gamma3s + Y379*gamma3s + Y414*gamma3s + Y449*gamma3s + Y484*gamma3s + Y519*gamma3s + Y554*gamma3s - Y1*mu - Y36*mu - Y71*mu - Y106*mu - Y141*mu - Y176*mu - Y211*mu - Y246*mu - Y281*mu - Y316*mu - Y351*mu - Y386*mu - Y421*mu - Y456*mu - Y491*mu - Y526*mu + Y8*rho_s + Y15*rho_s + Y22*rho_s + Y29*rho_s + Y43*rho_s + Y50*rho_s + Y57*rho_s + Y64*rho_s + Y78*rho_s + Y85*rho_s + Y92*rho_s + Y99*rho_s + Y8*rho_sc + Y8*rho_sg + Y8*rho_hs + Y113*rho_s + Y15*rho_sc + Y15*rho_sg + Y15*rho_hs + Y120*rho_s + Y22*rho_sc + Y22*rho_sg + Y22*rho_hs + Y127*rho_s + Y29*rho_sc + Y29*rho_sg + Y29*rho_hs + Y134*rho_s + Y43*rho_sc + Y43*rho_sg + Y43*rho_hs + Y148*rho_s + Y50*rho_sc + Y50*rho_sg + Y50*rho_hs + Y155*rho_s + Y57*rho_sc + Y57*rho_sg + Y57*rho_hs + Y162*rho_s + Y64*rho_sc + Y64*rho_sg + Y64*rho_hs + Y169*rho_s + Y78*rho_sc + Y78*rho_sg + Y78*rho_hs + Y183*rho_s + Y85*rho_sc + Y85*rho_sg + Y85*rho_hs + Y190*rho_s + Y92*rho_sc + Y92*rho_sg + Y92*rho_hs + Y197*rho_s + Y99*rho_sc + Y99*rho_sg + Y99*rho_hs + Y204*rho_s + Y8*rho_scg + Y8*rho_hsc + Y8*rho_hsg + Y113*rho_sg + Y15*rho_scg + Y113*rho_hs + Y15*rho_hsc + Y218*rho_s + Y15*rho_hsg + Y120*rho_sg + Y22*rho_scg + Y120*rho_hs + Y22*rho_hsc + Y225*rho_s + Y22*rho_hsg + Y127*rho_sg + Y29*rho_scg + Y127*rho_hs + Y29*rho_hsc + Y232*rho_s + Y29*rho_hsg + Y134*rho_sg + Y134*rho_hs + Y239*rho_s + Y43*rho_scg + Y43*rho_hsc + Y148*rho_sc + Y43*rho_hsg + Y148*rho_sg + Y50*rho_scg + Y148*rho_hs + Y50*rho_hsc + Y253*rho_s + Y155*rho_sc + Y50*rho_hsg + Y155*rho_sg + Y57*rho_scg + Y155*rho_hs + Y57*rho_hsc + Y260*rho_s + Y162*rho_sc + Y57*rho_hsg + Y162*rho_sg + Y64*rho_scg + Y162*rho_hs + Y64*rho_hsc + Y267*rho_s + Y169*rho_sc + Y64*rho_hsg + Y169*rho_sg + Y169*rho_hs + Y274*rho_s + Y78*rho_scg + Y78*rho_hsc + Y183*rho_sc + Y78*rho_hsg + Y183*rho_sg + Y85*rho_scg + Y183*rho_hs + Y85*rho_hsc + Y288*rho_s + Y190*rho_sc + Y85*rho_hsg + Y190*rho_sg + Y92*rho_scg + Y190*rho_hs + Y92*rho_hsc + Y295*rho_s + Y197*rho_sc + Y92*rho_hsg + Y197*rho_sg + Y99*rho_scg + Y197*rho_hs + Y99*rho_hsc + Y302*rho_s + Y204*rho_sc + Y99*rho_hsg + Y204*rho_sg + Y204*rho_hs + Y309*rho_s + Y8*rho_hscg + Y218*rho_sc + Y113*rho_hsg + Y15*rho_hscg + Y218*rho_sg + Y218*rho_hs + Y323*rho_s + Y225*rho_sc + Y120*rho_hsg + Y22*rho_hscg + Y225*rho_sg + Y225*rho_hs + Y330*rho_s + Y232*rho_sc + Y127*rho_hsg + Y29*rho_hscg + Y232*rho_sg + Y232*rho_hs + Y337*rho_s + Y239*rho_sc + Y134*rho_hsg + Y239*rho_sg + Y239*rho_hs + Y344*rho_s + Y43*rho_hscg + Y148*rho_scg + Y148*rho_hsc + Y148*rho_hsg + Y50*rho_hscg + Y253*rho_sg + Y155*rho_scg + Y253*rho_hs + Y155*rho_hsc + Y358*rho_s + Y155*rho_hsg + Y57*rho_hscg + Y260*rho_sg + Y162*rho_scg + Y260*rho_hs + Y162*rho_hsc + Y365*rho_s + Y162*rho_hsg + Y64*rho_hscg + Y267*rho_sg + Y169*rho_scg + Y267*rho_hs + Y169*rho_hsc + Y372*rho_s + Y169*rho_hsg + Y274*rho_sg + Y274*rho_hs + Y379*rho_s + Y78*rho_hscg + Y183*rho_scg + Y183*rho_hsc + Y288*rho_sc + Y183*rho_hsg + Y85*rho_hscg + Y288*rho_sg + Y190*rho_scg + Y288*rho_hs + Y190*rho_hsc + Y393*rho_s + Y295*rho_sc + Y190*rho_hsg + Y92*rho_hscg + Y295*rho_sg + Y197*rho_scg + Y295*rho_hs + Y197*rho_hsc + Y400*rho_s + Y302*rho_sc + Y197*rho_hsg + Y99*rho_hscg + Y302*rho_sg + Y204*rho_scg + Y302*rho_hs + Y204*rho_hsc + Y407*rho_s + Y309*rho_sc + Y204*rho_hsg + Y309*rho_sg + Y309*rho_hs + Y414*rho_s + Y218*rho_scg + Y218*rho_hsc + Y323*rho_sc + Y218*rho_hsg + Y323*rho_sg + Y225*rho_scg + Y323*rho_hs + Y225*rho_hsc + Y428*rho_s + Y330*rho_sc + Y225*rho_hsg + Y330*rho_sg + Y232*rho_scg + Y330*rho_hs + Y232*rho_hsc + Y435*rho_s + Y337*rho_sc + Y232*rho_hsg + Y337*rho_sg + Y239*rho_scg + Y337*rho_hs + Y239*rho_hsc + Y442*rho_s + Y344*rho_sc + Y239*rho_hsg + Y344*rho_sg + Y344*rho_hs + Y449*rho_s + Y148*rho_hscg + Y358*rho_sc + Y253*rho_hsg + Y155*rho_hscg + Y358*rho_sg + Y358*rho_hs + Y463*rho_s + Y365*rho_sc + Y260*rho_hsg + Y162*rho_hscg + Y365*rho_sg + Y365*rho_hs + Y470*rho_s + Y372*rho_sc + Y267*rho_hsg + Y169*rho_hscg + Y372*rho_sg + Y372*rho_hs + Y477*rho_s + Y379*rho_sc + Y274*rho_hsg + Y379*rho_sg + Y379*rho_hs + Y484*rho_s + Y183*rho_hscg + Y288*rho_scg + Y288*rho_hsc + Y288*rho_hsg + Y190*rho_hscg + Y393*rho_sg + Y295*rho_scg + Y393*rho_hs + Y295*rho_hsc + Y498*rho_s + Y295*rho_hsg + Y197*rho_hscg + Y400*rho_sg + Y302*rho_scg + Y400*rho_hs + Y302*rho_hsc + Y505*rho_s + Y302*rho_hsg + Y204*rho_hscg + Y407*rho_sg + Y309*rho_scg + Y407*rho_hs + Y309*rho_hsc + Y512*rho_s + Y309*rho_hsg + Y414*rho_sg + Y414*rho_hs + Y519*rho_s + Y218*rho_hscg + Y323*rho_scg + Y323*rho_hsc + Y428*rho_sc + Y323*rho_hsg + Y225*rho_hscg + Y330*rho_scg + Y428*rho_hs + Y330*rho_hsc + Y533*rho_s + Y435*rho_sc + Y330*rho_hsg + Y232*rho_hscg + Y337*rho_scg + Y435*rho_hs + Y337*rho_hsc + Y540*rho_s + Y442*rho_sc + Y337*rho_hsg + Y239*rho_hscg + Y344*rho_scg + Y442*rho_hs + Y344*rho_hsc + Y547*rho_s + Y449*rho_sc + Y344*rho_hsg + Y449*rho_hs + Y554*rho_s + Y358*rho_scg + Y358*rho_hsc + Y463*rho_sc + Y358*rho_hsg + Y365*rho_scg + Y463*rho_hs + Y365*rho_hsc + Y470*rho_sc + Y365*rho_hsg + Y372*rho_scg + Y470*rho_hs + Y372*rho_hsc + Y477*rho_sc + Y372*rho_hsg + Y379*rho_scg + Y477*rho_hs + Y379*rho_hsc + Y484*rho_sc + Y379*rho_hsg + Y484*rho_hs + Y288*rho_hscg + Y498*rho_sc + Y393*rho_hsg + Y295*rho_hscg + Y498*rho_hs + Y505*rho_sc + Y400*rho_hsg + Y302*rho_hscg + Y505*rho_hs + Y512*rho_sc + Y407*rho_hsg + Y309*rho_hscg + Y512*rho_hs + Y519*rho_sc + Y414*rho_hsg + Y519*rho_hs + Y323*rho_hscg + Y428*rho_hsc + Y330*rho_hscg + Y533*rho_hs + Y435*rho_hsc + Y337*rho_hscg + Y540*rho_hs + Y442*rho_hsc + Y344*rho_hscg + Y547*rho_hs + Y449*rho_hsc + Y554*rho_hs + Y358*rho_hscg + Y463*rho_hsc + Y365*rho_hscg + Y470*rho_hsc + Y372*rho_hscg + Y477*rho_hsc + Y379*rho_hscg + Y484*rho_hsc + Y498*rho_hsc + Y505*rho_hsc + Y512*rho_hsc + Y519*rho_hsc - b*ph
    
    Lambdas=0
    Lambdah=0
    mu=0
    gamma3s=0
    rho_s=0
    rho_h=0
    rho_hs=0
    
    simplify(SS_hs_hs - SS_hs_hscg)
    Lambdah*Y1 - Lambdas*X1 - Lambdah*X1 + Lambdas*Y1 + Lambdah*Y36 + Lambdas*Y36 + Lambdah*Y71 + Lambdas*Y71 + Lambdah*Y106 + Lambdas*Y106 + Lambdah*Y141 + Lambdas*Y141 + Lambdah*Y176 + Lambdas*Y176 + Lambdah*Y211 + Lambdas*Y211 + Lambdah*Y246 + Lambdas*Y246 + Lambdah*Y281 + Lambdas*Y281 + Lambdah*Y316 + Lambdas*Y316 + Lambdah*Y351 + Lambdas*Y351 + Lambdah*Y386 + Lambdas*Y386 + Lambdah*Y421 + Lambdas*Y421 + Lambdah*Y456 + Lambdas*Y456 + Lambdah*Y491 + Lambdas*Y491 + Lambdah*Y526 + Lambdas*Y526 + X29*gamma3s - Y29*gamma3s - Y64*gamma3s - Y99*gamma3s - Y134*gamma3s - Y169*gamma3s - Y204*gamma3s - Y239*gamma3s - Y274*gamma3s - Y309*gamma3s - Y344*gamma3s - Y379*gamma3s - Y414*gamma3s - Y449*gamma3s - Y484*gamma3s - Y519*gamma3s - Y554*gamma3s - X1*mu + Y1*mu + Y36*mu + Y71*mu + Y106*mu + Y141*mu + Y176*mu + Y211*mu + Y246*mu + Y281*mu + Y316*mu + Y351*mu + Y386*mu + Y421*mu + Y456*mu + Y491*mu + Y526*mu + X8*rho_s + X15*rho_s + X22*rho_s + X29*rho_s + X8*rho_hs + X15*rho_hs + X22*rho_hs + X29*rho_hs - Y8*rho_s - Y15*rho_s - Y22*rho_s - Y29*rho_s - Y43*rho_s - Y50*rho_s - Y57*rho_s - Y64*rho_s - Y78*rho_s - Y85*rho_s - Y92*rho_s - Y99*rho_s - Y8*rho_sc - Y8*rho_sg - Y8*rho_hs - Y113*rho_s - Y15*rho_sc - Y15*rho_sg - Y15*rho_hs - Y120*rho_s - Y22*rho_sc - Y22*rho_sg - Y22*rho_hs - Y127*rho_s - Y29*rho_sc - Y29*rho_sg - Y29*rho_hs - Y134*rho_s - Y43*rho_sc - Y43*rho_sg - Y43*rho_hs - Y148*rho_s - Y50*rho_sc - Y50*rho_sg - Y50*rho_hs - Y155*rho_s - Y57*rho_sc - Y57*rho_sg - Y57*rho_hs - Y162*rho_s - Y64*rho_sc - Y64*rho_sg - Y64*rho_hs - Y169*rho_s - Y78*rho_sc - Y78*rho_sg - Y78*rho_hs - Y183*rho_s - Y85*rho_sc - Y85*rho_sg - Y85*rho_hs - Y190*rho_s - Y92*rho_sc - Y92*rho_sg - Y92*rho_hs - Y197*rho_s - Y99*rho_sc - Y99*rho_sg - Y99*rho_hs - Y204*rho_s - Y8*rho_scg - Y8*rho_hsc - Y8*rho_hsg - Y113*rho_sg - Y15*rho_scg - Y113*rho_hs - Y15*rho_hsc - Y218*rho_s - Y15*rho_hsg - Y120*rho_sg - Y22*rho_scg - Y120*rho_hs - Y22*rho_hsc - Y225*rho_s - Y22*rho_hsg - Y127*rho_sg - Y29*rho_scg - Y127*rho_hs - Y29*rho_hsc - Y232*rho_s - Y29*rho_hsg - Y134*rho_sg - Y134*rho_hs - Y239*rho_s - Y43*rho_scg - Y43*rho_hsc - Y148*rho_sc - Y43*rho_hsg - Y148*rho_sg - Y50*rho_scg - Y148*rho_hs - Y50*rho_hsc - Y253*rho_s - Y155*rho_sc - Y50*rho_hsg - Y155*rho_sg - Y57*rho_scg - Y155*rho_hs - Y57*rho_hsc - Y260*rho_s - Y162*rho_sc - Y57*rho_hsg - Y162*rho_sg - Y64*rho_scg - Y162*rho_hs - Y64*rho_hsc - Y267*rho_s - Y169*rho_sc - Y64*rho_hsg - Y169*rho_sg - Y169*rho_hs - Y274*rho_s - Y78*rho_scg - Y78*rho_hsc - Y183*rho_sc - Y78*rho_hsg - Y183*rho_sg - Y85*rho_scg - Y183*rho_hs - Y85*rho_hsc - Y288*rho_s - Y190*rho_sc - Y85*rho_hsg - Y190*rho_sg - Y92*rho_scg - Y190*rho_hs - Y92*rho_hsc - Y295*rho_s - Y197*rho_sc - Y92*rho_hsg - Y197*rho_sg - Y99*rho_scg - Y197*rho_hs - Y99*rho_hsc - Y302*rho_s - Y204*rho_sc - Y99*rho_hsg - Y204*rho_sg - Y204*rho_hs - Y309*rho_s - Y8*rho_hscg - Y218*rho_sc - Y113*rho_hsg - Y15*rho_hscg - Y218*rho_sg - Y218*rho_hs - Y323*rho_s - Y225*rho_sc - Y120*rho_hsg - Y22*rho_hscg - Y225*rho_sg - Y225*rho_hs - Y330*rho_s - Y232*rho_sc - Y127*rho_hsg - Y29*rho_hscg - Y232*rho_sg - Y232*rho_hs - Y337*rho_s - Y239*rho_sc - Y134*rho_hsg - Y239*rho_sg - Y239*rho_hs - Y344*rho_s - Y43*rho_hscg - Y148*rho_scg - Y148*rho_hsc - Y148*rho_hsg - Y50*rho_hscg - Y253*rho_sg - Y155*rho_scg - Y253*rho_hs - Y155*rho_hsc - Y358*rho_s - Y155*rho_hsg - Y57*rho_hscg - Y260*rho_sg - Y162*rho_scg - Y260*rho_hs - Y162*rho_hsc - Y365*rho_s - Y162*rho_hsg - Y64*rho_hscg - Y267*rho_sg - Y169*rho_scg - Y267*rho_hs - Y169*rho_hsc - Y372*rho_s - Y169*rho_hsg - Y274*rho_sg - Y274*rho_hs - Y379*rho_s - Y78*rho_hscg - Y183*rho_scg - Y183*rho_hsc - Y288*rho_sc - Y183*rho_hsg - Y85*rho_hscg - Y288*rho_sg - Y190*rho_scg - Y288*rho_hs - Y190*rho_hsc - Y393*rho_s - Y295*rho_sc - Y190*rho_hsg - Y92*rho_hscg - Y295*rho_sg - Y197*rho_scg - Y295*rho_hs - Y197*rho_hsc - Y400*rho_s - Y302*rho_sc - Y197*rho_hsg - Y99*rho_hscg - Y302*rho_sg - Y204*rho_scg - Y302*rho_hs - Y204*rho_hsc - Y407*rho_s - Y309*rho_sc - Y204*rho_hsg - Y309*rho_sg - Y309*rho_hs - Y414*rho_s - Y218*rho_scg - Y218*rho_hsc - Y323*rho_sc - Y218*rho_hsg - Y323*rho_sg - Y225*rho_scg - Y323*rho_hs - Y225*rho_hsc - Y428*rho_s - Y330*rho_sc - Y225*rho_hsg - Y330*rho_sg - Y232*rho_scg - Y330*rho_hs - Y232*rho_hsc - Y435*rho_s - Y337*rho_sc - Y232*rho_hsg - Y337*rho_sg - Y239*rho_scg - Y337*rho_hs - Y239*rho_hsc - Y442*rho_s - Y344*rho_sc - Y239*rho_hsg - Y344*rho_sg - Y344*rho_hs - Y449*rho_s - Y148*rho_hscg - Y358*rho_sc - Y253*rho_hsg - Y155*rho_hscg - Y358*rho_sg - Y358*rho_hs - Y463*rho_s - Y365*rho_sc - Y260*rho_hsg - Y162*rho_hscg - Y365*rho_sg - Y365*rho_hs - Y470*rho_s - Y372*rho_sc - Y267*rho_hsg - Y169*rho_hscg - Y372*rho_sg - Y372*rho_hs - Y477*rho_s - Y379*rho_sc - Y274*rho_hsg - Y379*rho_sg - Y379*rho_hs - Y484*rho_s - Y183*rho_hscg - Y288*rho_scg - Y288*rho_hsc - Y288*rho_hsg - Y190*rho_hscg - Y393*rho_sg - Y295*rho_scg - Y393*rho_hs - Y295*rho_hsc - Y498*rho_s - Y295*rho_hsg - Y197*rho_hscg - Y400*rho_sg - Y302*rho_scg - Y400*rho_hs - Y302*rho_hsc - Y505*rho_s - Y302*rho_hsg - Y204*rho_hscg - Y407*rho_sg - Y309*rho_scg - Y407*rho_hs - Y309*rho_hsc - Y512*rho_s - Y309*rho_hsg - Y414*rho_sg - Y414*rho_hs - Y519*rho_s - Y218*rho_hscg - Y323*rho_scg - Y323*rho_hsc - Y428*rho_sc - Y323*rho_hsg - Y225*rho_hscg - Y330*rho_scg - Y428*rho_hs - Y330*rho_hsc - Y533*rho_s - Y435*rho_sc - Y330*rho_hsg - Y232*rho_hscg - Y337*rho_scg - Y435*rho_hs - Y337*rho_hsc - Y540*rho_s - Y442*rho_sc - Y337*rho_hsg - Y239*rho_hscg - Y344*rho_scg - Y442*rho_hs - Y344*rho_hsc - Y547*rho_s - Y449*rho_sc - Y344*rho_hsg - Y449*rho_hs - Y554*rho_s - Y358*rho_scg - Y358*rho_hsc - Y463*rho_sc - Y358*rho_hsg - Y365*rho_scg - Y463*rho_hs - Y365*rho_hsc - Y470*rho_sc - Y365*rho_hsg - Y372*rho_scg - Y470*rho_hs - Y372*rho_hsc - Y477*rho_sc - Y372*rho_hsg - Y379*rho_scg - Y477*rho_hs - Y379*rho_hsc - Y484*rho_sc - Y379*rho_hsg - Y484*rho_hs - Y288*rho_hscg - Y498*rho_sc - Y393*rho_hsg - Y295*rho_hscg - Y498*rho_hs - Y505*rho_sc - Y400*rho_hsg - Y302*rho_hscg - Y505*rho_hs - Y512*rho_sc - Y407*rho_hsg - Y309*rho_hscg - Y512*rho_hs - Y519*rho_sc - Y414*rho_hsg - Y519*rho_hs - Y323*rho_hscg - Y428*rho_hsc - Y330*rho_hscg - Y533*rho_hs - Y435*rho_hsc - Y337*rho_hscg - Y540*rho_hs - Y442*rho_hsc - Y344*rho_hscg - Y547*rho_hs - Y449*rho_hsc - Y554*rho_hs - Y358*rho_hscg - Y463*rho_hsc - Y365*rho_hscg - Y470*rho_hsc - Y372*rho_hscg - Y477*rho_hsc - Y379*rho_hscg - Y484*rho_hsc - Y498*rho_hsc - Y505*rho_hsc - Y512*rho_hsc - Y519*rho_hsc
    
    rho_sc=0
    rho_sg=0
    rho_scg=0
    rho_hsc=0
    rho_hscg=0
    rho_hsg=0
    - Y8*rho_sc - Y8*rho_sg - Y15*rho_sc - Y15*rho_sg - Y22*rho_sc - Y22*rho_sg - Y29*rho_sc - Y29*rho_sg - Y43*rho_sc - Y43*rho_sg - Y50*rho_sc - Y50*rho_sg - Y57*rho_sc - Y57*rho_sg - Y64*rho_sc - Y64*rho_sg - Y78*rho_sc - Y78*rho_sg - Y85*rho_sc - Y85*rho_sg - Y92*rho_sc - Y92*rho_sg - Y99*rho_sc - Y99*rho_sg - Y8*rho_scg - Y8*rho_hsc - Y8*rho_hsg - Y113*rho_sg - Y15*rho_scg - Y15*rho_hsc - Y15*rho_hsg - Y120*rho_sg - Y22*rho_scg - Y22*rho_hsc - Y22*rho_hsg - Y127*rho_sg - Y29*rho_scg - Y29*rho_hsc - Y29*rho_hsg - Y134*rho_sg - Y43*rho_scg - Y43*rho_hsc - Y148*rho_sc - Y43*rho_hsg - Y148*rho_sg - Y50*rho_scg - Y50*rho_hsc - Y155*rho_sc - Y50*rho_hsg - Y155*rho_sg - Y57*rho_scg - Y57*rho_hsc - Y162*rho_sc - Y57*rho_hsg - Y162*rho_sg - Y64*rho_scg - Y64*rho_hsc - Y169*rho_sc - Y64*rho_hsg - Y169*rho_sg - Y78*rho_scg - Y78*rho_hsc - Y183*rho_sc - Y78*rho_hsg - Y183*rho_sg - Y85*rho_scg - Y85*rho_hsc - Y190*rho_sc - Y85*rho_hsg - Y190*rho_sg - Y92*rho_scg - Y92*rho_hsc - Y197*rho_sc - Y92*rho_hsg - Y197*rho_sg - Y99*rho_scg - Y99*rho_hsc - Y204*rho_sc - Y99*rho_hsg - Y204*rho_sg - Y8*rho_hscg - Y218*rho_sc - Y113*rho_hsg - Y15*rho_hscg - Y218*rho_sg - Y225*rho_sc - Y120*rho_hsg - Y22*rho_hscg - Y225*rho_sg - Y232*rho_sc - Y127*rho_hsg - Y29*rho_hscg - Y232*rho_sg - Y239*rho_sc - Y134*rho_hsg - Y239*rho_sg - Y43*rho_hscg - Y148*rho_scg - Y148*rho_hsc - Y148*rho_hsg - Y50*rho_hscg - Y253*rho_sg - Y155*rho_scg - Y155*rho_hsc - Y155*rho_hsg - Y57*rho_hscg - Y260*rho_sg - Y162*rho_scg - Y162*rho_hsc - Y162*rho_hsg - Y64*rho_hscg - Y267*rho_sg - Y169*rho_scg - Y169*rho_hsc - Y169*rho_hsg - Y274*rho_sg - Y78*rho_hscg - Y183*rho_scg - Y183*rho_hsc - Y288*rho_sc - Y183*rho_hsg - Y85*rho_hscg - Y288*rho_sg - Y190*rho_scg - Y190*rho_hsc - Y295*rho_sc - Y190*rho_hsg - Y92*rho_hscg - Y295*rho_sg - Y197*rho_scg - Y197*rho_hsc - Y302*rho_sc - Y197*rho_hsg - Y99*rho_hscg - Y302*rho_sg - Y204*rho_scg - Y204*rho_hsc - Y309*rho_sc - Y204*rho_hsg - Y309*rho_sg - Y218*rho_scg - Y218*rho_hsc - Y323*rho_sc - Y218*rho_hsg - Y323*rho_sg - Y225*rho_scg - Y225*rho_hsc - Y330*rho_sc - Y225*rho_hsg - Y330*rho_sg - Y232*rho_scg - Y232*rho_hsc - Y337*rho_sc - Y232*rho_hsg - Y337*rho_sg - Y239*rho_scg - Y239*rho_hsc - Y344*rho_sc - Y239*rho_hsg - Y344*rho_sg - Y148*rho_hscg - Y358*rho_sc - Y253*rho_hsg - Y155*rho_hscg - Y358*rho_sg - Y365*rho_sc - Y260*rho_hsg - Y162*rho_hscg - Y365*rho_sg - Y372*rho_sc - Y267*rho_hsg - Y169*rho_hscg - Y372*rho_sg - Y379*rho_sc - Y274*rho_hsg - Y379*rho_sg - Y183*rho_hscg - Y288*rho_scg - Y288*rho_hsc - Y288*rho_hsg - Y190*rho_hscg - Y393*rho_sg - Y295*rho_scg - Y295*rho_hsc - Y295*rho_hsg - Y197*rho_hscg - Y400*rho_sg - Y302*rho_scg - Y302*rho_hsc - Y302*rho_hsg - Y204*rho_hscg - Y407*rho_sg - Y309*rho_scg - Y309*rho_hsc - Y309*rho_hsg - Y414*rho_sg - Y218*rho_hscg - Y323*rho_scg - Y323*rho_hsc - Y428*rho_sc - Y323*rho_hsg - Y225*rho_hscg - Y330*rho_scg - Y330*rho_hsc - Y435*rho_sc - Y330*rho_hsg - Y232*rho_hscg - Y337*rho_scg - Y337*rho_hsc - Y442*rho_sc - Y337*rho_hsg - Y239*rho_hscg - Y344*rho_scg - Y344*rho_hsc - Y449*rho_sc - Y344*rho_hsg - Y358*rho_scg - Y358*rho_hsc - Y463*rho_sc - Y358*rho_hsg - Y365*rho_scg - Y365*rho_hsc - Y470*rho_sc - Y365*rho_hsg - Y372*rho_scg - Y372*rho_hsc - Y477*rho_sc - Y372*rho_hsg - Y379*rho_scg - Y379*rho_hsc - Y484*rho_sc - Y379*rho_hsg - Y288*rho_hscg - Y498*rho_sc - Y393*rho_hsg - Y295*rho_hscg - Y505*rho_sc - Y400*rho_hsg - Y302*rho_hscg - Y512*rho_sc - Y407*rho_hsg - Y309*rho_hscg - Y519*rho_sc - Y414*rho_hsg - Y323*rho_hscg - Y428*rho_hsc - Y330*rho_hscg - Y435*rho_hsc - Y337*rho_hscg - Y442*rho_hsc - Y344*rho_hscg - Y449*rho_hsc - Y358*rho_hscg - Y463*rho_hsc - Y365*rho_hscg - Y470*rho_hsc - Y372*rho_hscg - Y477*rho_hsc - Y379*rho_hscg - Y484*rho_hsc - Y498*rho_hsc - Y505*rho_hsc - Y512*rho_hsc - Y519*rho_hsc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tout ce qui est en dessous, n'a pas été vraiment testé et adapté à la
%%% version v5_bis de P1234, ODE...
%%% Ca vient de script_SEIISm_SICTP_SEIIIS_v5_bis.suite_tests.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(0)
    if(0)
        
        %% Comparaison entre fsolve, ode45 et knitro
        % Initial condition
        Y0 = ones(560,1)/560;
        solveur = "ode45"; %fsolve, ode45 knitro, knitro_ampl
        tspan=[0,50];
        
        N_th=b/mu;
        restart=true; iterNo=0;
        while restart && iterNo<15
            tic;
            iterNo=iterNo+1
            %     if iterNo>10
            %         solveur="ode45";
            %     end
            if solveur=="fsolve"
                options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-6,'MaxFunctionEvaluations',100000,...
                    'Algorithm','trust-region','SubproblemAlgorithm','cg');
                [ES,fval,exitflag,output] = fsolve(@(Y) ODE_SICTPSEIIISSEIIS2_v5_bis(0,Y,allParameters4d{:}),Y0,options);
                Y0 = rand(560,1)*b/mu; %pour le coup d'après
                
            elseif solveur=="knitro"
                my_f = @(Y) ODE_SICTPSEIIISSEIIS2_v5_bis(0,Y,allParameters4d{:});
                options = knitro_options('maxtime_real',100.0);
                %options = knitro_options('maxtime_real',10.0, 'ms_enable', 1, 'ms_maxsolves', 5, 'ms_maxtime_real', 500);
                [ES,fval,exitflag,output] = knitro_nlneqs(myf, Y0, {}, options);
                objfunc = @(Y) ([1]);
                confunc = @(Y) (disperse2([], myf(Y)));
                %objfunc = @(Y) ([norm(myf(Y))]);
                %confunc = @(Y) (disperse2([], []));
                lb=zeros(size(Y0));
                [ES,fval,exitflag,output,lambda,grad,hessian] = knitro_nlp(objfunc, Y0,...
                    [], [], [], [],...
                    lb, [], confunc, [], options);
                disp(['error=', num2str(max(abs(myf(ES))))])
                Y0 = rand(560,1)*b/mu; %pour le coup d'apres
                
            elseif solveur=="ode45"
                [res] = ode45(@(t,Y) ODE_SICTPSEIIISSEIIS2_v5_bis(t,Y,allParameters4d{:}),tspan, Y0);
                ES = res.y(:,end);
                Y0 = ES; %pour continuer la résolution
            end
            
            %only HIV:
            res_4_h=[sum(ES(1:7:554)), sum(ES(2:7:555)), sum(ES(3:7:556)), sum(ES(4:7:557)),...
                sum(ES(5:7:558)), sum(ES(6:7:559)), sum(ES(7:7:560))];
            disp(['4dis,HIV: ' ,num2str(res_4_h)])
            %from the sictp model
            disp(['SICTP: ' ,num2str(ES_sictp)])
            
            %only syph:
            res_4_s = [sum(ES(reshape(repmat((1:35:554),7,1)+[0:6]',[1,560/5]))),...
                sum(ES(reshape(repmat((8:35:561),7,1)+[0:6]',[1,560/5]))),...
                sum(ES(reshape(repmat((15:35:568),7,1)+[0:6]',[1,560/5]))),...
                sum(ES(reshape(repmat((22:35:569),7,1)+[0:6]',[1,560/5]))),...
                sum(ES(reshape(repmat((29:35:576),7,1)+[0:6]',[1,560/5])))];
            disp(['4dis,s: ' ,num2str(res_4_s)])
            
            %from the seiiis model
            disp(['SEIIIS: ' ,num2str(ES_seiiis)])
            
            % seiis model Ct
            %only Ct:
            res_4_Ct = [sum(ES(reshape(repmat((1:140:421),35,1)+[0:34]',[1,560/4]))),...
                sum(ES(reshape(repmat((36:140:456),35,1)+[0:34]',[1,560/4]))),...
                sum(ES(reshape(repmat((71:140:491),35,1)+[0:34]',[1,560/4]))),...
                sum(ES(reshape(repmat((106:140:527),35,1)+[0:34]',[1,560/4])))];
            disp(['4dis,Ct: ' ,num2str(res_4_Ct)])
            %from the Ct model
            disp(['SEIIS,Ct: ' ,num2str(ES_ct)])
            
            % seiis model Ng
            res_4_Ng = [sum(ES(1:140)),...
                sum(ES(141:280)),...
                sum(ES(281:420)),...
                sum(ES(421:560))];
            disp(['4dis,Ng: ' ,num2str(res_4_Ng)])
            %from the Ng model
            disp(['SEIIS,Ng: ' ,num2str(ES_ng)])
            
            myTol=1e-4;
            popTot = sum(ES(1:560));
            if (abs(popTot-N_th)>myTol || sum(ES<-myTol)~=0)
                tspan=[tspan(end),tspan(end)+50];
            else
                restart=false;
            end
            toc
        end
        
        
        %% Comparaison avec AMPL
        
        % Lancer le système d'ODE avec AMPL et comparer avec celui de Matlab
        % en terme de résultats et vitesse de calcul
        
        % Insérer le code : comparaison_avec_ampl_test
        
    end
    
    
    
    
    
    %% Comparison of the different methods: fsolve (matlab), ode45 (matlab), knitro (ampl)
    clear all;
    addpath 'C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ODEsystems\'
    f=1; nbSimu=1; outputSolve = zeros(nbSimu,4,3); maxIter=15;
    b=100; P_inf_knitro_tot=[];
    pHIV=0;
    for simu=1:nbSimu
        %------------------------%
        [paramTab,mu,~] = sampleParameters_v3_extent(true,true,true,true,b,pHIV);   %Ct,Ng,HIV,syph
        createParamRho;
        %------------------------%
        
        tic;
        %[ES_fsolve,~,P_inf_fsolve,lastIter,optimFound] = P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,f,'fsolve',maxIter);
        %outputSolve(simu,1,1) = toc; outputSolve(simu,2,1) = lastIter; outputSolve(simu,3,1) = ~optimFound;
        tic
        %[ES_ode45,~,P_inf_ode45,lastIter,optimFound] = P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,f,'ode45',maxIter);
        %outputSolve(simu,1,2) = toc; outputSolve(simu,2,2) = lastIter; outputSolve(simu,3,2) = ~optimFound;
        tic;
        [ES_knitro,~,P_inf_knitro,lastIter,optimFound,changeSolver] = P1234_SICTPSEIIISSEIIS2_v5_bis(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,...
            paramRho,f,'knitroampl',maxIter);
        outputSolve(simu,1,3) = toc; outputSolve(simu,2,3) = lastIter; outputSolve(simu,3,3) = ~optimFound; outputSolve(simu,4,3) = changeSolver;
        P_inf_knitro_tot = [P_inf_knitro_tot;P_inf_knitro];
        
        %     if (max(abs(ES_ode45(1:560)-ES_knitro))>0.01*b/mu)
        %         break;
        %     end
        
        %     disp(P_inf_ode45);
        disp(P_inf_knitro);
        
        [ES_sictp,ES_seiiis,ES_ct,ES_ng] = singleDisES(paramTab,mu,b);
        P_sictp = 1 - sum(ES_sictp([1,4,7])/(b/mu));
        P_seiiis = 1 - sum(ES_seiiis(1)/(b/mu));
        P_ct = 1 - sum(ES_ct([1,4])/(b/mu));
        P_ng = 1 - sum(ES_ng([1,4])/(b/mu));
        disp([P_sictp,P_seiiis,P_ct,P_ng])
    end
    
    createTabComp;
    tabComp.ES = round(ES_knitro/sum(ES_knitro),4)*100;
    tabComp(tabComp.ES>0,:)
    
    
    if(1)
        %% The utility maximization problem with AMPL/knitro
        %clear all;
        addpath('C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/matlab/');
        addpath('C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/lib/');
        addpath('C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab');
        addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\AMPL_models')
        setupOnce;
        
        ampl = AMPL;
        ampl.read('C:/Users/Moi/Documents/IPLESP/These/Codes/AMPL_models/main_prev.mod'); %pour representer U en fonction de rho
        
        paramHIV=paramTab{3};paramS=paramTab{4};paramC=paramTab{1};paramG=paramTab{2};
        assigningParameters;
        createTabComp
        
        %% 1. Show U for some specific value of c
        c=0;
        vecRho = 0:0.01:0.5;
        i=0;
        knitro_options = 'ms_enable=1 ms_maxsolves=5 feastol=1e-6 maxtime_real=20 ms_maxtime_real=20 ms_outsub=1 ms_numthreads=4 ncvx_qcqp_init=0';
        
        ES = zeros(length(vecRho),560);
        
        for rho=vecRho
            disp(' ')
            disp('############################')
            disp(['rho=',num2str(rho),' iter=', num2str(i-1),' on ', num2str(length(vecRho))])
            i=i+1;
            rho_to_change   = ampl.getParameter('rho_g'); rho_to_change.setValues([rho]);
            ampl.setOption('knitro_options', [knitro_options]);
            ampl.solve();
            status = ampl.getValue("solve_result_num");
            
            if status ~= 0
                disp("/n#### Non-optimal status, check multi-start procedure. ###/n");
                ampl.close();
                break;
            end
            
            ES_var = ampl.getVariable('Y');
            df = ES_var.getValues;
            a = df.val;
            ES(i,:) = cell2mat(a)';
            current_ES = ES(i,:);
            
            disp(' ')
            disp(['Pop. total_num - Pop tot theoric.:'])
            disp(abs(sum(current_ES) - b/mu))
            
            % Comparaison avec les modèles à 1 infection
            paramTab{1}.rhob = rho;
            paramTab{2}.rhob = 0;
            [ES_sictp,ES_seiiis,ES_ct,ES_ng] = singleDisES(paramTab,mu,b);
            
            % SICTP
            disp(['HIV pop., with the model of 4-dis.:'])
            disp([sum(current_ES(tabComp(tabComp.HIV=="S",:).no));...
                sum(current_ES(tabComp(tabComp.HIV=="I",:).no));...
                sum(current_ES(tabComp(tabComp.HIV=="C",:).no));...
                sum(current_ES(tabComp(tabComp.HIV=="P",:).no));...
                sum(current_ES(tabComp(tabComp.HIV=="Ip",:).no));...
                sum(current_ES(tabComp(tabComp.HIV=="Cp",:).no));...
                sum(current_ES(tabComp(tabComp.HIV=="T",:).no))]'/(b/mu))
            disp(['HIV pop., with the model SICTP:'])
            disp(ES_sictp/(b/mu));
            % SEIIIS
            disp(['Syphilis pop., with the model of 4-dis.:'])
            disp([sum(current_ES(tabComp(tabComp.syph=="S",:).no));...
                sum(current_ES(tabComp(tabComp.syph=="E",:).no));...
                sum(current_ES(tabComp(tabComp.syph=="I1",:).no));...
                sum(current_ES(tabComp(tabComp.syph=="I2",:).no));...
                sum(current_ES(tabComp(tabComp.syph=="I3",:).no))]'/(b/mu))
            disp(['Syphilis pop., with the model SEIIIS.:'])
            disp(ES_seiiis/(b/mu));
            % SEIIS (Ct)
            disp(['Ct pop., with the model of 4-dis.:'])
            disp([sum(current_ES(tabComp(tabComp.Ct=="S",:).no));...
                sum(current_ES(tabComp(tabComp.Ct=="E",:).no));...
                sum(current_ES(tabComp(tabComp.Ct=="IA",:).no));...
                sum(current_ES(tabComp(tabComp.Ct=="IS",:).no))]'/(b/mu))
            disp(['Ct pop., with the model SEIIS:'])
            disp(ES_ct/(b/mu));
            % SEIIS (Ng)
            disp(['Ng pop., with the model of 4-dis.:'])
            disp([sum(current_ES(tabComp(tabComp.Ng=="S",:).no));...
                sum(current_ES(tabComp(tabComp.Ng=="E",:).no));...
                sum(current_ES(tabComp(tabComp.Ng=="IA",:).no));...
                sum(current_ES(tabComp(tabComp.Ng=="IS",:).no))]'/(b/mu))
            disp(['Ng pop., with the model SEIIS:'])
            disp(ES_ng/(b/mu));
            
        end
        
        disp(['Pop. total_num - Pop tot theoric.:'])
        disp(max(abs(sum(ES,2) - b/mu)))
        ampl.close();
        
        
        %dY = ODE_SICTPSEIIISSEIIS2_v5(1,ES,allParameters4d{:});
        %max(abs(dY))
        
        %
        figure(1)
        prev_g=[];
        for i=1:length(vecRho)
            prev_g(i) = sum([sum(ES(i,tabComp(tabComp.Ng=="E",:).no),2);...
                sum(ES(i,tabComp(tabComp.Ng=="IA",:).no),2);...
                sum(ES(i,tabComp(tabComp.Ng=="IS",:).no),2)])/(b/mu);
        end
        plot(vecRho,prev_g)
        
        figure(2)
        U_g = vecRho.*(prev_g-c);
        plot(vecRho, U_g)
    end
    
    if(0)
        
        %% 2. Find the argmax of U
        addpath('C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/matlab/');
        addpath('C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/lib/');
        addpath('C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab');
        addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\AMPL_models')
        setupOnce;
        
        ampl = AMPL;
        ampl.read('C:/Users/Moi/Documents/IPLESP/These/Codes/AMPL_models/main_obj_g.mod'); %pour resoudre le probleme d'optim, cad trouver rho_hat
        
        up_bnd_alpha  = ampl.getParameter('up_bnd_alpha'); up_bnd_alpha.setValues(max([paramTab{1}.alpha,paramTab{2}.alpha]));
        %assigningParameters;
        ampl = assigningParametersToAMPL(paramTab,paramRho,mu,b,ampl,{'Ng'});
        
        knitro_options = 'outmode=1 ms_enable=1 ms_maxsolves=4 feastol=1e-5 maxtime_real=20 ms_maxtime_real=60 ms_outsub=1 ms_numthreads=4 ncvx_qcqp_init=0';
        ampl.setOption('solver', 'knitroampl');
        ampl.setOption('solver_msg','0');
        ampl.setOption('knitro_options', [knitro_options]);
        
        log_path='C:/Users/Moi/Documents/IPLESP/These/Codes/AMPL_models/knitro_out/';
        mkdir(log_path);
        vec_c = -0.1:0.01:0.12;
        i=0; vec_rho_max=zeros(length(vec_c),1); vec_all_status = [];
        for c_val = vec_c
            i=i+1;
            disp("/n######################################################/n");
            disp([i,vec_c(i)])
            %cost
            c = ampl.getParameter('c'); c.setValues(c_val);
            outdir=[log_path, 'kn_out_', num2str(c_val)];
            mkdir(outdir);
            ampl.setOption('knitro_options', [knitro_options, ' outdir=', outdir, ' outmode=2', ' outname=knitro.log']);
            ampl.solve();
            status = ampl.getValue("solve_result_num");
            if status ~= 0
                disp("/n#### Non-optimal status, check multi-start procedure. ###/n");
                break;
            end
            vec_all_status(i) = status;
            
            %rho_hat
            %rho = ampl.getVariable('rho_cg');
            rho = ampl.getVariable('rho_g');
            df = rho.getValues;
            a = df.val;
            vec_rho_max(i) = a{1};
        end
        if size(vec_all_status, 1) > 0
            figure()
            plot(vec_c,vec_rho_max)
        end
        ampl.close();
        
    end
    
    %%
    clear all; close all;
    
    b=100;
    [paramTab,mu,vecAlphas] = sampleParameters_v4(true,true,true,true,b,0.1);   %Ct,Ng,HIV,syph
    createParamRho;
    paramRho.rho_hsc = 0.1544;
    [ES,~,~,~,~,~] = ...
        P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,1,'ode45',20);
    
    kit={'HIV','syphilis','Ct'};
    [vecAlpha,ES] = findAlpha(paramTab,paramRho,mu,b,kit,1);
    
    vecAlpha
    P   = P_kit(ES,kit)
    P_h = P_kit(ES,{'HIV'})
    P_s = P_kit(ES,{'syphilis'})
    P_c = P_kit(ES,{'Ct'})
    %P_g = P_kit(ES,{'Ng'})
    
    
    
    
    
    
    
    
    %%
    vec_c = -1:0.01:2;
    
    clear vec_rho_hat;
    %vecAlpha = 0.3542467;
    
    [vec_rho_hat] = findRhohat_kit(paramTab,mu,b,paramRho,kit,vec_c,vecAlpha);
    
    figure()
    plot(vec_c,vec_rho_hat)
    hold on
    
    %%
    b=100;
    [paramTab,mu,vecAlphas] = sampleParameters_v4(true,true,true,true,b,0.1);   %Ct,Ng,HIV,syph
    createParamRho;
    
    aprioriBndsC.sup = 1; aprioriBndsC.inf = -1;
    %vecAlpha = findAlpha(paramTab,paramRho,mu,b,kit);
    clear tab_res cnn;
    [cnn, tab_res, msg] = find_Cnn_kit(paramTab,mu,b,paramRho,1,kit,aprioriBndsC,vecAlpha);
    
    struct2table(tab_res)
    figure()
    plot(tab_res.c,tab_res.rhohat,'-*')
    cnn
    
    figure()
    plot(tab_res.c,tab_res.kit, '-*','DisplayName','$\hat\rho_{hc}$')
    hold on
    plot(tab_res.c,tab_res.HIV, '-*', 'DisplayName','$\hat\rho_h$')
    plot(tab_res.c,tab_res.Ct, '-*','DisplayName','$\hat\rho_c$')
    legend('Interpreter','latex')
    
    figure(1)
    plot(vec_c_res,vec_rho_hat_res,'*')
    
    [paramTab,paramRho] = updateParamRho(paramTab,paramRho,kit,0.348);
    [ES,~,~,~,~,~] = ...
        P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,1,'knitroampl',20);
    P = P_kit(ES,kit)
    %--------------------------------------------------------------------------%
    
    %%
    vecP=[]; vecP_s=vecP;vecP_h=vecP;vecP_c=[];
    vecRHO=0.5:0.01:0.8;
    for rho=vecRHO
        [paramTab,paramRho] = updateParamRho(paramTab,paramRho,kit,rho);
        [ES,~,~,~,~,~] = ...
            P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,1,'knitroampl',15);
        vecP   = [vecP,P_kit(ES,kit)];
        vecP_s = [vecP_s,P_kit(ES,{'syphilis'})];
        vecP_h = [vecP_h,P_kit(ES,{'HIV'})];
        vecP_c = [vecP_c,P_kit(ES,{'Ct'})];
    end
    figure(2)
    plot(vecRHO,vecP(end-length(vecRHO)+1:end))
    figure(3)
    hold on
    c=-0.0%0.296875000000000;
    plot(vecRHO,vecRHO.*(vecP - c))
end


%% sensitivityAnalyis_v2.m; (voir le code)





