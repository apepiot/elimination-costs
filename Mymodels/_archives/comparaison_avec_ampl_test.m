

%% --- Utility plotting --- %%
%clear all;
b=2;
paramTab{1}.beta = 0.5%1%1.3759;
paramTab{1}.gamma = 24.8922;
paramTab{1}.nu = 0.5246;
paramTab{1}.eps = 0.0486;
paramTab{1}.sigma=40.2979;

paramTab{2}.beta=27.9007;
paramTab{2}.gamma=34.9018;
paramTab{2}.nu=1.2615;
paramTab{2}.eps=0.9444;
paramTab{2}.sigma=35.6416;

paramTab{3}.betaI = 0.1%13.194799683714759;
paramTab{3}.betaC = 0.1%1.406761238295781;
paramTab{3}.sigma = 6.890538634865351;
paramTab{3}.theta = 0.204544709146265; %%%%% verifier que c'est bien ça et pas gamma
paramTab{3}.zeta = 0.631739291427463;
paramTab{3}.eta = 4;
paramTab{3}.p = 0.%0.6000;
paramTab{3}.rhob = 0;

paramTab{4}.beta  = 0.05%0.9130;
paramTab{4}.sigma = 9.4303;
paramTab{4}.gamma1=0.;
paramTab{4}.gamma3=0.0388;
paramTab{4}.tau = 18.5060;
paramTab{4}.theta=1.7839;
paramTab{4}.nu =0;

mu=1/35;

paramRho.eta_c_prep = 0; 
paramRho.eta_s_prep = 0;
paramRho.eta_h_prep = 0;
paramRho.rho_h = 0;
paramRho.rho_s = 0;
paramRho.rho_c = 0;
paramRho.rho_g = 0;
paramRho.rho_hs = 0;
paramRho.rho_hc = 0;
paramRho.rho_hg = 0;
paramRho.rho_sc = 0;
paramRho.rho_sg = 0;
paramRho.rho_cg = 0;
paramRho.rho_hsc = 0;
paramRho.rho_hsg = 0;
paramRho.rho_hcg = 0;
paramRho.rho_scg = 0;
paramRho.rho_hscg = 0;
paramRho.eta_s_prep = 0;
paramRho.eta_c_prep = 0;
paramRho.eta_g_prep = 0;
paramRho.eta_s_art = 0;
paramRho.eta_c_art = 0;
paramRho.eta_g_art = 0;

%Calul des rho_prime
addpath('C:/Users/Moi/Documents/IPLESP/These/Codes/Mymodels')
[paramTab{3}.R,~,paramTab{3}.alpha] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
        paramTab{3}.theta,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
        paramTab{3}.p,mu,b,paramTab{3}.rhob);
[paramTab{4}.R,~,paramTab{4}.alpha] = Rp_SEIIIS_v4(paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.tau,...
        paramTab{4}.nu,paramTab{4}.gamma1,paramTab{4}.theta,paramTab{4}.gamma3,mu,b,0);
[paramTab{1}.R,~,paramTab{1}.alpha] = Rp_SEIIS_v4(paramTab{1}.beta,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,paramTab{1}.gamma,mu,b,0);
[paramTab{2}.R,~,paramTab{2}.alpha] = Rp_SEIIS_v4(paramTab{2}.beta,paramTab{2}.nu,paramTab{2}.eps,paramTab{2}.sigma,paramTab{2}.gamma,mu,b,0);

%%
if (1)
c=0; f=1;
i=0; 
%vecU=[]; 
vecU_g=[];
vecRho=0:0.5:10;
Prev_g_seiis=[];
Prev_g_4=[];
for rho = vecRho
    paramRho.rho_cg = rho;
    i=i+1
    [~,~,~,~,ES] = U1234_SICTPSEIIISSEIIS2_v3(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,c,f,'ode45');
    Prev_c_et_g(i) = sum(ES(36:560))./(b/mu);
    vecU(i)= paramRho.rho_cg*(Prev_c_et_g(i)-c);
    Prev_g_4(i) = sum(ES(141:560))./(b/mu);
    
    [~,ES_g] = ode45(@(t,Y) ODE_SEIIS_v4(t,Y,paramTab{2}.beta,paramTab{2}.nu,paramTab{2}.gamma,paramTab{2}.sigma,paramTab{2}.eps,rho,mu,b), 1:150, ones(4,1));
    Prev_g_seiis(i) = sum(ES_g(end,2:4))/(b/mu);
    vecU_g(i) = rho*(Prev_g_seiis(i)-c);
end

%%
figure()
hold on
plot(vecRho,vecU)


[U] = U_SEIIS_v4(paramTab{2},mu,b,vecRho,c,f)
hold on;
plot(vecRho,vecU_g)

end


%% systeme d'ode simplifie sur gono (recherceh d'erreurs)
syms betas sigmas gamma3s taus thetas
syms betaX gammaX nuX epsX sigmaX
syms betaY gammaY nuY epsY sigmaY
syms betaIh betaCh sigmah thetah zetah eta_h_prep ph
syms rho_g
syms mu
syms b
paramRho.rho_g = rho_g;
syms Y [1 560]
Y([2:140,142:280,282:420,422:560]) =0
dY = ODE_SICTPSEIIISSEIIS2_v4(0,Y,betaIh,betaCh,sigmah,thetah,zetah,eta_h_prep,ph,...
                            betas,sigmas,gamma3s,taus,thetas,...
                            betaX,gammaX,nuX,epsX,sigmaX,...
                            betaY,gammaY,nuY,epsY,sigmaY,...
                            0,mu,b,...
                            paramRho.rho_h,paramRho.rho_s,paramRho.rho_c,paramRho.rho_g,...
                            paramRho.rho_hs,paramRho.rho_hc,paramRho.rho_hg,...
                            paramRho.rho_sc,paramRho.rho_sg,...
                            paramRho.rho_cg,...
                            paramRho.rho_hsc,paramRho.rho_hsg, paramRho.rho_hcg, paramRho.rho_scg,...
                            paramRho.rho_hscg,...
                            paramRho.eta_s_prep,paramRho.eta_c_prep,paramRho.eta_g_prep,...
                            paramRho.eta_s_art,paramRho.eta_c_art,paramRho.eta_g_art)
syms Y [560,1]
Y([2:3,5:140,142:143,145:280,282:283,285:420,422:423,425:560])=0
%% ---------------------
%clear all;
% Load from file the ampl model
addpath('C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/matlab/');
addpath('C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/lib/');
addpath('C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab');
setupOnce;

ampl = AMPL;
ampl.read('C:/Users/Moi/Documents/IPLESP/These/Codes/AMPL_models/U_4dis.mod');

% Reassign data - all instances
% HIV data
betaIh = ampl.getParameter('betaIh'); betaIh.setValues([paramTab{3}.betaI]);
betaCh = ampl.getParameter('betaCh'); betaCh.setValues([paramTab{3}.betaC]);
sigmah = ampl.getParameter('sigmah'); sigmah.setValues([paramTab{3}.sigma]);
thetah = ampl.getParameter('thetah'); thetah.setValues([paramTab{3}.theta]);
zetah = ampl.getParameter('zetah'); zetah.setValues([paramTab{3}.zeta]);
eta_h_prep = ampl.getParameter('eta_h_prep'); eta_h_prep.setValues([paramTab{3}.eta]);
ph = ampl.getParameter('ph'); ph.setValues([paramTab{3}.p]);

% syphilis data
betas = ampl.getParameter('betas'); betas.setValues([paramTab{4}.beta]);
sigmas = ampl.getParameter('sigmas'); sigmas.setValues([paramTab{4}.sigma]);
gamma3s = ampl.getParameter('gamma3s'); gamma3s.setValues([paramTab{4}.gamma3]);
taus = ampl.getParameter('taus'); taus.setValues([paramTab{4}.tau]);
thetas = ampl.getParameter('thetas'); thetas.setValues([paramTab{4}.theta]);

% Ct data
betaX = ampl.getParameter('betaX'); betaX.setValues([paramTab{1}.beta]);
gammaX = ampl.getParameter('gammaX'); gammaX.setValues([paramTab{1}.gamma]);
nuX = ampl.getParameter('nuX'); nuX.setValues([paramTab{1}.nu]);
epsX = ampl.getParameter('epsX'); epsX.setValues([paramTab{1}.eps]);
sigmaX = ampl.getParameter('sigmaX'); sigmaX.setValues([paramTab{1}.sigma]);

%Ng
betaY = ampl.getParameter('betaY'); betaY.setValues([paramTab{2}.beta]);
gammaY = ampl.getParameter('gammaY'); gammaY.setValues([paramTab{2}.gamma]);
nuY = ampl.getParameter('nuY'); nuY.setValues([paramTab{2}.nu]);
epsY = ampl.getParameter('epsY'); epsY.setValues([paramTab{2}.eps]);
sigmaY = ampl.getParameter('sigmaY'); sigmaY.setValues([paramTab{2}.sigma]);

%General parameters
mu = ampl.getParameter('mu'); mu.setValues([0.0286]);
b = ampl.getParameter('b'); b.setValues([100]);

%Testing parameters (baseline testing rate)
rho_h = ampl.getParameter('rho_h'); rho_h.setValues([paramTab{3}.rhob]);
rho_s = ampl.getParameter('rho_s'); rho_s.setValues([0]);
rho_c = ampl.getParameter('rho_c'); rho_c.setValues([0]);
rho_g = ampl.getParameter('rho_g'); rho_g.setValues([0]);

%Testing rates (kit)
rho_hs   = ampl.getParameter('rho_hs'); rho_hs.setValues([0]);
rho_hc   = ampl.getParameter('rho_hc'); rho_hc.setValues([0]);
rho_hg   = ampl.getParameter('rho_hg'); rho_hg.setValues([0]);
rho_sc   = ampl.getParameter('rho_sc'); rho_sc.setValues([0]);
rho_sg   = ampl.getParameter('rho_sg'); rho_sg.setValues([0]);
%rho_cg   = ampl.getParameter('rho_cg'); rho_cg.setValues([0]);
rho_hsc  = ampl.getParameter('rho_hsc'); rho_hsc.setValues([0]);
rho_hsg  = ampl.getParameter('rho_hsg'); rho_hsg.setValues([0]);
rho_hcg  = ampl.getParameter('rho_hcg'); rho_hcg.setValues([0]);
rho_scg  = ampl.getParameter('rho_scg'); rho_scg.setValues([0]);
rho_hscg = ampl.getParameter('rho_hscg'); rho_hscg.setValues([0]);

eta_s_prep = ampl.getParameter('eta_s_prep'); eta_s_prep.setValues([0]);
eta_c_prep = ampl.getParameter('eta_c_prep'); eta_c_prep.setValues([0]);
eta_g_prep = ampl.getParameter('eta_g_prep'); eta_g_prep.setValues([0]);
eta_s_art  = ampl.getParameter('eta_s_art'); eta_s_art.setValues([0]);
eta_c_art  = ampl.getParameter('eta_c_art'); eta_c_art.setValues([0]);
eta_g_art  = ampl.getParameter('eta_g_art'); eta_g_art.setValues([0]);
    
%up_bnd_alpha  = ampl.getParameter('up_bnd_alpha'); up_bnd_alpha.setValues(max([paramTab{1}.alpha,paramTab{2}.alpha,paramTab{3}.alpha,paramTab{4}.alpha])); %pas sure que ca marche tout le temps
up_bnd_alpha  = ampl.getParameter('up_bnd_alpha'); up_bnd_alpha.setValues(max([paramTab{1}.alpha,paramTab{2}.alpha])); 


ampl.setOption('solver', 'knitroampl');
%ampl.setOption('knitro_options', 'ms_enable=0 feastol=1e-6 maxtime_real=20 ncvx_qcqp_init=0');
%ampl.setOption('presolve', '0');
knitro_options = 'ms_enable=1 ms_maxsolves=10 feastol=1e-6 maxtime_real=20 ms_maxtime_real=120 ms_outsub=1 ms_numthreads=4 ncvx_qcqp_init=0';
log_path='C:/Users/Moi/Documents/IPLESP/These/Codes/AMPL_models/knitro_out';
mkdir(log_path);
vec_c = -1:0.01:1;
i=0; vec_rho_max=zeros(length(vec_c),1); vec_all_status = [];
for c_val = vec_c
    i=i+1;
    disp("/n###############################/n");
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
    rho = ampl.getVariable('rho_cg');
    df = rho.getValues;
    a = df.val;
    vec_rho_max(i) = a{1};
end
if size(vec_all_status, 1) > 0
    figure()
    plot(vec_c,vec_rho_max)
end
ampl.close();

