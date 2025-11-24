%% Test du systeme d'ODE (avec ode45)
clear all; 
close all;

%% 0. Parameters
%clear all;
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\')
b=100;

newSet=1; defaultSet=0;
if newSet 
    pHIV = 0;
    [paramTab,mu,~] = sampleParameters_v3_extent(true,true,true,true,b,pHIV);   %Ct,Ng,HIV,syph    
elseif defaultSet
    mu=1/35;
    paramTab{1}.beta = 1.3759;
    paramTab{1}.gamma = 24.8922;
    paramTab{1}.nu = 0.5246;
    paramTab{1}.eps = 0.0486;
    paramTab{1}.sigma=40.2979;
    paramTab{1}.rhob = 0;
    [~,~,paramTab{1}.alpha] = Rp_SEIIS_v4(paramTab{1}.beta,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,paramTab{1}.gamma,mu,b,0);

    paramTab{2}.beta=27.9007;
    paramTab{2}.gamma=34.9018;
    paramTab{2}.nu=1.2615;
    paramTab{2}.eps=0.9444;
    paramTab{2}.sigma=35.6416;
    paramTab{2}.rhob = 0;
    [~,~,paramTab{2}.alpha] = Rp_SEIIS_v4(paramTab{2}.beta,paramTab{2}.nu,paramTab{2}.eps,paramTab{2}.sigma,paramTab{2}.gamma,mu,b,0);

    paramTab{3}.betaI = 13.194799683714759;
    paramTab{3}.betaC = 1.406761238295781;
    paramTab{3}.sigma = 6.890538634865351;
    paramTab{3}.theta = 0.204544709146265; %%%%% verifier que c'est bien ça et pas gamma
    paramTab{3}.zeta = 0.631739291427463;
    paramTab{3}.eta = 4;
    paramTab{3}.p = 0%0.5%0.6000;
    paramTab{3}.rhob = 0;

    paramTab{4}.beta  = 0.9130;
    paramTab{4}.sigma = 9.4303;
    paramTab{4}.gamma1=0.;
    paramTab{4}.gamma3=0.0388;
    paramTab{4}.tau = 18.5060;
    paramTab{4}.theta=1.7839;
    paramTab{4}.nu = 0;
    paramTab{4}.rhob = 0;

end
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
paramRho.VunderART=1;

% Input parameters (all in one)
allParameters4d = {paramTab{3}.betaI,paramTab{3}.betaC,paramTab{3}.sigma,paramTab{3}.theta0,paramTab{3}.zeta,paramTab{3}.eta,paramTab{3}.p,...
    paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.gamma30,paramTab{4}.tau,paramTab{4}.theta,...
    paramTab{1}.beta,paramTab{1}.gamma,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,...
    paramTab{2}.beta,paramTab{2}.gamma,paramTab{2}.nu,paramTab{2}.eps,paramTab{2}.sigma,...
    [],mu,b,...
    paramTab{3}.rhob, paramTab{4}.rhob, paramTab{1}.rhob, paramTab{2}.rhob,...
    paramRho.rho_hs, paramRho.rho_hc,paramRho.rho_hg,...
    paramRho.rho_sc, paramRho.rho_sg,...
    paramRho.rho_cg,...
    paramRho.rho_hsc, paramRho.rho_hsg, paramRho.rho_hcg, paramRho.rho_scg,...
    paramRho.rho_hscg,...
    paramRho.eta_s_prep, paramRho.eta_c_prep, paramRho.eta_g_prep,...
    paramRho.eta_s_art, paramRho.eta_c_art, paramRho.eta_g_art, paramRho.VunderART};

%--------------------------------------------------------------------------%

%% Recreating tabComp
createTabComp;

%--------------------------------------------------------------------------%

if (1)
%% ODE system with the 4-disease model
addpath 'C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ODEsystems\'
Y0 = ones(560,1); tspan=[1,600];
[res4dis] = ode45( @(t,Y) ODE_SICTPSEIIISSEIIS2_v5_bis(t,Y,allParameters4d{:}),tspan, Y0);

%Total of the population
disp(['Pop. total_num - Pop tot theoric.:'])
disp(sum(res4dis.y(:,end)) - b/mu)

%% Comparaison avec les modèles à 1 infection
[ES_sictp,ES_seiiis,ES_ct,ES_ng] = singleDisES(paramTab,mu,b);
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tout ce qui est en dessous, n'a pas été vraiment testé et adapté à la
%%% version v5_bis de P1234, ODE...
%%% Ca vient de script_SEIISm_SICTP_SEIIIS_v5_bis.suite_tests.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



%% sensitivityAnalyis_v2.m; (voir le code)





