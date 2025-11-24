function [P,P12,PHIV,PIST,ES] = P12_SICTPSEIIS_v7(param1,param2,mu,b,paramRho_all,vecRho,f,solveMethod,optSolver)
%tolP0=10^(-5);
tolP0 = optSolver.tolP0;
verbose = optSolver.verbose;

%SICTP
theta1 = param1.theta0; betaI1 = param1.betaI;   betaC1 = param1.betaC;
sigma1 = param1.sigma; zeta1 = param1.zeta;  p1 = param1.p;

%SEIIS
beta2 = param2.beta; gamma2 = param2.gamma; nu2 = param2.nu; sigma2 = param2.sigma; eps2 = param2.eps;

paramRho.eta_X_prep = paramRho_all.(['eta_',param2.mini_d,'_prep']);
paramRho.eta_h_prep = paramRho_all.eta_h_prep;
paramRho.eta_X_art  = paramRho_all.(['eta_',param2.mini_d,'_art']);
paramRho.rho_h      = paramRho_all.rho_h;
paramRho.rho_X      = paramRho_all.(['rho_',param2.mini_d]);
paramRho.VTunderART = paramRho_all.VTunderART;

mySeed=123;
% cm = lines(20); k=1;
% figure()
% for p1=0:0.1:0.5
    i=1;
    P = zeros(length(vecRho),1);      P12=zeros(length(vecRho),1);
    PHIV = zeros(length(vecRho),1);   resteP=zeros(length(vecRho),1);
    PIST = zeros(length(vecRho),1);   ES = zeros(28,length(vecRho));
    for rho=vecRho
        paramRho.rho_hX = rho;
        %rho = min([RHO,alpha1,alpha2]);
        Y0 = ones(4*7,1)*2;
        if(isequal(solveMethod,'ode45'))
            tspan = optSolver.tspan;
            options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
            [~,Ys] = ode45(@(t,Y) ODE_SICTPSEIIS_v7(t,Y,betaI1,betaC1,theta1,sigma1,zeta1,p1,...
                beta2, gamma2, nu2, eps2, sigma2,...
                paramRho.rho_h,paramRho.rho_X,paramRho.rho_hX,...
                paramRho.eta_h_prep,paramRho.eta_X_prep,paramRho.eta_X_art,...
                paramRho.VTunderART,mu,b),...
                tspan,Y0,options);
            ES(:,i)=Ys(end,:);
            
            if min(ES(:,i))<-1
                warning('pb de convergence de ES')
            end
            
        elseif (isequal(solveMethod,'fsolve'))
            options = optimoptions('fsolve','Display','none','TolFun',optSolver.TolFun);
            ES(:,i) = fsolve(@(Y)  ODE_SICTPSEIIS_v7(1,Y,betaI1,betaC1,theta1,sigma1,zeta1,p1,...
                beta2, gamma2, nu2, eps2, sigma2,...
                paramRho.rho_h,paramRho.rho_X,paramRho.rho_hX,...
                paramRho.eta_h_prep,paramRho.eta_X_prep,paramRho.eta_X_art,...
                paramRho.VTunderART,mu,b),...
                Y0,options);
            
        elseif (isequal(solveMethod,'knitro-ampl'))
            
            addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
            log_path = 'C:\Users\Moi\Desktop\Temporaire\tests';
            modele = 'main_sictp_seiis.mod';
            
            ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
            ampl = AMPL;
            ampl.reset();
            ampl.cd(ampl_models_dir);
            ampl.read([ampl_models_dir, modele])
            
            %--- assigningParameters ---%;
            
            up_bnd_P0  = ampl.getParameter('bnd_sup_0');
            up_bnd_P0.setValues(tolP0);
            
            betaX = ampl.getParameter('betaX');   betaX.setValues([beta2]);
            gammaX = ampl.getParameter('gammaX'); gammaX.setValues([gamma2]);
            nuX = ampl.getParameter('nuX');       nuX.setValues([nu2]);
            epsX = ampl.getParameter('epsX');     epsX.setValues([eps2]);
            sigmaX = ampl.getParameter('sigmaX'); sigmaX.setValues([sigma2]);
            betaIh = ampl.getParameter('betaIh'); betaIh.setValues([betaI1]);
            betaCh = ampl.getParameter('betaCh'); betaCh.setValues([betaC1]);
            sigmah = ampl.getParameter('sigmah'); sigmah.setValues([sigma1]);
            thetah = ampl.getParameter('thetah'); thetah.setValues([theta1]);
            zetah = ampl.getParameter('zetah');   zetah.setValues([zeta1]);
            ph = ampl.getParameter('ph');         ph.setValues([p1]);
            
            %General parameters
            mu_apml = ampl.getParameter('mu'); mu_apml.setValues([mu]);
            b_ampl = ampl.getParameter('b');   b_ampl.setValues([b]);
            
            %Routine testing under PrEP/ART
            eta_X_prep = ampl.getParameter('eta_X_prep');     eta_X_prep.setValues([paramRho.eta_X_prep]);
            eta_h_prep = ampl.getParameter('eta_h_prep');     eta_h_prep.setValues([paramRho.eta_h_prep]);
            eta_X_art  = ampl.getParameter('eta_X_art');      eta_X_art.setValues([paramRho.eta_X_art]);
            vt_under_art  = ampl.getParameter('VTunderART');  vt_under_art.setValues([paramRho.VTunderART]);
            
            rho_h = ampl.getParameter('rho_h'); rho_h.setValues([paramRho.rho_h]);
            rho_X = ampl.getParameter('rho_X'); rho_X.setValues([paramRho.rho_X]);
            rho_hX   = ampl.getParameter('rho_hX'); rho_hX.setValues([paramRho.rho_hX]);

            % --- %
            
            [~,~]=mkdir(log_path);
            knitro_options = '';
            knitro_options = [knitro_options, 'outmode=2 ms_enable=1 ms_maxsolves=1 feastol=1e-7 maxtime_real=20 ms_maxtime_real=20 ms_outsub=1 ms_numthreads=4 ncvx_qcqp_init=0 ms_seed=',num2str(mySeed)];
            
            %store outputs
            outdir=[log_path];
            [~,~]=mkdir(outdir);
            ampl.setOption('knitro_options', [knitro_options, ' outdir=', outdir,' outname=knitro.log']);
            ampl.setOption('solver', 'knitroampl');
            
            if verbose
                ampl.solve();
            else
                output = evalc("ampl.solve()");
            end
            
            status = ampl.getValue("solve_result_num");
            if status ~= 0 && status ~=401
                disp("#### Non-optimal status, check multi-start procedure. ###");
                %ampl.close();
                %return;
                %break;
            end

            ES_ampl = ampl.getVariable('Y');
            df = ES_ampl.getValues;
            a = df.val;
            ES(1:28,i) = cell2mat(a);
            
            ampl.close();
        end
        
        
        %prevalence of asymptomatic and untreated:
        PHIV(i) = f*sum(ES([2,3,5,6,9,10,12,13,16,17,19,20],i))/(b/mu);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % retirer les gens sous prep !
        % VTunderART à ceux sous ART !
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        resteP(i) = sum(ES([8,11,14,15,18,21],i))/(b/mu);
        PIST(i) = sum(ES(8:28,i))/(b/mu);
        %PIST_prep(i) = sum(ES(7:28,i))/(b/mu);
        P12(i)  = PHIV(i)+resteP(i);
        % table2array(tabComp((tabComp.HIV1=="I" | tabComp.HIV1=="Ip" | ...
        % tabComp.HIV1=="C" | tabComp.HIV1=="Cp" | tabComp.STI1=="E" | ...
        % tabComp.STI1=="IA") & tabComp.STI1~="IS","no"))'
        P(i) = min(max(P12(i),0),1);
        %sousPrep(i) = sum(ES([2,3,5,6,9,10,12,13,16,17,19,20],i))/(b/mu);
        i=i+1;
    end
    
%     plot(vecRho,PHIV,'--','DisplayName',['$P_h$ p=',num2str(p1)],'Color',cm(k,:))
%     hold on
%     plot(vecRho,PIST,'-','DisplayName',['$P_c$ p=',num2str(p1)],'Color',cm(k,:))
%     legend('Interpreter','latex')
%     k=k+1;
% end


%     HIV1    STI1    no
%     ____    ____    __
%
%     "S"     "S"      1
%     "I"     "S"      2
%     "C"     "S"      3
%     "P"     "S"      4
%     "Ip"    "S"      5
%     "Cp"    "S"      6
%     "T"     "S"      7
%     "S"     "E"      8
%     "I"     "E"      9
%     "C"     "E"     10
%     "P"     "E"     11
%     "Ip"    "E"     12
%     "Cp"    "E"     13
%     "T"     "E"     14
%     "S"     "IA"    15
%     "I"     "IA"    16
%     "C"     "IA"    17
%     "P"     "IA"    18
%     "Ip"    "IA"    19
%     "Cp"    "IA"    20
%     "T"     "IA"    21
%     "S"     "IS"    22
%     "I"     "IS"    23
%     "C"     "IS"    24
%     "P"     "IS"    25
%     "Ip"    "IS"    26
%     "Cp"    "IS"    27
%     "T"     "IS"    28
end