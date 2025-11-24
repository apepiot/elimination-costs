function [vec_rho_hat,ES,status] = findRhohat_kit(paramTab,mu,b,paramRho,kit,vec_c,vecAlpha,verbose,mySeed,log_path,solver)

%log_path='C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/knitro_out/';
%tolP0 = 1e-3;
k = indexKit(kit);

if strcmp(solver,'knitroampl')
    warning('le systeme d ODE pour le modele a 4 infections n a pas du etre mis a jour')
    modele = ['main_obj_',k,'.mod'];
    
    ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
    ampl = AMPL;
    ampl.reset();
    ampl.cd(ampl_models_dir);
    ampl.read([ampl_models_dir, modele])
    
    up_bnd_alpha  = ampl.getParameter('up_bnd_alpha');
    up_bnd_alpha.setValues(max(vecAlpha));
    
    up_bnd_P0  = ampl.getParameter('bnd_sup_0');
    up_bnd_P0.setValues(tolP0);
    
    %assigningParameters;
    ampl = assigningParametersToAMPL(paramTab,paramRho,mu,b,ampl,kit);
    
    [~,~]=mkdir(log_path);
    
    knitro_options = '';
    knitro_options = [knitro_options, 'outmode=2 ms_enable=1 ms_maxsolves=4 feastol=1e-8 maxtime_real=20 ms_maxtime_real=60 ms_outsub=1 ms_numthreads=4 ncvx_qcqp_init=0 ms_seed=',num2str(mySeed)];
    
    i=0; vec_rho_hat=nan(length(vec_c),1); vec_all_status = ones(length(vec_c),1);
    ES = nan(560,length(vec_c));
    for c_val = vec_c
        i=i+1;
        
        %cost
        c = ampl.getParameter('c'); c.setValues(c_val);
        
        %store outputs
        outdir=[log_path, 'kn_out_', num2str(c_val),'/'];
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
            disp(['status=',num2str(status),' c=',num2str(c_val)]);
            %ampl.close();
            %return;
            break;
        end
        
        vec_all_status(i) = status;
        
        rho = ampl.getVariable(['rho_',k]);
        df = rho.getValues;
        a = df.val;
        vec_rho_hat(i) = a{1};
        
        ES_ampl = ampl.getVariable('Y');
        df = ES_ampl.getValues;
        a = df.val;
        ES(1:560,i) = cell2mat(a);
    end
    ampl.close();
    
elseif strcmp(solver,'matlab-fsolve')
elseif strcmp(solver,'matlab-ode45')
end
end

