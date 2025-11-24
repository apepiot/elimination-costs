function [P,ES,status] = P_mod_v7(paramTab,paramRho,b,mu,f,mod,verbose,paramSolver,log_path,ampl_models_dir)
tolP0         = paramSolver.tolP0;
nbRelancesTot = paramSolver.nbRelanceMax;
timeLimit     = paramSolver.timeLimit; %seconds

model   = ['mainP_',mod,'_v2.mod'];
sub_mod = mod;
%k=indexKit(kit);

relancerCalcul=1; 
nbRelances=0;

if ~isfolder(log_path)
    mkdir(log_path);
end

while relancerCalcul && nbRelances<nbRelancesTot
    disp([' Calcul de P_',sub_mod, ' avec le modele ',sub_mod])

    if nbRelances>0
        disp(['  Relance numero ',num2str(nbRelances)])
    end
    
    ampl = AMPL;
    ampl.reset();
    ampl.cd(ampl_models_dir);
    ampl.read([ampl_models_dir, model])
    outdir=[log_path, ['kn_out_P_',mod]];
    
    if ~isfolder(outdir)
        mkdir(outdir);
    end
    
    %assigningParameters;
    ampl = assigningParametersToAMPL_v3(paramTab,paramRho,mu,b,ampl,{'HIV','syphilis','Ct','Ng'});
        
    up_bnd_prevalence  = ampl.getParameter('rho_hscg');
    up_bnd_prevalence.setValues(paramRho.rho_hscg);
    up_bnd_prevalence  = ampl.getParameter('bnd_sup_0');
    up_bnd_prevalence.setValues(tolP0);
    
    
    %% Sub-model
    ps_0_ampl  = ampl.getParameter('p_s_0');
    ps_0_ampl.setValues(int8(contains(sub_mod,'s')));
    pc_0_ampl  = ampl.getParameter('p_c_0');
    pc_0_ampl.setValues(int8(contains(sub_mod,'c')));
    pg_0_ampl  = ampl.getParameter('p_g_0');
    pg_0_ampl.setValues(int8(contains(sub_mod,'g')));
    ph_0_ampl  = ampl.getParameter('p_h_0');
    ph_0_ampl.setValues(int8(contains(sub_mod,'h')));
    
    %%
    knitro_options = '';
    knitro_options = [knitro_options, 'outmode=2 ms_enable=1 ms_maxsolves=4 feastol_abs=1e-6 maxtime_real=',num2str(timeLimit),' ms_maxtime_real=60 ms_outsub=1 ms_numthreads=4 ncvx_qcqp_init=0 ms_seed=',num2str(nbRelances)];
    ampl.setOption('knitro_options', [knitro_options, ' outdir=', outdir,' outname=knitro.log']);
    ampl.setOption('solver', 'knitroampl');
    
    if verbose
        ampl.solve();
    else
        output = evalc("ampl.solve()");
    end
    
    status = ampl.getValue("solve_result_num");
    
    ES_ampl = ampl.getVariable('Y');
    df = ES_ampl.getValues;
    a = df.val;
    ES = cell2mat(a);
    
%     Pmod_ampl = ampl.getVariable(['Prevalence_',upper(mod)]);
%     df = Pmod_ampl.getValues;
%     a = df.val;
%     P.(mod) = cell2mat(a);
     
    Ph_ampl = ampl.getVariable(['Prevalence_H']);
    temp = Ph_ampl.getValues;
    a = temp.val;
    P.h = cell2mat(a);
    
    Ps_ampl = ampl.getVariable(['Prevalence_S']);
    temp = Ps_ampl.getValues;
    a = temp.val;
    P.s = cell2mat(a);
    
    Pc_ampl = ampl.getVariable(['Prevalence_C']);
    temp = Pc_ampl.getValues;
    a = temp.val;
    P.c = cell2mat(a);
    
    Pg_ampl = ampl.getVariable(['Prevalence_G']);
    temp = Pg_ampl.getValues;
    a = temp.val;
    P.g = cell2mat(a);
    
    %P.(k) = P_kit(ES,kit);
    
    ampl.close();
   
    
    %disp(['  P_',mod,'=',num2str(P.(mod))])
    
    %Pk = P_kit(ES,kit);
    relancerCalcul = 0;
    
    if P.h <= tolP0 && contains(sub_mod,'h')
        sub_mod = erase(sub_mod,'h');
        relancerCalcul = 1;
    end

    if relancerCalcul
        nbRelances = nbRelances+1;
    end
end

end

