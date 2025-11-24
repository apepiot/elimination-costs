function [alpha,ES,msg] = findAlpha(paramTab,paramRho,mu,b,f,kit,verbose,log_path,paramSolver,ampl_models_dir)
%find alpha_kit, i.e. rho_kit such that Pkit(rho_kit)=0

paramHIV = getParams('h',paramTab);

if paramHIV.R_prep_0<1 && isequal(kit{:},'HIV')
    disp('HIV already eliminated')
    alpha = 0;
    ES = zeros(560,1);
    msg = '00';
else
    tolP0         = paramSolver.tolP0;
    maxBndAlpha   = paramSolver.maxBndAlpha;
    nbRelancesTot = paramSolver.nbRelanceMax;
    timeLimit     = paramSolver.timeLimit; %in seconds
    
    relancerCalculAlpha=1;
    k=indexKit(kit);
    nbRelances=0;
    output=0;
    
%     if any('s'==k)
%         timeLimit=30;
%     end
    
    modele = ['main_obj_P_alpha_',k,'.mod'];
    if ~isfolder(log_path)
        mkdir(log_path);
    end
    
    disp([' Calcul de alpha_',k])
    while relancerCalculAlpha && nbRelances<nbRelancesTot
        if nbRelances>0
            disp(['  Relance numero ',num2str(nbRelances)])
        end
        
        ampl = AMPL;
        ampl.reset();
        ampl.cd(ampl_models_dir);
        ampl.read([ampl_models_dir, modele])
        outdir=[log_path, ['kn_out_alpha_',k]];
        
        if ~isfolder(outdir)
            mkdir(outdir);
        end
        
        %assigningParameters;
        ampl = assigningParametersToAMPL(paramTab,paramRho,mu,b,ampl,kit);
        up_bnd_alpha  = ampl.getParameter('up_bnd_alpha');
        up_bnd_alpha.setValues(maxBndAlpha);
        
        up_bnd_prevalence  = ampl.getParameter('bnd_sup_0');
        up_bnd_prevalence.setValues(tolP0);
        
        knitro_options = '';
        knitro_options = [knitro_options, 'outmode=2 ms_enable=1 ms_maxsolves=6 feastol=1e-5 maxtime_real=',num2str(timeLimit),' ms_maxtime_real=60 ms_outsub=1 ms_numthreads=4 ncvx_qcqp_init=0 ms_seed=',num2str(nbRelances)];
        ampl.setOption('knitro_options', [knitro_options, ' outdir=', outdir,' outname=knitro.log']);
        ampl.setOption('solver', 'knitroampl');
        
        if verbose
            ampl.solve();
        else
            output = evalc("ampl.solve()");
        end
        
        status = ampl.getValue("solve_result_num");
        msg = num2str(status);
        
        rho = ampl.getVariable(['rho_',k]);
        df = rho.getValues;
        a = df.val;
        alpha = a{:};
        
        ES_ampl = ampl.getVariable('Y');
        df = ES_ampl.getValues;
        a = df.val;
        ES = cell2mat(a);
        
        ampl.close();
        
        disp(['  alpha_',k,'=',num2str(alpha)])
        
        Pk = P_kit(ES,kit);
        relancerCalculAlpha = 0;
        
        if abs(alpha-maxBndAlpha)<1e-3
            disp(['le alpha trouvé converge vers la borne superieure et Pk=',num2str(Pk)])
            relancerCalculAlpha=1;
            msg = [msg,'_-3'];
        end
        if Pk>=tolP0*1.05
            output
            disp(['  La prevalence du kit ',k,' evaluee en alpha_',k,'=',num2str(alpha),' n est pas nulle, Pk=',num2str(Pk)])
            msg = [msg, '_-1'];
            relancerCalculAlpha=1;
            if status==0 && abs(alpha-maxBndAlpha)<1e-3
                relancerCalculAlpha=0;
            end
        elseif status~=0 
            disp(['  status=',num2str(status)])
            msg = [msg, '_-2'];
            relancerCalculAlpha=1;
        end
        
        if nbRelances>=nbRelancesTot && relancerCalculAlpha
            disp(['  alpha_',k,' n a pas été trouvé.'])
            msg = [msg,'_-3',];
        end
        
        %if status==-401
        %   relancerCalculAlpha=0;
        %end
        
        if relancerCalculAlpha
            nbRelances = nbRelances+1;
        end
    end
end
end
