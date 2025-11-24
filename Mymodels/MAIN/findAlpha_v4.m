function [alpha,Palpha,ES,msg,elim_i] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir)
%find alpha_kit, i.e. rho_kit such that Pkit(rho_kit)=0
%pareil que v4, mais on initialise ampl qu'une fois
tolP0         = paramSolver.tolP0;
maxBndAlpha   = paramSolver.maxBndAlpha;
nbRelancesTot = paramSolver.nbRelanceMax;
timeLimit     = paramSolver.timeLimit; %seconds
nbIterMax     = paramSolver.iterMaxDicho;
tolAlpha      = paramSolver.tolAlpha;
msg           = '0';
% if checkFeas
%      [msgSolver,msgSol] = checkFeasability_v2(paramTab,paramRho,b,mu,f,kit,mod,verbose,paramSolver);
% end
%addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
% log_path = 'C:/Users/Moi/Desktop/Temporaire/tests';
% ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';

k = indexKit(kit);

%--- Ampl model to call ---%

% if contains(k,'h')
%     mod = k;
% else
%     mod = ['h',k];
% end


%Name of the VT rate to maximize U
nameRho = ['rho_',k];
%--------------------------%

relancerCalculAlpha=1;
nbRelances=0;
tabRecap_sim = [];

if ~isfolder(log_path)
    mkdir(log_path);
end

disp([' Recherche de alpha_',k, ' avec le modele ',mod])
if ~isequal(paramSolver.method_alpha,'dicho')
    model = ['mainFindAlpha_',mod,'-2.mod'];

    ampl = AMPL;
    ampl.reset();
    ampl.cd(ampl_models_dir);
    ampl.read([ampl_models_dir, model])
    outdir=[log_path, ['kn_out_alpha_',k]];
    
    if ~isfolder(outdir)
        mkdir(outdir);
    end
    
    %assigningParameters;
    ampl = assigningParametersToAMPL_v3(paramTab,paramRho,mu,b,ampl,kit);
    
    up_bnd_alpha  = ampl.getParameter('up_bnd_alpha');
    up_bnd_alpha.setValues(maxBndAlpha);
    
    up_bnd_prevalence  = ampl.getParameter('bnd_sup_0');
    up_bnd_prevalence.setValues(tolP0);
    
    %% Model
    ps_0_ampl  = ampl.getParameter('p_s_0');
    ps_0_ampl.setValues(int8(contains(mod,'s')));
    pc_0_ampl  = ampl.getParameter('p_c_0');
    pc_0_ampl.setValues(int8(contains(mod,'c')));
    pg_0_ampl  = ampl.getParameter('p_g_0');
    pg_0_ampl.setValues(int8(contains(mod,'g')));
    ph_0_ampl  = ampl.getParameter('p_h_0');
    ph_0_ampl.setValues(int8(contains(mod,'h')));
    
    while relancerCalculAlpha && nbRelances<nbRelancesTot
        if nbRelances>0
            disp(['  Relance numero ',num2str(nbRelances)])
        end
        
        knitro_options = '';
        knitro_options = [knitro_options, 'outmode=2 ms_enable=0 ms_maxsolves=6 feastol_abs=1e-6 maxtime_real=',num2str(timeLimit),' ms_maxtime_real=60 ms_outsub=1 ms_numthreads=4 ncvx_qcqp_init=0 ms_seed=',num2str(nbRelances)];
        ampl.setOption('knitro_options', [knitro_options, ' outdir=', outdir,' outname=knitro.log']);
        ampl.setOption('solver', 'knitroampl');
        
        %ampl.eval("write gmodelfile.nl;");
        if verbose
            ampl.solve();
        else
            output = evalc("ampl.solve()");
        end
        
        status = ampl.getValue("solve_result_num");
        msg = num2str(status);
        
        rho = ampl.getVariable([nameRho]);
        df = rho.getValues;
        a = df.val;
        alpha = a{:};
        
        ES_ampl = ampl.getVariable('Y');
        df = ES_ampl.getValues;
        a = df.val;
        ES = cell2mat(a);
        
        Pk_ampl = ampl.getVariable(['Prevalence_',upper(k)]);
        df = Pk_ampl.getValues;
        a = df.val;
        P.k = cell2mat(a);
        
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
        
        
        disp(['   ',nameRho,'=',num2str(alpha)])
        disp(['   P_',k,'=',num2str(P.k)])
        disp(['   P_s=',num2str(P.s)])
        disp(['   P_h=',num2str(P.h)])
        disp(['   P_c=',num2str(P.c)])
        disp(['   P_g=',num2str(P.g)])
        disp(['   status=',num2str(status)])
        
        relancerCalculAlpha = 0;
        if abs(alpha-maxBndAlpha)<1
            disp([' le alpha trouvé converge vers la borne superieure et Pk=',num2str(P.k)])
            relancerCalculAlpha=1;
            msg = [msg,'_-3'];
        end
        
        if nbRelances>=nbRelancesTot && relancerCalculAlpha
            disp(['  alpha_',k,' n a pas été trouvé.'])
            msg = [msg,'_-3',];
        end
        
        %Looking for the infection that has been eliminated first in the kit
        elim_i = [];
        if length(k) == 1
            elim_i = k;
        else
            for inf=k
                if abs(P.(inf)-tolP0)< tolP0 %10*tolP0
                    elim_i = [elim_i,inf];
                end
            end
        end
        
        if length(elim_i)~=1
            warning('il semblerait que alpha n a pas été bien trouvé')
            relancerCalculAlpha = 1;
        end
        if status~=0 && status ~=411
            relancerCalculAlpha = 1;
        end
        
        % Reverification
%         paramRho_temp = paramRho;
%         paramRho_temp.(['rho_',k]) = (1-tolAlpha)*alpha;
%         [Palpha,ES,status,ampl2] = P_mod_v8(paramTab,paramRho_temp,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,0);
%         ampl2.close();
%         disp([Palpha])
%         vecP=[];
%         for inf=k
%             vecP = [vecP,Palpha.(inf)];
%         end
%         
%         if sum(vecP<=tolP0)~=0
%             relancerCalculAlpha=1;
%         end
        
        if abs(alpha)<1e-8
            alpha=0;
        end
        
        tabRecap_sim(nbRelances+1,1) = alpha;
        tabRecap_sim(nbRelances+1,2:6) = table2array(struct2table(P));
        tabRecap_sim(nbRelances+1,7) = status;
        
        if relancerCalculAlpha
            nbRelances = nbRelances+1;
        end
    end
    ampl.close();
    
    tabRecap_sim2 = array2table(tabRecap_sim,'VariableNames',[convertCharsToStrings(nameRho),convertCharsToStrings(['P_',k]),"Ph","Ps","Pc","Pg","status"]);
    tabRecap_sim3 = sortrows(tabRecap_sim2,1);
    Palpha = P;
else
    tabRecap_sim3 = table('Size',[0,7],...
        'VariableNames',[convertCharsToStrings(nameRho),convertCharsToStrings(['P_',k]),"Ph","Ps","Pc","Pg","status"],...
        'VariableTypes',{'double','double','double','double','double','double','double'});
end
%%
% DICHOTOMIE
if nbRelances>=nbRelancesTot || isequal(paramSolver.method_alpha,'dicho')
    disp('  ...par dichotomie')
    affinerAlpha = 1;
    i = 0; nbRelances = 1;
    borne_sup = paramSolver.maxBndAlpha;
    borne_inf = 0;
    paramRho_temp = paramRho;
    new_alpha = 1;
    ampl = 0;
    
    while affinerAlpha && i<nbIterMax
        vecP = [];
        alpha_apriori = new_alpha;
        paramRho_temp.(['rho_',k]) = alpha_apriori;
        paramSolver.varToChange = ['rho_',k];
        disp(['    ',nameRho,'=',num2str(alpha_apriori)])
        
        [Palpha,ES,status,ampl] = P_mod_v8(paramTab,paramRho_temp,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl);
        %disp(Palpha)
        %PK = P_kit_v2(ES,kit);
        PK = P_kit_v3(ES,kit,f);
        tabRecap_sim3(nbRelances+i,:) = array2table(([alpha_apriori,PK,table2array(struct2table(Palpha)),status]),...
            'VariableNames',tabRecap_sim3.Properties.VariableNames);
        if status==0
            for inf=k
                vecP = [vecP,Palpha.(inf)];
            end
            if min(vecP)>tolP0
                new_alpha = (alpha_apriori+borne_sup)/2;
                borne_inf = alpha_apriori;
            else
                new_alpha = (alpha_apriori+borne_inf)/2;
                borne_sup = alpha_apriori;
            end
            i=i+1;
            if abs(new_alpha-alpha_apriori)<tolAlpha
                affinerAlpha = 0;
            end
            tabRecap_sim4 = sortrows(tabRecap_sim3,1);
            
            
            % Inconstitency check
            % for the infection j that has been considered eliminated first, we check
            % that there is a monotony in P_j around the rho of elimination,
            % should be decreasing
            tabRecap_sim5 = tabRecap_sim4(tabRecap_sim4.status==0,:);
            idx = find(tabRecap_sim5.(nameRho)==alpha_apriori);
            subTable_mono = table2array(tabRecap_sim5(unique([max(idx-1,1),idx,min(idx+1,size(tabRecap_sim5,1))]),{['P_',k],'Ph','Ps','Pc','Pg'}));
            monotonyRes = checkMonotonicity(subTable_mono,'decreasing',tolP0*10);
            if sum(monotonyRes)~=size(monotonyRes,2)
                warning(['probleme autour de rho=',num2str(alpha_apriori)])
                %si ça arrive au moment ou hiv est elimine c est normal,
                %car il peut y avoir un rebond des prevalences des autres
                %ist
                %sinon faut voir si on peut relancer ou mettre un msg d
                %erreur a ce alpha
            end
        else
            tabRecap_sim4 = sortrows(tabRecap_sim3,1);
            alpha_inf = tabRecap_sim4(tabRecap_sim4.(nameRho) < alpha_apriori,1);
            if isempty(alpha_inf)
                alpha_inf = table(0,'VariableNames',{nameRho});
            end
            alpha_sup = tabRecap_sim4(tabRecap_sim4.(nameRho) > alpha_apriori,1);
            
            if isempty(alpha_sup)
                new_alpha = alpha_apriori*1.1;
            else
                new_alpha = (alpha_sup.(nameRho)(1)-alpha_inf.(nameRho)(end))/(1+9*rand()) + alpha_inf.(nameRho)(end);
            end
            i=i+1;
            tabRecap_sim5 = tabRecap_sim4(tabRecap_sim4.status==0,:);
        end
    end
    ampl.close();
    
    %Looking for the infections that have been eliminated first in the kit
    idx = find(tabRecap_sim5.(nameRho)>=alpha_apriori);
    subTable = tabRecap_sim5(idx(1:2),[convertCharsToStrings(nameRho),"Ph","Ps","Pc","Pg"]);
    elim_i=[];
    for inf=k
        if (Palpha.(inf)-tolP0)<tolP0
            elim_i = [elim_i,inf];
        end
    end
    if isempty(elim_i)
        for inf=k
            if (subTable(2,:).(['P',inf])-tolP0)<tolP0
                elim_i = [elim_i,inf];
            end
        end
        alpha = subTable(2,:).(nameRho);
    else
        alpha=alpha_apriori;
    end
    
    if (affinerAlpha && i==nbIterMax)
        msg='2';
    end
    
    alpha_before_inf_stat = tabRecap_sim5(tabRecap_sim5.(nameRho)<alpha,"status");
    alpha_before_sup_stat = tabRecap_sim5(tabRecap_sim5.(nameRho)>alpha,"status");
    
    if isempty(alpha_before_inf_stat)
        alpha=0;
        msg='3';
    else
        if alpha_before_inf_stat.status(end)~=0 || alpha_before_sup_stat.status(1)~=0
            msg='1';
        end
    end
    disp('')
    disp(tabRecap_sim4);
end


disp('-----------------------------')




end
