function [rhohat,Cval,ES,msg] = fminU_knitro_v7(kit,paramTab,paramRho,mu,b,c,optSolver)
mySeed=123;
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
log_path = 'C:/Users/Moi/Desktop/Temporaire/tests';
ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';

k = indexKit(kit); %kit a tester
mod = k; % modele a faire tourner, a voir

%On minimise selon rho_k


% Ampl model to call
model = 'main_obj';
modelTypes = {'SICTP','SEIIIS','SEIIS'};
for i=1:length(paramTab)
    for j=1:length(modelTypes)
        thismodel = modelTypes{j};
        if strcmp(paramTab{i}.modelType,thismodel)
            model = [model,'_',thismodel];
        end
    end
end
model = [model,'.mod'];

%Name of the VT rate to maximize U
nameRho = k;
if contains(k,'c')
    nameRho = strrep(nameRho,'c','X');
    if contains(k,'g')
        nameRho = strrep(nameRho,'g','Y');
    end
else
    if contains(k,'g')
        nameRho = strrep(nameRho,'g','Y');
    end
end

ampl = AMPL;
ampl.reset();
ampl.cd(ampl_models_dir);
ampl.read([ampl_models_dir, model])

%--- assigningParameters ---%;
up_bnd_P0  = ampl.getParameter('bnd_sup_0');
up_bnd_P0.setValues(optSolver.tolP0);
up_bnd_alpha = ampl.getParameter('up_bnd_alpha');
up_bnd_alpha.setValues(optSolver.up_bnd_alpha);

ampl = assigningParametersToAMPL_v2(paramTab,paramRho,mu,b,ampl,k,mod);


% --- %

c_ampl  = ampl.getParameter('c');  c_ampl.setValues([c]);


[~,~]=mkdir(log_path);
knitro_options = '';
knitro_options = [knitro_options, 'outmode=2 ms_enable=1 ms_maxsolves=4 feastol=1e-6 maxtime_real=20 ms_maxtime_real=20 ms_outsub=1 ms_numthreads=4 ncvx_qcqp_init=0 ms_seed=',num2str(mySeed)];

%store outputs
outdir=[log_path];
[~,~]=mkdir(outdir);
ampl.setOption('knitro_options', [knitro_options, ' outdir=', outdir,' outname=knitro.log']);
ampl.setOption('solver', 'knitroampl');

if optSolver.verbose
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
ES = cell2mat(a);

Ph_ampl = ampl.getVariable(['Prevalence_H']);
temp = Ph_ampl.getValues;
a = temp.val;
Ph = cell2mat(a);

Pc_ampl = ampl.getVariable(['Prevalence_X']);
temp = Pc_ampl.getValues;
a = temp.val;
Pc = cell2mat(a);

rho_ampl = ampl.getVariable(['rho_',nameRho]);
temp = rho_ampl.getValues;
a = temp.val;
rhohat = cell2mat(a);

z = ampl.getObjective('Cost');
Cval = z.value;

msg = num2str(status);
ampl.close();

        
end