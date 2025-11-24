function [rhohat,Cval,ES,P,msg,ampl] = fmaxU_knitro_v7_bis(kit,mod,paramTab,paramRho,mu,b,c,optSolver,ampl,newmod)
mySeed = optSolver.seed;
addpath 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';
log_path = 'C:/Users/Moi/Desktop/Temporaire/tests';
ampl_models_dir = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/AMPL_models/';

k = indexKit(kit); %kit a tester

% Ampl model to call
model = ['mainFindArgmaxU_',mod,'.mod'];

%Name of the VT rate to maximize U
nameRho = ['rho_',k];

if ampl==0  || newmod==1          %initializing
    if newmod && ampl~=0
         ampl.close();
    end
    ampl = AMPL;
    ampl.reset();
    ampl.cd(ampl_models_dir);
    ampl.read([ampl_models_dir, model])
    ampl = assigningParametersToAMPL_v3(paramTab,paramRho,mu,b,ampl,k);
end

%--- assigningParameters ---%;
ps_0_ampl  = ampl.getParameter('p_s_0');
ps_0_ampl.setValues(int8(contains(mod,'s')));
pc_0_ampl  = ampl.getParameter('p_c_0');
pc_0_ampl.setValues(int8(contains(mod,'c')));
pg_0_ampl  = ampl.getParameter('p_g_0');
pg_0_ampl.setValues(int8(contains(mod,'g')));
ph_0_ampl  = ampl.getParameter('p_h_0');
ph_0_ampl.setValues(int8(contains(mod,'h')));
up_bnd_P0  = ampl.getParameter('bnd_sup_0');
up_bnd_P0.setValues(optSolver.tolP0);
up_bnd_alpha = ampl.getParameter('up_bnd_alpha');
up_bnd_alpha.setValues(optSolver.up_bnd_alpha);
up_bnd_alpha = ampl.getParameter('inf_bnd_alpha');
up_bnd_alpha.setValues(optSolver.inf_bnd_alpha);

% --- %

c_ampl  = ampl.getParameter('c');  c_ampl.setValues([c]);


[~,~]=mkdir(log_path);
knitro_options = '';
knitro_options = [knitro_options, 'outmode=2 ms_enable=1 ms_maxsolves=4 feastol_abs=1e-6 maxtime_real=20 ms_maxtime_real=20 ms_outsub=1 ms_numthreads=4 ncvx_qcqp_init=0 ms_seed=',num2str(mySeed)];

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
    disp("   #### Non-optimal status, check multi-start procedure. ###");
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

rho_ampl = ampl.getVariable([nameRho]);
temp = rho_ampl.getValues;
a = temp.val;
rhohat = cell2mat(a);

z = ampl.getObjective('Cost');
Cval = z.value;

msg = num2str(status);
ampl.close();
ampl=0;        
end