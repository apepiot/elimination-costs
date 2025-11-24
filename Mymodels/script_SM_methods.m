%clear all; 
close all;
addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
addpath('MAIN')
log_path = 'C:/Users/Moi/Desktop/Temporaire/tests/';
fig_path = 'C:\Users\Moi\Documents\IPLESP\These\Articles\Article 2\figures\';
ampl_models_dir = [pwd,'\MAIN\AMPL_models\'];
setupOnce;

set(groot,'defaultAxesTickLabelInterpreter','latex');  

%couleurs def
rougeHIV = [215, 0, 0]/255;
jauneS   = [250, 215, 0]/255;
bleuCt   = [56, 57, 186]/255;
vertNg   = [43, 152, 38]/255;

%colormap
% [WGB_log_cm] = subdivisedColormap([[1,1,1];[0.08,0.08,0.08];[0,0,0]], 100, 'log');
% size(parula)
% PARULA_log  = subdivisedColormap([[0.9769,0.9839,0.0805];[0.2422,0.1504,0.6603]],100,'log');
% PARULA_flip_quad  = subdivisedColormap(flipud(parula),5,'quad')
% map = flipud([0.1 0.1 0.1
%     0.15 0.1 0.2
%     0.2 0.1 0.7
%     0.1 0.5 0.8
%     0.2 0.7 0.6
%     0.8 0.7 0.3
%     0.9 1 0
%     1 1 1]);
% new_map_1  = subdivisedColormap(map(1:3,:),6,'quad');
% new_map_2  = subdivisedColormap(map(4:5,:),5,'quad');
% new_map_3  = subdivisedColormap(map(6:7,:),3,'quad');
% new_map = subdivisedColormap([new_map_1;new_map_2;new_map_3],5,'quad');
% 
% n1=500; n2=5;
% k_b = [linspace(0.1,0.2422,n1);linspace(0.1,0.1504,n1);linspace(0.1,0.6603,n1)]';
% y_w = [linspace(0.9769,1,n2);linspace(0.9839,1,n2);linspace(0.0805,1,n2)]';
% parula_w = [k_b;parula;y_w];
n=size(parula,1);
n1=500; n2=50;
k_b = [linspace(0.1,0.2422,n1);linspace(0.1,0.1504,n1);linspace(0.1,0.6603,n1)]';
y_w = [linspace(0.9769,0.985,n2);linspace(0.9839,0.985,n2);linspace(0.0805,0.985,n2)]';
parula_w = [k_b;parula;y_w;1,1,1];


%% Simple évolution de la prévalence (exemple 1, hscg)
%close all;
clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm
%----------------------------------%
pHIV    = 0.;%1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit     = {'HIV','syphilis','Ct','Ng'};
paramSolver.varToChange = 'rho_hscg';
vecRho  =  [0,0.05,0.1,0.11,0.12,0.13,0.14,0.2,0.3,0.35,0.39,0.4,0.5,0.6,0.65,0.7,1];
vecRho = sort(unique([vecRho,linspace(0,1,100)]));
fontSize = 18;
getParameters_avg; disp(num2str(ID_ech))
%----------------------------------%
%[alphas,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit,mod,0,log_path,paramSolver,ampl_models_dir);

paramTab{3}.p = pHIV;
Ph = [];Pc = [];Pg = [];Ps = [];Pk_r = zeros(1,length(vecRho));
ampl = 0;
i=1; paramRho_bis = paramRho;

for rho = vecRho
    disp([num2str(i),'/',num2str(length(vecRho))])
    disp(rho)
    paramRho_bis.(paramSolver.varToChange) = rho;
    [P,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl); 
    Ph(i) = P.h;
    disp(Ph(i));
    Pc(i) = P.c;
    Ps(i) = P.s;
    Pg(i) = P.g;
    [Pk_r(i),~,~] = P_kit_v3(ES,kit,f); %prevalence of asymptomatic and undiagnosed
    disp(Pk_r(i))
    i=i+1;
end
ampl.close();

fig = figure();
plot(vecRho,Pk_r,'k-','LineWidth',2)
hold on
plot(vecRho,Ph,'Color',rougeHIV,'LineWidth',2)
plot(vecRho,Ps,'Color',jauneS,'LineWidth',2)
plot(vecRho,Pc,'Color',bleuCt,'LineWidth',2)
plot(vecRho,Pg,'Color',vertNg,'LineWidth',2)
lgd = legend(['$\Pi_{',indexKit(kit),'}$'],'$\Pi_h$','$\Pi_s$','$\Pi_c$','$\Pi_g$','Interpreter','latex','FontSize',fontSize);
xlabel(['$\rho_{',indexKit(kit),'}$'],'FontSize',fontSize,'Interpreter','latex')
ylabel('Prevalence','Interpreter','latex','FontSize',fontSize)
xlim([0,0.7])
legend('boxoff') 
x1=18 ; x2=18;
lgd.ItemTokenSize = [x1,x2];

ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];

set(fig,'position',[100,100,400,340])
thisFigPath = [fig_path,'prevalence_',indexKit(kit),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

temp = [vecRho;Pk_r;Ph;Ps;Pc;Pg]'
temp=round(temp,5)


%%

figure()
vecC=-0.5:0.05:0.1;
for c=vecC
    U = vecRho.*(Pk_r-c);
    plot(vecRho,U,'DisplayName',num2str(c));
    hold on
end

%% Simple évolution de la prévalence (exemple 2, sc)
close all; clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm
%----------------------------------%
pHIV    = 0.;
f       = 1;
roundNo = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit     = {'HIV','Ct'};
paramSolver.varToChange = 'rho_hc';
fontSize = 18;
getParameters_avg; disp(num2str(ID_ech))
%[alphas,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit,mod,0,log_path,paramSolver,ampl_models_dir);
vecRho  = [0,0.05,0.1:0.005:0.2,0.4:0.005:0.42,1];
vecRho  = sort([0,0.01,0.05,0.08,0.1:0.005:0.2,0.4:0.01:0.67,0.655,0.645,0.68,1]);
%vecRho  = sort([0,0.01,0.05,0.08,0.1:0.005:0.2])%,0.4:0.025:0.67,0.68,1]);
vecRho = sort(unique([vecRho,linspace(0,1,100)]));

%----------------------------------%

Ph = [];Pc = [];Pg = [];Ps = [];Ph_p = [];Pk_r = zeros(1,length(vecRho));
ampl = 0;
i=1; paramRho_bis = paramRho;

for rho = vecRho
    paramRho_bis.(paramSolver.varToChange) = rho;
    [P,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl); 
    Ph(i) = P.h;
    Pc(i) = P.c;
    Ps(i) = P.s;
    Pg(i) = P.g;
    [Pk_r(i),~,~] = P_kit_v3(ES,kit,f);
    disp([num2str(i),'/',num2str(length(vecRho))])
    i=i+1;
end
ampl.close();

fig = figure();
plot(vecRho,Pk_r,'k-','LineWidth',2)
hold on
plot(vecRho,Ph,'Color',rougeHIV,'LineWidth',2)
plot(vecRho,Ps,'Color',jauneS,'LineWidth',2)
plot(vecRho,Pc,'Color',bleuCt,'LineWidth',2)
plot(vecRho,Pg,'Color',vertNg,'LineWidth',2)
lgd = legend(['$\Pi_{',indexKit(kit),'}$'],'$\Pi_h$','$\Pi_s$','$\Pi_c$','$\Pi_g$','Interpreter','latex','FontSize',fontSize);
xlabel(['$\rho_{',indexKit(kit),'}$'],'FontSize',fontSize,'Interpreter','latex')
ylabel('Prevalence','Interpreter','latex','FontSize',fontSize)
xlim([0,max(vecRho)])
legend('boxoff') 
x1=18 ; x2=18;
lgd.ItemTokenSize = [x1,x2];

ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];

set(fig,'position',[100,100,390,330])
thisFigPath = [fig_path,'prevalence_',indexKit(kit),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

temp = [vecRho;Pk_r;Ph;Ps;Pc;Pg]'
temp=round(temp,5);

%% Evolution de la prévalence percue
close all; clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm parula_w
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 3;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit     = {'HIV','syphilis'};
paramSolver.varToChange = ['rho_',indexKit(kit)];
vecRho  = [0:0.01:0.38,0.385:0.001:0.40, 0.41:0.01:0.66,0.661:0.001:0.68,1];
fontSize = 18;
getParameters_avg; disp(num2str(ID_ech))
%----------------------------------%

Pk_r = zeros(1,length(vecRho)); Pk_p = Pk_r;
Ph = [];Pc = [];Pg = [];Ps = [];Ph_p = [];
ampl = 0;
i=1;
paramRho_bis = paramRho;
for rho = vecRho
    paramRho_bis.(paramSolver.varToChange) = rho;
    [P,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl); 
    Ph(i) = P.h;
    Pc(i) = P.c;
    Ps(i) = P.s;
    Pg(i) = P.g;
    [Pk_r(i),Pk_p(i),Ph_p(i)] = P_kit_v3(ES,kit,f);
    disp([num2str(i),'/',num2str(length(vecRho))])
    i=i+1;
end
ampl.close();

%% (prevalence)
fig = figure();
plot(vecRho,Pk_r,'k-','LineWidth',2)
hold on
plot(vecRho,Pk_p,'k--','LineWidth',2)
plot(vecRho,Ph,'Color',rougeHIV,'LineWidth',2)
plot(vecRho,Ph_p,'Color',rougeHIV,'LineWidth',2,'LineStyle','--')
plot(vecRho,Ps,'Color',jauneS,'LineWidth',2)
plot(vecRho,Pc,'Color',bleuCt,'LineWidth',2)
plot(vecRho,Pg,'Color',vertNg,'LineWidth',2)
lgd = legend(['$\Pi_{',indexKit(kit),'}$'],['$\Pi_{',indexKit(kit),'}^b$'],'$\Pi_h$','$\Pi_{h}^b$','$\Pi_s$','$\Pi_c$','$\Pi_g$',...
    'Interpreter','latex','FontSize',fontSize);
xlabel(['$\rho_{',indexKit(kit),'}$'],'FontSize',fontSize,'Interpreter','latex')
ylabel('Prevalence','Interpreter','latex','FontSize',fontSize)
xlim([0,max(vecRho)])
legend('boxoff')
lgd.NumColumns = 2;
x1=18 ; x2=18;
lgd.ItemTokenSize = [x1,x2];

ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];

set(fig,'position',[100,100,390,330])
thisFigPath = [fig_path,'prevalence_',indexKit(kit),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

%% Prevalence (surf/3D graphs)
close all; clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm new_map parula_w
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit1    = {'syphilis','Ct'}; k_mod1='hscg_hscg';
kit2    = {'HIV','Ng'}; k_mod2 = 'hg_hg';
getParameters_avg; disp(num2str(ID_ech))
paramSolver.tolAlpha=1e-4; 
paramSolver.tolP0 = 5e-5; %%default:0.5e-4
[alphas,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit2,k_mod2,0,log_path,paramSolver,ampl_models_dir);
a2 = alphas.k;
[alphas,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit1,k_mod1,0,log_path,paramSolver,ampl_models_dir);
a1 = alphas.k; %quand rho_hg=0
%Récupérer les valeurs deja calculees avec load("var_Phg_sc.mat")
vecRho1 = sort([linspace(0,0.35,7),0.38,0.39,0.4,0.41,0.5,0.7,1]);%[linspace(0,a1-0.01,10),a1-1e-4,a1+1e-4,linspace(a1+0.01,1,3)];
vecRho2 = [linspace(0,a2-1.2e-3,12),linspace(a2-1e-3,a2+1e-3,9),linspace(a2+1.1*1e-3,1,3)];
%----------------------------------%
n = length(vecRho1);
m = length(vecRho2);
P1 = zeros(n,m);
P2 = zeros(n,m);
Ph  = zeros(n,m);
Ps  = zeros(n,m);
Pc  = zeros(n,m);
Pg  = zeros(n,m);
paramRho_bis = paramRho;
for i=1:n
    rho1 = vecRho1(i);
    %if rho1<=0.5 %sinon, 0 par defaut
        disp([num2str(i),'/',num2str(n)]); tic
        for j=1:m
            rho2 = vecRho2(j);
            paramRho_bis.(['rho_',indexKit(kit1)]) = rho1;
            disp(rho1)
            paramRho_bis.(['rho_',indexKit(kit2)]) = rho2;
            disp(rho2)
            ampl=0;
            [~,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,'hscg_hscg',verbose,paramSolver,log_path,ampl_models_dir,ampl); 
            ampl.close();
            P1(i,j) = P_kit_v3(ES,kit1,f);
            
            ampl=0;
            [~,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,'hg_hg',verbose,paramSolver,log_path,ampl_models_dir,ampl); 
            ampl.close();
            P2(i,j) = P_kit_v3(ES,kit2,f);
            Ps(i,j)  = P_kit_v3(ES,{'syphilis'},f);
            Pc(i,j)  = P_kit_v3(ES,{'Ct'},f);
            Ph(i,j)  = P_kit_v3(ES,{'HIV'},f);
            Pg(i,j)  = P_kit_v3(ES,{'Ng'},f);
        end
    %end
    toc
end
%save("var_Phg_sc04.mat")
%%
%load("var_Phg_sc0.mat")
%%
close all
idx_a2 = max(find(vecRho2<a2));
fontSize = 18;

%Colormap
colormap(flipud(parula_w))
cmin = 0.00001;
cmax = 0.9;

% Figure 1 : Psc
fig = figure();
subplot(1,2,1)
s1=surf(vecRho2(1:idx_a2),vecRho1,P1(:,1:idx_a2));
s1.EdgeColor = 'none';
s1.FaceColor = 'interp';
shading interp
hold on
s2=surf(vecRho2(idx_a2+1:end),vecRho1,P1(:,idx_a2+1:end));
s2.EdgeColor = 'none';
s2.FaceColor = 'interp';
shading interp
view(2)

ax = gca;
ax.FontSize = 0.6*fontSize; 
xticks(0:0.2:1)
yticks(0:0.2:1)
xlabel(['$\rho_{',indexKit(kit2),'}$'],'interpreter','latex','FontSize',1.2*fontSize)
ylabel(['$\rho_{',indexKit(kit1),'}$'],'interpreter','latex','FontSize',1.2*fontSize)
title(['$\Pi_{', indexKit(kit1),'}(','\rho_{',indexKit(kit2),'},','\rho_{',indexKit(kit1),'})$'],'interpreter','latex','FontSize',fontSize)
colorbar;
caxis([cmin, cmax]);

%xline(a2,'k:')

ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = 1.2*fontSize;
ax.YLabel.FontSize = 1.2*fontSize;
ax.Title.FontSize = 1.2*fontSize;
ax.YAxis.TickLength = [0,0];

% Figure 2 : Phg
subplot(1,2,2)
s=surf(vecRho2,vecRho1,P2);
s.EdgeColor = 'none';
s.FaceColor = 'interp';
shading interp
ax = gca;
ax.FontSize = 0.6*fontSize; 
xticks(0:0.2:1)
yticks(0:0.2:1)
xlabel(['$\rho_{',indexKit(kit2),'}$'],'interpreter','latex','FontSize',1.2*fontSize)
ylabel(['$\rho_{',indexKit(kit1),'}$'],'interpreter','latex','FontSize',1.2*fontSize)
title(['$\Pi_{', indexKit(kit2),'}(','\rho_{',indexKit(kit2),'},','\rho_{',indexKit(kit1),'})$'],'interpreter','latex','FontSize',fontSize)
colorbar;
caxis([cmin, cmax]);
view(2)

ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = 1.2*fontSize;
ax.YLabel.FontSize = 1.2*fontSize;
ax.Title.FontSize = 1.2*fontSize;
ax.YAxis.TickLength = [0,0];

set(fig,'position',[100,100,850,300])
colormap(flipud(parula_w))
thisFigPath = [fig_path,'prev_',indexKit(kit1),'_',indexKit(kit2),'_strat_',indexKit(kit1),'_',indexKit(kit2),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'col.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)
thisFigPath = [fig_path,'prev_',indexKit(kit1),'_',indexKit(kit2),'_strat_',indexKit(kit1),'_',indexKit(kit2),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'col.fig'];
%exportgraphics(fig,thisFigPath,'Resolution',300)


%% autre représentation

figure(3)
x = reshape(repmat(vecRho2,[length(vecRho1),1]),[1,length(vecRho2)*length(vecRho1)]);
y = repmat(vecRho1,[1,length(vecRho2)]);
z = reshape(P2,1,length(vecRho2)*length(vecRho1));
scatter3(x,y,z,12,1-z'./max(z)*[1,1,1])
view(2)

%cla
%patch([x nan],[y nan],[z nan],[z nan],'EdgeColor','interp','FaceColor','none')

%% utilité à la suite de ce qui précède (Pk_real utilise i.e. f=1)
%----------
f=1;
c=0.;
%----------
U1 = vecRho1'.*(P1-c);
U2 = vecRho2.*(P2-c);

cmin = 0.001;
cmax = max([U1(:)']);

fig1=figure(3);
s=surf(vecRho2,vecRho1,U1);
s.EdgeColor = 'none';
xlabel(['$\rho_{',indexKit(kit2),'}$'],'interpreter','latex','FontSize',fontSize)
ylabel(['$\rho_{',indexKit(kit1),'}$'],'interpreter','latex','FontSize',fontSize)
title(['$U_{', indexKit(kit1),'}(','\rho_{',indexKit(kit2),'},','\rho_{',indexKit(kit1),'})$'],'interpreter','latex','FontSize',fontSize)
colorbar;
view(2)
colormap(flipud(gray));
caxis([cmin, cmax]);
%
fig2=figure(4);
cmin = 0.001;
cmax = max([U2(:)']);
s=surf(vecRho2,vecRho1,U2);
s.EdgeColor = 'none';
xlabel(['$\rho_{',indexKit(kit2),'}$'],'interpreter','latex','FontSize',fontSize)
ylabel(['$\rho_{',indexKit(kit1),'}$'],'interpreter','latex','FontSize',fontSize)
title(['$U_{', indexKit(kit2),'}(','\rho_{',indexKit(kit2),'},','\rho_{',indexKit(kit1),'})$'],'interpreter','latex','FontSize',fontSize)
colorbar;
view(2)
colormap(flipud(gray));
%colormap(WGB_log_cm);
caxis([cmin, cmax]);

set(fig1,'position',[100,100,400,330])
%shading interp
thisFigPath = [fig_path,'utility_',indexKit(kit1),'_strat_',indexKit(kit1),'_',indexKit(kit2),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_c_',num2str(c),'_',num2str(ID_ech),'.pdf'];
exportgraphics(fig1,thisFigPath,'Resolution',300)

set(fig2,'position',[100,100,400,330])
%shading interp
thisFigPath = [fig_path,'utility_',indexKit(kit2),'_strat_',indexKit(kit1),'_',indexKit(kit2),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_c_',num2str(c),'_',num2str(ID_ech),'.pdf'];
% exportgraphics(fig2,thisFigPath,'Resolution',300)

%% Prevalences h,s,c,g as a function of rho_h and rho_s (surf/3D graphs) 
close all; clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm parula_w
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit1    = {'HIV'}; k_mod1 = 'hscg_hscg'; %{'Ct'} 
kit2    = {'syphilis'}; k_mod2 = 'hscg_hscg';%{'Ng'}
getParameters_avg; disp(num2str(ID_ech))
%----------------------------------%
[alphas,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit1,k_mod1,0,log_path,paramSolver,ampl_models_dir);
a1 = alphas.k; 
[alphas,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit2,k_mod2,0,log_path,paramSolver,ampl_models_dir);
a2 = alphas.k; 
%Récupérer les valeurs deja calculees avec load("var_Phg_sc.mat")
%vecRho2 = [linspace(0,0.5,20),linspace(0.52,1,5)];
%vecRho1 = [linspace(0,a1-1e-2,20),linspace(a1-1e-3,a1-1e-4,3),...
%   linspace(a1-5e-5,a1+5e-5,3), linspace(a1+1e-4,1,10)];
vecRho1 = [linspace(0,a1-1.1*1e-2,5), linspace(a1-1e-2,a1+1e-2,11), linspace(a1+1.1*1e-2,1,5)];%[linspace(0,a1-1.1*1e-2,5),a1-1e-3, a1-1e-4,a1+1e-4,a1+1e-3, linspace(a1+1.1*1e-2,1,5)];
vecRho2 = [linspace(0,a2-1.1*1e-2,5), linspace(a2-1e-2,a2+1e-2,5), linspace(a2+1.1*1e-2,1,5)];

n = length(vecRho1);
m = length(vecRho2);
P.h  = zeros(n,m);
P.s  = zeros(n,m);
P.c  = zeros(n,m);
P.g  = zeros(n,m);

disp(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'))

paramRho_bis = paramRho;
for i=1:n
    rho1 = vecRho1(i);
    disp([num2str(i),'/',num2str(n)]); tic
    for j=1:m
        rho2 = vecRho2(j);
        paramRho_bis.(['rho_',indexKit(kit2)]) = rho2;
        paramRho_bis.(['rho_',indexKit(kit1)]) = rho1;
        ampl=0;
        [~,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl); 
        ampl.close();
        P.s(i,j)  = P_kit_v3(ES,{'syphilis'},f);
        P.c(i,j)  = P_kit_v3(ES,{'Ct'},f);
        P.h(i,j)  = P_kit_v3(ES,{'HIV'},f);
        P.g(i,j)  = P_kit_v3(ES,{'Ng'},f);
    end
    toc
end
disp(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'))

%load("var_Ph_s_c_g_hs0.mat")
warning('pas sauvegarde')
%%
%load("var_Ph_s_c_g_hs0.mat")
close all;
fontSize = 18;
cmin = 0.00001;%(0.005);
cmax = 0.9;%(max([P.h(:)',P.s(:)',P.c(:)',P.g(:)']));
k=0; %letters = {'A','B','C','D'};
for inf={'h','s','c','g'}
    k=k+1;
    d = inf{:};
    fig = figure(k)
    Pi = (P.(d))';
    s=surf(vecRho1,vecRho2,Pi);
    %s=surf(vecRho1(vecRho1<a1),vecRho2,Pi(:,vecRho1<a1));
    s.EdgeColor = 'none';
    ax = gca;
    ax.FontSize = 0.6*fontSize; 
    %text(0.07,0.85,0.1,letters{k},'Interpreter','latex','FontSize',18)
    xlabel(['$\rho_{',indexKit(kit1),'}$'],'interpreter','latex','FontSize',1.2*fontSize)
    ylabel(['$\rho_{',indexKit(kit2),'}$'],'interpreter','latex','FontSize',1.2*fontSize)
    title(['$\Pi_{',d,'}(','\rho_{',indexKit(kit1),'},','\rho_{',indexKit(kit2),'})$'],'interpreter','latex','FontSize',1.2*fontSize)
    cbh = colorbar ; %Create Colorbar
    %cbh.Ticks = log(linspace(cmin, cmax, 8)) ; %Create 8 ticks from zero to 1
    %cbh.TickLabels = num2cell(exp(cbh.Ticks));    %Replace the labels of th
    xticks(0:0.2:1)
    view(2)
    colormap(flipud(parula_w))
    caxis([cmin, cmax]);
    set(fig,'position',[100,100,400,320])
    shading interp
    thisFigPath = [fig_path,'prev_',d,'_strat_',indexKit(kit1),'_',indexKit(kit2),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'col.pdf'];
    warning('pas de sauvegarde')
    %exportgraphics(fig,thisFigPath,'Resolution',300)
end

%% Figure représentant la prévalence du kit hscg, recursif
close all; clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit     = {'HIV','syphilis','Ct','Ng'};
paramSolver.varToChange = 'rho_hscg';
fontSize = 18;
getParameters_avg; disp(num2str(ID_ech))
%----------------------------------%
paramSolver.tolP0 = 0.5e-4;

kit1 = {'HIV','syphilis','Ct','Ng'}; k_mod1 = 'hscg_hscg';
[a1,P1,~,~,elim1] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit1,k_mod1,verbose,log_path,paramSolver,ampl_models_dir);
%
kit2 = {'HIV','syphilis','Ng'}; k_mod2 = 'hsg_hsg';
[a2,P2,~,~,elim2] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit2,k_mod2,verbose,log_path,paramSolver,ampl_models_dir);
%
kit3 = {'HIV','syphilis'}; k_mod3 = 'hs_hs';
[a3,P3,~,~,elim3] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit3,k_mod3,verbose,log_path,paramSolver,ampl_models_dir);
%
kit4 = {'HIV'}; k_mod4 = 'h_h';
[a4,P4,~,~,elim4] = findAlpha_v4(paramTab,paramRho,mu,b,f,kit4,k_mod4,verbose,log_path,paramSolver,ampl_models_dir);

%
vecRho  = sort(unique([0:0.01:1,reshape([a1;a2;a3;a4]*[0.99,1,1.01],1,12),0.8]));

Ph = [];Pc = [];Pg = [];Ps = [];Ph_p = [];Pk_r = zeros(1,length(vecRho));
i=1; 
ampl = 0;
for rho = vecRho
    paramRho_bis = paramRho;
    paramRho_bis.(paramSolver.varToChange) = rho;
    [P,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl); 
    Ph(i) = P.h;
    Pc(i) = P.c;
    Ps(i) = P.s;
    Pg(i) = P.g;
    [Pk_r(i),~,~] = P_kit_v3(ES,kit,f);
    disp([num2str(i),'/',num2str(length(vecRho))])
    i=i+1;
end
ampl.close();

%%
close all
fig = figure();
plot(vecRho,Pk_r,'LineWidth',3)
hold on

i1 = max(find(vecRho<=a1));
i2 = max(find(vecRho<=a2));
i3 = max(find(vecRho<=a3));
i4 = max(find(vecRho<=a4));

xline(a1,':','LineWidth',1,'Color',[0,0,0])
xline(a2,':','LineWidth',1,'Color',[0,0,0])
xline(a3,':','LineWidth',1,'Color',[0,0,0])
xline(a4,':','LineWidth',1,'Color',[0,0,0])

%plot(vecRho(1:i1),Pk_r(1:i1),'LineWidth',3)
plot(vecRho(i1:i2),Pk_r(i1:i2),'LineWidth',3)
plot(vecRho(i2:i3),Pk_r(i2:i3),'LineWidth',3)
plot(vecRho(i3:i4),Pk_r(i3:i4),'LineWidth',3)
plot(vecRho(i4:end),Pk_r(i4:end),'LineWidth',3)


ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];

title(['$\Pi_{',indexKit(kit1),'}(\rho_{',indexKit(kit1),'})$'],'Interpreter','latex','FontSize',fontSize)
xlabel(['$\rho_{',indexKit(kit1),'}$'],'FontSize',fontSize,'Interpreter','latex')
ylabel('Prevalence','Interpreter','latex','FontSize',fontSize)
xlim([0,0.8])
set(fig,'position',[100,100,400,340])
thisFigPath = [fig_path,'prevalence_',indexKit(kit1),'_rec_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

%% (utilite, f=1)
%------------------%
close all
f = 1; %ne pas changer par defaut
c = -0.2;
in_color = 0;
%------------------%
Uk = vecRho.*(Pk_r-c);
fig = figure();

if (in_color)
    plot(vecRho,Uk,'LineWidth',3)
    hold on
    %plot(vecRho(1:i1),Uk(1:i1),'LineWidth',3)
    plot(vecRho(i1:i2),Uk(i1:i2),'LineWidth',3)
    plot(vecRho(i2:i3),Uk(i2:i3),'LineWidth',3)
    plot(vecRho(i3:i4),Uk(i3:i4),'LineWidth',3)
    plot(vecRho(i4:end),Uk(i4:end),'LineWidth',3)
else
    plot(vecRho,Uk,'k-','LineWidth',3)
end

xline(a1,':','LineWidth',1,'Color',[0,0,0])
xline(a2,':','LineWidth',1,'Color',[0,0,0])
xline(a3,':','LineWidth',1,'Color',[0,0,0])
xline(a4,':','LineWidth',1,'Color',[0,0,0])

title(['${\rm U}_{',indexKit(kit1),'}(\rho_{',indexKit(kit1),'})$'],'Interpreter','latex','FontSize',fontSize)
xlabel(['$\rho_{',indexKit(kit1),'}$'],'FontSize',fontSize,'Interpreter','latex')
ylabel('Utility','Interpreter','latex','FontSize',fontSize)
xlim([0,0.7])
ylim([0,max(Uk)*1.01])
ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];
set(fig,'position',[100,100,400,350])
thisFigPath = [fig_path,'utility_',indexKit(kit1),'_rec_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'c_',num2str(c),'_',num2str(ID_ech),'.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

%% Evolution de rho_sc' en fonction de rho_hg
close all; clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit1  = {'HIV','Ng'};
kit2     = {'syphilis','Ct'};
paramSolver.varToChange = 'rho_hg';
fontSize = 18;
getParameters_avg; disp(num2str(ID_ech))
%----------------------------------%
vecAlpha=[];
paramRho_bis = paramRho;
paramSolver.tolP0 = 1e-4;

vecRho1 = [linspace(0,0.645,5),0.655:0.003:0.67,linspace(0.68,1,5)]; %pour rho_sc'(rho_hg)
%
%vecRhobis = 0.62;
%vecRho1 = [vecRho1,vecRhobis]
for rho1=vecRho1
    disp(rho1)
    paramRho_bis.(['rho_',indexKit(kit1)]) =  rho1;
    [alpha,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho_bis,mu,b,f,kit2,mod,verbose,log_path,paramSolver,ampl_models_dir)
    disp(alpha.k)
    vecAlpha = [vecAlpha,alpha.k];
end
%
%[aa,bb] = sort(vecRho1);
%vecRho1 = aa;
%vecAlpha = vecAlpha(bb);
%save("rhoprime_sc_b_1_p_10_04.mat")
%%
%load("rhoprime_sc_b_1_p_10_04.mat")
close all
fig = figure()
plot(vecRho1([1:6,11:end]), vecAlpha([1:6,11:end]), 'k--','LineWidth',0.5)
hold on
plot(vecRho1(1:6), vecAlpha(1:6), 'k', 'LineWidth',1.5)
plot(vecRho1(11:end), vecAlpha(11:end), 'k','LineWidth',1.5)
title(['$\rho_{',indexKit(kit2),'}^\prime(\rho_{',indexKit(kit1),'})$'],'Interpreter','latex','FontSize',fontSize)
xlabel(['$\rho_{',indexKit(kit1),'}$'],'interpreter','latex','FontSize',fontSize)
ylabel(['$\rho_{',indexKit(kit2),'}^\prime$'],'interpreter','latex','FontSize',fontSize)
yticks(0.38:0.02:0.45)
ylim([0.38,0.45])
xticks(0:0.2:1)
ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];
set(fig,'position',[100,100,390,330])
thisFigPath = [fig_path,'rhoprime_',indexKit(kit2),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

%% Evolution de rho_hg' en fonction de rho_sc
close all; clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit2  = {'HIV','Ng'};
kit1     = {'syphilis','Ct'};
paramSolver.varToChange = 'rho_hg';
fontSize = 18;
getParameters_avg; disp(num2str(ID_ech))
%----------------------------------%
vecAlpha=[];
paramRho_bis = paramRho;
vecRho1 = [0:0.1:1]; %pour rho_sc'(rho_hg)

%vecRhobis = 0.62;
%vecRho1 = [vecRho1,vecRhobis]
for rho1=vecRho1
    disp(rho1)
    paramRho_bis.(['rho_',indexKit(kit1)]) =  rho1;
    [alpha,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho_bis,mu,b,f,kit2,mod,verbose,log_path,paramSolver,ampl_models_dir)
    disp(alpha)
    vecAlpha = [vecAlpha,alpha.k];    
end

%
%[aa,bb] = sort(vecRho1);
%vecRho1 = aa;
%vecAlpha = vecAlpha(bb);
%save("rhoprime_hg_b_1_p_10_04.mat")
%%
close all
fig = figure()
plot(vecRho1, round(vecAlpha,3), 'k-','LineWidth',0.5)
title(['$\rho_{',indexKit(kit2),'}^\prime(\rho_{',indexKit(kit1),'})$'],'Interpreter','latex','FontSize',fontSize)
xlabel(['$\rho_{',indexKit(kit1),'}$'],'interpreter','latex','FontSize',fontSize)
ylabel(['$\rho_{',indexKit(kit2),'}^\prime$'],'interpreter','latex','FontSize',fontSize)
ylim([0,1.6*max(vecAlpha)])
xticks(0:0.2:1)
ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];
set(fig,'position',[100,100,390,330])
thisFigPath = [fig_path,'rhoprime_',indexKit(kit2),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)


%% rhohat_hscg as a function of c_hscg
close all; clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit     = {'HIV','syphilis','Ct','Ng'};
fontSize= 18;
verbose = 0; paramSolver.verbose=0;
getParameters_avg; disp(['ID=',num2str(ID_ech)])
%----------------------------------%
paramSolver.tolP0 = 0.5e-4;
%[alphas,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir);
paramsC         = paramSolver;
paramsC.sup     = 1;
paramsC.inf     = -5;
paramsC.tolC    = 1e-3;
paramsC.iterMax = 60;
afficherOutput  = 0;
paramSolver.verbose       = 0;
paramSolver.sup_bnd_alpha = 20;
paramSolver.timeSolver    = 20;
paramSolver.tolP0         = paramSolver.tolP0;

alphas.elim_h=0;
TABalphas = struct2table(alphas);
[cnn,tabP,msg] = find_Cnn_kit_2(paramTab,mu,b,paramRho,f,kit,paramsC,TABalphas,0,log_path,paramSolver,[],ampl_models_dir)

%%
vecC   = [linspace(-0.8,cnn.kit-1e-2,3),...
    linspace(cnn.kit-1e-3,cnn.kit+1e-3,3),...
    linspace(cnn.kit+1e-2,cnn.Ng-1e-3,5),...
    linspace(cnn.Ng+1e-3,cnn.Ct-1e-3,5),...
    linspace(cnn.Ct+1e-3,1.2,10)];
n = length(vecC);
Ph  = zeros(1,n);
Ps  = zeros(1,n);
Pc  = zeros(1,n);
Pg  = zeros(1,n);
Pk_r = zeros(1,n);
newmod = 1; %to iniitialize ampl
ampl_c = 0; vecRhohat = []; i=1;
for c=vecC
    mySeed = 1;
    [rhohat,Cval,ES,msg_tot,ampl_c] = findRhohat_kit_v9(paramTab,mu,b,f,paramRho,kit,c,TABalphas,verbose,mySeed,log_path,paramSolver,ampl_c,newmod,ampl_models_dir);
    vecRhohat = [vecRhohat,rhohat.tot];
    Ph(i) = P_kit_v3(ES.tot,{'HIV'},f);
    Pc(i) = P_kit_v3(ES.tot,{'Ct'},f);
    Ps(i) = P_kit_v3(ES.tot,{'syphilis'},f);
    Pg(i) = P_kit_v3(ES.tot,{'Ng'},f);
    [Pk_r(i),~,~] = P_kit_v3(ES.tot,kit,f);
    newmod = 0;
    i=i+1;
end
ampl_c.close();
%save("rhohat_hscg_chscg04.mat")

%%
%load("rhohat_hscg_chscg.mat")
close all
idx = max(find(vecRhohat==alphas.h_h));
fig = figure()
%set(groot,'defaultAxesTickLabelInterpreter','latex');  

plot(vecC, vecRhohat, 'k--','LineWidth',0.5)
hold on
yline(alphas.hscg_hscg,':','LineWidth',1.,'Color',[0,0,0])
yline(alphas.hsg_hsg,':','LineWidth',1.,'Color',[0,0,0])
yline(alphas.hs_hs,':','LineWidth',1.,'Color',[0,0,0])
yline(alphas.h_h,':','LineWidth',1.,'Color',[0,0,0])
%hold on
plot(vecC(1:idx),vecRhohat(1:idx),'k-','LineWidth',3)
plot(vecC(idx+1:end),vecRhohat(idx+1:end),'k-','LineWidth',3)
title(['$\hat\rho_{',indexKit(kit),'}(c_{',indexKit(kit),'})$'],'Interpreter','latex','FontSize',fontSize)
xlabel(['$c_{',indexKit(kit),'}$'],'interpreter','latex','FontSize',fontSize)
ylabel(['$\hat\rho_{',indexKit(kit),'}$'],'interpreter','latex','FontSize',fontSize)

yticks([0,0.92*alphas.hscg_hscg,1.3*alphas.hsg_hsg,0.25,alphas.hs_hs,0.5,alphas.h_h])
yticklabels({'0','$\rho_{hscg,c}^\prime$','$\rho_{hscg,g}^\prime$','0.25','$\rho_{hscg,s}^\prime$','0.5','$\rho_{hscg,h}^\prime$'})
xlim([vecC(1),vecC(end)])
ylim([0,1.2*max(vecRhohat)])
plot(min(vecC)+[0,0.03],[alphas.hscg_hscg,alphas.hscg_hscg],'k-','LineWidth',0.1)
plot(min(vecC)+[0,0.03],[alphas.hsg_hsg,alphas.hsg_hsg],'k-','LineWidth',0.1)
plot(min(vecC)+[0,0.03],[alphas.hs_hs,alphas.hs_hs],'k-','LineWidth',0.1)
plot(min(vecC)+[0,0.03],[alphas.h_h,alphas.h_h],'k-','LineWidth',0.1)
plot(min(vecC)+[0,0.03],[0.5,0.5],'k-','LineWidth',0.1)
plot(min(vecC)+[0,0.03],[0.25,0.25],'k-','LineWidth',0.1)

ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];

set(fig,'position',[100,100,400,340])

thisFigPath = [fig_path,'rhohat_',indexKit(kit),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'2.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)
%%
paramsC = paramSolver;
paramsC.sup=1;
paramsC.inf=-5;
paramsC.tolC=1e-3;
paramsC.iterMax=30;
afficherOutput=0;

createParamRho;
[cnn,tabP,msg] = find_Cnn_kit_2(paramTab,mu,b,paramRho,f,kit,paramsC,TABalphas,afficherOutput,log_path,paramSolver,[],ampl_models_dir);
%%
close all
fig = figure();
plot(vecC,Pk_r,'k--','LineWidth',0.5,'HandleVisibility','off')
hold on
plot(vecC(1:idx),Pk_r(1:idx),'k-','LineWidth',2,'HandleVisibility','off')
plot(vecC(idx+1:end),Pk_r(idx+1:end),'k-','LineWidth',2)

xline(cnn.HIV,':','LineWidth',1,'Color',[0,0,0],'HandleVisibility','off')
xline(cnn.syphilis,':','LineWidth',1,'Color',[0,0,0],'HandleVisibility','off')
xline(cnn.Ct,':','LineWidth',1,'Color',[0,0,0],'HandleVisibility','off')
xline(cnn.Ng,':','LineWidth',1,'Color',[0,0,0],'HandleVisibility','off')

plot(vecC,Ph,'--','Color',rougeHIV,'LineWidth',0.5,'HandleVisibility','off')
plot(vecC(1:idx),Ph(1:idx),'Color',rougeHIV,'LineWidth',2,'HandleVisibility','off')
plot(vecC(idx+1:end),Ph(idx+1:end),'Color',rougeHIV,'LineWidth',2)

plot(vecC,Ps,'--','Color',jauneS,'LineWidth',0.5,'HandleVisibility','off')
plot(vecC(1:idx),Ps(1:idx),'Color',jauneS,'LineWidth',2,'HandleVisibility','off')
plot(vecC(idx+1:end),Ps(idx+1:end),'Color',jauneS,'LineWidth',2)

%plot(vecC,Pc,'--','Color',bleuCt,'LineWidth',0.5,'HandleVisibility','off')
%plot(vecC(1:idx),Pc(1:idx),'Color',bleuCt,'LineWidth',2,'HandleVisibility','off')
plot(vecC(1:end),Pc(1:end),'Color',bleuCt,'LineWidth',2)

%plot(vecC,Pg,'--','Color',vertNg,'LineWidth',0.5,'HandleVisibility','off')
%plot(vecC(1:idx),Pg(1:idx),'Color',vertNg,'LineWidth',2,'HandleVisibility','off')
plot(vecC(1:end),Pg(1:end),'Color',vertNg,'LineWidth',2)

xlim([vecC(1),vecC(end)])
lgd = legend(['$\hat\Pi_{',indexKit(kit),'}$'],'$\hat\Pi_h$','$\hat\Pi_s$','$\hat\Pi_c$','$\hat\Pi_g$',...
    'Interpreter','latex','FontSize',fontSize,'Location','Northwest');
xlabel(['$c_{',indexKit(kit),'}$'],'FontSize',fontSize,'Interpreter','latex')
ylabel('Prevalence','Interpreter','latex','FontSize',fontSize)
%title('','Interpreter','latex','FontSize',fontSize)
legend('boxoff') 
x1=18 ; x2=18;
lgd.ItemTokenSize = [x1,x2];
ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];
set(fig,'position',[100,100,400,340])
thisFigPath = [fig_path,'prevalence_',indexKit(kit),'_c_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

%% Utility concave par morceaux
close all; clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hcg_hcg';
kit     = {'HIV','Ct','Ng'};
paramSolver.varToChange = 'rho_hcg';
vecRho  = [0:0.01:0.11, 0.105:0.005,0.14, 0.145:0.01:0.55, 0.56:0.005:0.64, 0.65:0.01:1]; 
fontSize = 18;
getParameters_avg; disp(num2str(ID_ech))
c=0%-0.01;
%----------------------------------%

Ph = [];Pc = [];Pg = [];Ps = [];Ph_p = [];Pk_r = zeros(1,length(vecRho));
ampl = 0;
i=1; paramRho_bis = paramRho;

for rho = vecRho
    paramRho_bis.(paramSolver.varToChange) = rho;
    [P,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl); 
    Ph(i) = P.h;
    Pc(i) = P.c;
    Ps(i) = P.s;
    Pg(i) = P.g;
    [Pk_r(i),~,~] = P_kit_v3(ES,kit,f); %prevalence of asymptomatic and undiagnosed
    %[Pk_r(i),~,~] = P_kit_v3(ES,kit,f); %prevalence of asymptomatic and undiagnosed
    disp([num2str(i),'/',num2str(length(vecRho))])
    i=i+1;
end
ampl.close();

%%
[alphas,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir);
%%
c=-0.01;
Uk = vecRho.*(Pk_r-c);
fig = figure();
plot(vecRho,Uk,'k-','LineWidth',3)

xline(alphas.hcg_hcg,':','LineWidth',1,'Color',[0,0,0])
xline(alphas.hg_hg,':','LineWidth',1,'Color',[0,0,0])
xline(alphas.h_h,':','LineWidth',1,'Color',[0,0,0])

title(['${\rm U}_{',indexKit(kit),'}(\rho_{',indexKit(kit),'})$'],'Interpreter','latex','FontSize',fontSize)
xlabel(['$\rho_{',indexKit(kit),'}$'],'FontSize',fontSize,'Interpreter','latex')
ylabel('Utility','Interpreter','latex','FontSize',fontSize)
%xlim([0,0.7])
yticks([0,0.004,0.008,0.012])
ylim([0,max(Uk)*1.2])
ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];
set(fig,'position',[100,100,400,350])
thisFigPath = [fig_path,'utility_',indexKit(kit),'_max_loc_1_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'c_',num2str(c),'_',num2str(ID_ech),'.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

%see animationU.m
%% Reproduction number and prevalences of the single-disease models as a function of rho
close all; clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm
%----------------------------------%
sictp=0; seiiis=1; seiis=0;
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
fontSize = 18;
getParameters_avg; disp(num2str(ID_ech))
addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\INITIALISATION\parametersSets')
%----------------------------------%
vecRho  = 0:0.01:1.5; n=length(vecRho);
Rp = zeros(1,n);
Ptot = zeros(1,n); IS = Ptot;
Pun = zeros(1,n);
mu = paramTab{4}.mu; b = 10;
for i=1:n
    if sictp==1
        betaI = paramTab{3}.betaI; betaC = paramTab{3}.betaC; theta0 = paramTab{3}.theta0; 
        gamma0 = paramTab{3}.gamma0; sigma = paramTab{3}.sigma; zeta = paramTab{3}.zeta;
        eta = paramTab{3}.eta;
        [Rp(i),~,alpha,Ptot(i),Pun(i)] = Rp_SICTP(betaI,betaC,theta0,gamma0,sigma,zeta,eta,pHIV,mu,b,vecRho(i));
        folderName = 'SICTP';
    else 
        if seiiis==1
            beta = paramTab{4}.beta; sigma=paramTab{4}.sigma; tau=paramTab{4}.tau;
            nu=paramTab{4}.nu; gamma10=paramTab{4}.gamma10; theta=paramTab{4}.theta;
            gamma30=paramTab{4}.gamma30; 
            [Rp(i),~,alpha] = Rp_SEIIIS_v4(beta,sigma,tau,nu,gamma10,theta,gamma30,mu,b,vecRho(i));
            Ptot(i) = 1-1./Rp(i);
            Pun(i)  = Ptot(i);
            folderName = 'SEIIIS';
        else
            if seiis==1
                beta = paramTab{1}.beta; nu=paramTab{1}.nu; epsilon = paramTab{1}.eps;
                sigma=paramTab{1}.sigma; gamma0=paramTab{1}.gamma; 
                [Rp(i),~,alpha,Pun(i)] = Rp_SEIIS_v4(beta,nu,epsilon,sigma,gamma0,mu,b,vecRho(i));
                Ptot(i) = 1-1./Rp(i);
                folderName = 'SEIIS';
            end
        end
    end
end

myData = [vecRho;Rp;Ptot;Pun]';
myData = round(myData,4)
%writematrix(myData, ['C:\Users\Moi\Documents\IPLESP\These\Articles\Drafts\Article2_v5\supplementary_material\figures\',folderName, '/Rp_Ptot_Pun.dat']) 

%%
close all
fig1 = figure()
plot(vecRho,Rp,'k','LineWidth',3)
ylim([0 max(Rp)*1.3])
yline(1,':')
xlabel('$\rho$','Interpreter','latex','FontSize',fontSize)
ylabel('$\mathtt R(\rho)$','Interpreter','latex','FontSize',fontSize)
set(fig1,'position',[100,100,300,250])
xlim([vecRho(1),vecRho(end)])

fig2 = figure()
plot(vecRho,Ptot,'k','LineWidth',3)
ylim([0 max(Ptot)*1.3])
yline(0,':')
xlabel('$\rho$','Interpreter','latex','FontSize',fontSize)
ylabel('$\Pi^{\rm all}(\rho)$','Interpreter','latex','FontSize',fontSize)
set(fig2,'position',[100+400,100,300,250])
xlim([vecRho(1),vecRho(end)])

fig3=figure()
plot(vecRho,Pun,'k','LineWidth',3)
ylim([0 max(Pun)*1.3])
yline(0,':')
xlabel('$\rho$','Interpreter','latex','FontSize',fontSize)
ylabel('$\Pi(\rho)$','Interpreter','latex','FontSize',fontSize)
set(fig3,'position',[100+800,100,300,250])
xlim([vecRho(1),vecRho(end)])




%% Evolution de la prévalence percue
%close all; 
clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm parula_w
%----------------------------------%
vecC = [-0.5,-0.25,0,0.25,0.5];
vecP = 0:0.1:0.8;
vecRho  = 0:0.01:0.7;%[0:0.01:0.38,0.385:0.001:0.40, 0.41:0.01:0.66,0.661:0.001:0.68,1];


roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;%not useful

mod     = 'hscg_hscg';
kit     = {'HIV','syphilis'};
paramSolver.varToChange = ['rho_',indexKit(kit)];
fontSize = 18;
%----------------------------------%
c=0;
Ptot = zeros(length(vecRho), length(vecP));
k=1;
for pHIV=vecP
    getParameters_avg; disp(num2str(ID_ech))
    Pk_r = zeros(1,length(vecRho));
    ampl = 0;
    i=1;
    paramRho_bis = paramRho;
    for rho = vecRho
        paramRho_bis.(paramSolver.varToChange) = rho;
        [P,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl); 
        Ph(i) = P.h;
        Pc(i) = P.c;
        Ps(i) = P.s;
        Pg(i) = P.g;
        [Pk_r(i),Pk_p(i),Ph_p(i)] = P_kit_v3(ES,kit,f);
        disp([num2str(i),'/',num2str(length(vecRho))])
        i=i+1;
    end
    ampl.close();

    % (prevalence)
    figure(1);
    U = vecRho.*(Pk_r-c);
    plot(vecRho,U,'LineWidth',2)
    hold on
    Ptot(:,k) = Pk_r; k=k+1;
end

save('Ptot.mat','Ptot');

xlabel(['$\rho_{',indexKit(kit),'}$'],'FontSize',fontSize,'Interpreter','latex')
ylabel('Prevalence','Interpreter','latex','FontSize',fontSize)
xlim([0,max(vecRho)])
legend('boxoff')
lgd.NumColumns = 2;
x1=18 ; x2=18;
lgd.ItemTokenSize = [x1,x2];

ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];

set(fig,'position',[100,100,390,330])
thisFigPath = [fig_path,'prevalence_',indexKit(kit),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

%%
% p=0.4
% Pi = 
% for p=vecP
%     for c=vecC
%         U=vecRho.*(P
%         figure(2)
%         plot(vecRho, U)
%     end
% end


%% HIV syphilis versus HIV sans syphilis

close all;
clearvars -except log_path fig_path ampl_models_dir rougeHIV jauneS bleuCt vertNg WGB_log_cm
%----------------------------------%
pHIV    = 0.;%1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hscg_hscg';
kit     = {'HIV'};
paramSolver.varToChange = 'rho_h';
vecRho  =  [0,0.05,0.1,0.11,0.12,0.13,0.14,0.2,0.3,0.35,0.39,0.4,0.5,0.6,0.65,0.7,1];
vecRho = sort(unique([vecRho,linspace(0,1,100)]));
fontSize = 18;
getParameters_avg; disp(num2str(ID_ech))
%----------------------------------%
%[alphas,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit,mod,0,log_path,paramSolver,ampl_models_dir);

paramTab{3}.p = pHIV;
Ph = [];Pc = [];Pg = [];Ps = [];Pk_r = zeros(1,length(vecRho));
ampl = 0;
i=1; paramRho_bis = paramRho;

for rho = vecRho
    disp([num2str(i),'/',num2str(length(vecRho))])
    disp(rho)
    paramRho_bis.(paramSolver.varToChange) = rho;
    [P,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl); 
    Ph(i) = P.h;
    disp(Ph(i));
    Pc(i) = P.c;
    Ps(i) = P.s;
    Pg(i) = P.g;
    [Pk_r(i),~,~] = P_kit_v3(ES,kit,f); %prevalence of asymptomatic and undiagnosed
    disp(Pk_r(i))
    i=i+1;
end
ampl.close();

fig = figure();
plot(vecRho,Pk_r,'k-','LineWidth',2)
hold on
plot(vecRho,Ph,'Color',rougeHIV,'LineWidth',2)
plot(vecRho,Ps,'Color',jauneS,'LineWidth',2)
plot(vecRho,Pc,'Color',bleuCt,'LineWidth',2)
plot(vecRho,Pg,'Color',vertNg,'LineWidth',2)
lgd = legend(['$\Pi_{',indexKit(kit),'}$'],'$\Pi_h$','$\Pi_s$','$\Pi_c$','$\Pi_g$','Interpreter','latex','FontSize',fontSize);
xlabel(['$\rho_{',indexKit(kit),'}$'],'FontSize',fontSize,'Interpreter','latex')
ylabel('Prevalence','Interpreter','latex','FontSize',fontSize)
xlim([0,0.7])
legend('boxoff') 
x1=18 ; x2=18;
lgd.ItemTokenSize = [x1,x2];

ax=gca;
ax.XAxis.FontSize = 0.75*fontSize;
ax.YAxis.FontSize = 0.75*fontSize;
ax.XLabel.FontSize = fontSize;
ax.YLabel.FontSize = fontSize;
ax.Title.FontSize = fontSize;
ax.YAxis.TickLength = [0,0];

set(fig,'position',[100,100,400,340])
thisFigPath = [fig_path,'prevalence_',indexKit(kit),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str(ID_ech),'.pdf'];
%exportgraphics(fig,thisFigPath,'Resolution',300)

temp = [vecRho;Pk_r;Ph;Ps;Pc;Pg]'
temp=round(temp,5)





