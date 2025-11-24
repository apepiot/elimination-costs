%% Utility concave par morceaux
close all;
addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
addpath('MAIN')
log_path = 'C:/Users/Moi/Desktop/Temporaire/tests/';
fig_path = 'C:\Users\Moi\Documents\IPLESP\These\Articles\Article 2\figures\';
ampl_models_dir = [pwd,'\MAIN\AMPL_models\'];
setupOnce;
close all; clearvars -except log_path fig_path ampl_models_dir 
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
getParameters; disp(num2str(ID_ech))
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

%
[alphas,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir);

%%
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'h_h';
kit     = {'HIV'};
paramSolver.varToChange = 'rho_h';
vecRho  = [0:0.01:0.11, 0.105:0.005,0.14, 0.145:0.01:0.55, 0.56:0.005:0.64, 0.65:0.01:1]; 
fontSize = 18;
getParameters; disp(num2str(ID_ech))
%----------------------------------%

Ph = [];Pc = [];Pg = [];Ps = [];Ph_p = [];Pk_r_hiv = zeros(1,length(vecRho));
ampl = 0;
i=1; paramRho_bis = paramRho;

for rho = vecRho
    paramRho_bis.(paramSolver.varToChange) = rho;
    [P,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl); 
    Ph(i) = P.h;
    Pc(i) = P.c;
    Ps(i) = P.s;
    Pg(i) = P.g;
    [Pk_r_hiv(i),~,~] = P_kit_v3(ES,kit,f); %prevalence of asymptomatic and undiagnosed
    %[Pk_r(i),~,~] = P_kit_v3(ES,kit,f); %prevalence of asymptomatic and undiagnosed
    disp([num2str(i),'/',num2str(length(vecRho))])
    i=i+1;
end
ampl.close();

%
[alphas_h,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir);

%%
close all
Umax = vecRho.*(Pk_r+0.02);
%h = animatedline('Marker','o');
fig=figure();
set(fig,'position',[100,100,400,350])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
%gifFile = 'myAnimation.gif'; %%%
%exportgraphics(fig, gifFile);
set(0,'defaultfigurecolor',[1 1 1])
pause(2)
for c=[0.3:-0.01:0.029,0.03:-0.001:-0.03]
    clf;
    Uk = vecRho.*(Pk_r-c);
    Uk_h = vecRho.*(Pk_r_hiv-c);
    plot(vecRho,Uk,'k-','LineWidth',2)
    hold on
    plot(vecRho,Uk_h,'--','LineWidth',2)
    xline(alphas.hcg_hcg,':','LineWidth',0.5,'Color',[0.5,0.5,0.5])
    xline(alphas.hg_hg,':','LineWidth',0.5,'Color',[0.5,0.5,0.5])
    xline(alphas.h_h,':','LineWidth',0.5,'Color',[0.5,0.5,0.5])

    title(['${\rm U}_{',indexKit(kit),'}(\rho_{',indexKit(kit),'})$ with $c_{hcg}=',num2str(c),'$'],'Interpreter','latex','FontSize',0.9*fontSize)
    xlabel(['$\rho_{',indexKit(kit),'}$'],'FontSize',fontSize,'Interpreter','latex')
    ylabel('Utility','Interpreter','latex','FontSize',fontSize)
    ylim([0,max(Umax)])
    xlim([0,max(vecRho)])
    
    %text((min(vecRho)+max(vecRho))/4, max(Umax)*0.8,['$c=',num2str(c),'$'],'Interpreter','latex','FontSize',16)
    
    [maxU,imax] = max(Uk(vecRho<=alphas.h_h));
    rhohat = vecRho(imax);
    xticks(sort(unique([0,rhohat,1])))
    
    if abs(rhohat-alphas.h_h)<1e-2
        xticklabels({'0','$\hat\rho_{hcg}=\rho_{hcg,h}^\prime$','1'})
        xtickangle(0)
    else
        xticklabels({'0','$\hat\rho_{hcg}$','1'})
    end
    xline(rhohat,'r')
    yline(maxU,'r')
    ax=gca;
    ax.XAxis.FontSize = 0.75*fontSize;
    ax.YAxis.FontSize = 0.75*fontSize;
    ax.XLabel.FontSize = fontSize;
    ax.YLabel.FontSize = fontSize;
    ax.Title.FontSize = fontSize;
    %ax.YAxis.TickLength = [0,0];
    %drawnow limitrate
    pause(0.1)
    %drawnow
    %frame = getframe(fig);
   	%exportgraphics(obj, gifFile, Append=true);
end
%close;



%% Utility concave par morceaux
close all;
addpath 'C:/Program Files (x86)/Artelys/AMPL 13.1.20220703/amplapi/examples/matlab';
addpath('MAIN')
log_path = 'C:/Users/Moi/Desktop/Temporaire/tests/';
fig_path = 'C:\Users\Moi\Documents\IPLESP\These\Articles\Article 2\figures\';
ampl_models_dir = [pwd,'\MAIN\AMPL_models\'];
setupOnce;
close all; clearvars -except log_path fig_path ampl_models_dir 
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'hsc_hsc';
kit     = {'HIV','syphilis','Ct'};
paramSolver.varToChange = 'rho_hsc';
vecRho  = [0:0.01:0.11, 0.105:0.005,0.14, 0.145:0.01:0.55, 0.56:0.005:0.64, 0.65:0.01:1]; 
fontSize = 18;
getParameters; disp(num2str(ID_ech))
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

%
[alphas,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir);

%%
%----------------------------------%
pHIV    = 0.1;
roundNo = 1;
f       = 1;
paramNo = 4;
nbEch   = 8;
mod     = 'h_h';
kit     = {'HIV'};
paramSolver.varToChange = 'rho_h';
vecRho  = [0:0.01:0.11, 0.105:0.005,0.14, 0.145:0.01:0.55, 0.56:0.005:0.64, 0.65:0.01:1]; 
fontSize = 18;
getParameters; disp(num2str(ID_ech))
%----------------------------------%

Ph = [];Pc = [];Pg = [];Ps = [];Ph_p = [];Pk_r_hiv = zeros(1,length(vecRho));
ampl = 0;
i=1; paramRho_bis = paramRho;

for rho = vecRho
    paramRho_bis.(paramSolver.varToChange) = rho;
    [P,ES,~,ampl] = P_mod_v8(paramTab,paramRho_bis,b,mu,mod,verbose,paramSolver,log_path,ampl_models_dir,ampl); 
    Ph(i) = P.h;
    Pc(i) = P.c;
    Ps(i) = P.s;
    Pg(i) = P.g;
    [Pk_r_hiv(i),~,~] = P_kit_v3(ES,kit,f); %prevalence of asymptomatic and undiagnosed
    %[Pk_r(i),~,~] = P_kit_v3(ES,kit,f); %prevalence of asymptomatic and undiagnosed
    disp([num2str(i),'/',num2str(length(vecRho))])
    i=i+1;
end
ampl.close();

%
[alphas_h,Palpha,ES,msg,elim_i] = findAlpha_kit(paramTab,paramRho,mu,b,f,kit,mod,verbose,log_path,paramSolver,ampl_models_dir);

%%
close all
Umax = vecRho.*(Pk_r+0.5);
%h = animatedline('Marker','o');
fig=figure();
set(fig,'position',[100,100,800,350])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
%gifFile = 'myAnimation.gif'; %%%
%exportgraphics(fig, gifFile);
set(0,'defaultfigurecolor',[1 1 1])
rougeHIV = [215, 0, 0]/255;
pause(2)
vecC = [0.5:-0.01:-0.4]; vecRhohat=[]; vecRhohat_h=[];
k=1;
for c=vecC
    clf;
    subplot(1,2,1)
    Uk   = vecRho.*(Pk_r-c);
    Uk_h = vecRho.*(Pk_r_hiv-c);
    plot(vecRho,Uk,'k-','LineWidth',2)
    hold on
    plot(vecRho,Uk_h,'--','Color',rougeHIV,'LineWidth',2)
    xline(alphas.hsc_hsc,':','LineWidth',0.5,'Color',[0.5,0.5,0.5],'HandleVisibility','off');
    xline(alphas.hs_hs,':','LineWidth',0.5,'Color',[0.5,0.5,0.5],'HandleVisibility','off');
    xline(alphas.h_h,':','LineWidth',0.5,'Color',[0.5,0.5,0.5],'HandleVisibility','off');

    title(['$c=',num2str(c),'$'],'Interpreter','latex','FontSize',0.9*fontSize)
    xlabel(['$\rho$'],'FontSize',fontSize,'Interpreter','latex')
    ylabel('Utility','Interpreter','latex','FontSize',fontSize)
    ylim([0,max(Umax)])
    xlim([0,max(vecRho)])
    legend('${\rm U}_{hsc}(\rho_{hsc})$', '${\rm U}_h(\rho_h)$','Interpreter','latex',...
        'FontSize',0.8*fontSize,'Location','NorthWest','Color', [1,1,1])
    %legend('boxoff')
    
    %text((min(vecRho)+max(vecRho))/4, max(Umax)*0.8,['$c=',num2str(c),'$'],'Interpreter','latex','FontSize',16)
    
    [maxU,imax] = max(Uk(vecRho<=alphas.h_h));
    rhohat = vecRho(imax); vecRhohat = [vecRhohat,rhohat];
    xline(rhohat,'k','HandleVisibility','off');
    yline(maxU,'k','HandleVisibility','off');   
    
    [maxU_h,imax] = max(Uk_h(vecRho<=alphas_h.h_h));
    rhohat_h = vecRho(imax); vecRhohat_h = [vecRhohat_h,rhohat_h];
    
    xline(rhohat_h,'--','Color',rougeHIV,'HandleVisibility','off');
    yline(maxU_h,'--','Color',rougeHIV,'HandleVisibility','off');
    
    [rhohat_ticks,idxs] = sort([0,rhohat_h,rhohat,alphas.h_h,1]);
    rhohat_labels = {'0','$\hat\rho_{h}$','$\hat\rho_{hsc}$','$\rho_{h}^\prime$','1'};
    rhohat_labels_sorted = rhohat_labels(idxs);
    
    [rhohat_ticks,idxs_2]= unique(rhohat_ticks);
    rhohat_labels_sun    = rhohat_labels_sorted(idxs_2);
    xticks(rhohat_ticks)
    xticklabels(rhohat_labels_sun);
    

    ax=gca;
    ax.XAxis.FontSize = 0.75*fontSize;
    ax.YAxis.FontSize = 0.75*fontSize;
    ax.XLabel.FontSize = fontSize;
    ax.YLabel.FontSize = fontSize;
    ax.Title.FontSize = fontSize;
    
    subplot(1,2,2)
    plot(vecC(1:k),vecRhohat(1:k),'k-','LineWidth',2)
    hold on
    plot(vecC(1:k),vecRhohat_h(1:k),'r--','LineWidth',2)
    xlim([vecC(end),vecC(1)])
    ylim([0,alphas.h_h*1.2])
    ylabel('$\hat\rho(c)$','Interpreter','latex')
    xlabel('$c$','Interpreter','latex')
    legend('${\hat\rho}_{hsc}(c_{hsc})$', '${\hat\rho}_{h}(c_{h})$','Interpreter','latex',...
        'FontSize',0.8*fontSize,'Location','NorthEast','Color', [1,1,1])
    ax=gca;
    ax.XAxis.FontSize = 0.75*fontSize;
    ax.YAxis.FontSize = 0.75*fontSize;
    ax.XLabel.FontSize = fontSize;
    ax.YLabel.FontSize = fontSize;
    ax.Title.FontSize = fontSize;
%     if length(unique(round(rhohat_ticks,3)))==length(rhohat_ticks)
%         if rhohat_h<rhohat
%             xticklabels({'0','$\hat\rho_{h}$','$\hat\rho_{hsc}$','$\rho_{h}^\prime$','1'})
%         else
%             xticklabels({'0','$\hat\rho_{hsc}$','$\hat\rho_{h}$','$\rho_{hsc,h}^\prime$','1'})
%         end
%     else
%         
%     end
% 
%     if abs(rhohat-alphas.h_h)<1e-2 & abs()
%         xticklabels({'0','$\hat\rho_{hsc}$','$\hat\rho_{hsc}=\rho_{hsc,h}^\prime$','1'})
%         xtickangle(0)
%     else
%         xticklabels({'0','$\hat\rho_{hsc}$','1'})
%    end
    
    
    
    %ax.YAxis.TickLength = [0,0];
    %drawnow limitrate
    pause(0.001)
    %drawnow
    %frame = getframe(fig);
   	%exportgraphics(obj, gifFile, Append=true);
    
    
   
    
    k=k+1;
    
end
%close;