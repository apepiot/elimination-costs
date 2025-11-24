%% plot
% functions plot
close all
FontSizeGlobal = 14;
bleu = [0,0.4470,0.7410];jaune=[0.9290,0.6940,0.1250];rouge=[0.900,0.3250,0.0980];violet=[0.4940, 0.1840, 0.5560];vert=[0.4660, 0.6740, 0.1880];

modelTypes = repelem({'SEIIS', 'SICR', 'SEIIIS'}, [nSEIIS nSICR nSEIIIS]);
infections = repelem({'STI', 'HIV', 'syphilis'}, [nSEIIS nSICR nSEIIIS]);
infectionsMini = repelem({'STI', 'h', 's.'}, [nSEIIS nSICR nSEIIIS]);
if nSEIIS==2
    infections{1} = 'Ct'; %we replace STI by Ct
    infectionsMini{1} = 'c';
    infections{2} = 'Ng';
    infectionsMini{2} = 'g';
elseif nSEIIS==1
    infections{1} = paramTab{1}.disease;
    infectionsMini{1} = paramTab{1}.disease;
end

set(groot,'defaultAxesTickLabelInterpreter','latex');
figure()
maxalpha = max(vecAlphas);
limy=1.4*maxalpha;
ax = gca;

if N==1
    plot_paper1_procedure3(tab.rhohat.rhohat, cs.cs1dis, cs.c0, paramTab{1}.alpha, vecC, tab.rhohat.disease)
    yticks([ 0 vecAlphas(1)])
    yticklabels({'0',['$\rho_{',infectionsMini{1},'}^\prime$']})
elseif N>=2
    hold on
    %p_rho1  = plot(vecC,tab.one(1).rhohat,':','Linewidth',2,'Color','k','DisplayName',['$\hat\rho_{',infectionsMini{1},'}$']); %bleu/SEIIS 1
    %test
    p_rho1  = plot(vecC,tab.one(1).rhohat,':','Linewidth',2,'Color','k','DisplayName',['\fontsize{22}{0}\selectfont$\hat\rho$\fontsize{14}{0}\selectfont$_{',paramTab{1}.disease,'}$']); %bleu/SEIIS 1    
    %p_rho2  = plot(vecC,tab.one(2).rhohat,'--','Linewidth',2,'Color',rouge,'DisplayName',['$\hat\rho_{',infectionsMini{2},'}$']); %jaune/SEIIS 1
    p_rho2  = plot(vecC,tab.one(2).rhohat,'--','Linewidth',2,'Color',rouge,'DisplayName',['\fontsize{22}{0}\selectfont$\hat\rho$\fontsize{14}{0}\selectfont$_{',paramTab{2}.disease,'}$']); %jaune/SEIIS 1
    
    if N==2
        %%
        %p_rhohat = plot(vecC(vecC<=cs.cs1dis),tab.rhohat.rhohat(vecC<=cs.cs1dis),'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$'); %bleu/SEIIS 1
        p_rhohat = plot(vecC(vecC<=cs.cs1dis),tab.rhohat.rhohat(vecC<=cs.cs1dis),'-','Linewidth',3,'Color','k',...
            'DisplayName','\fontsize{22}{0}\selectfont$\hat\rho$'); %bleu/SEIIS 1
        hold on
        %plot(vecC(vecC>cs.cs1dis),tab.rhohat.rhohat(vecC>cs.cs1dis),'-','Linewidth',3,'Color','k',...
        %    'DisplayName','$\hat\rho$');
        plot(vecC(vecC>cs.cs1dis),tab.rhohat.rhohat(vecC>cs.cs1dis),'-','Linewidth',3,'Color','k',...
            'DisplayName','\fontsize{22}{0}\selectfont$\hat\rho$');
        
        %p_rhohat = plot(vecC,tab.rhohat,'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$'); %bleu/SEIIS 1
        lgd = legend([p_rhohat,p_rho1,p_rho2],'Interpreter','latex');
        %yticks([ 0  vecAlphas(cs.order(1)) vecAlphas(cs.order(2)) ])
        %yticklabels({'0',['$\rho_{',infectionsMini{cs.order(1)},'}^\prime$'],...
        %    ['$\rho_{',infectionsMini{cs.order(2)},'}^\prime$']})
        ylim([0,limy])
        if tabcn.one(1)<tabcn.one(2) %%%a revoir%%%
            %xticks([-1 tabcn.one(1) tabcn.one(2) 0 0.5 1])
            %xticklabels({'-1',['$c_{',infectionsMini{1},'}^\prime$'],...
            %['$c_{',infectionsMini{2},'}^\prime$'],' ','0.5','1'})
        else
            %xticks([-1 tabcn.one(2) tabcn.one(1) 0 0.5 1])
            %xticklabels({'-1',['$c_{',infectionsMini{2},'}^\prime$'],...
            %['$c_{',infectionsMini{1},'}^\prime$'],' ','0.5','1'})
        end
    end
    if (N>=3)
        p_rho3   = plot(vecC,tab.one(3).rhohat,':','Linewidth',2,'Color',bleu,'DisplayName',['$\hat\rho_{',infectionsMini{3},'}$']);
        if N==3
            %p_rhohat = plot(vecC,tab.rhohat,'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$');
            p_rhohat = plot(vecC(vecC<=cs.cs2dis),tab.rhohat.rhohat(vecC<=cs.cs2dis),'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$'); %bleu/SEIIS 1
            plot(vecC(vecC<=cs.cs1dis & vecC>cs.cs2dis),tab.rhohat.rhohat(vecC<=cs.cs1dis & vecC>cs.cs2dis),'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$'); %bleu/SEIIS 1
            hold on
            plot(vecC(vecC>cs.cs1dis),tab.rhohat.rhohat(vecC>cs.cs1dis),'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$'); %bleu/SEIIS 1
            
            legend([p_rhohat,p_rho1,p_rho2,p_rho3],'Interpreter','latex')
            yticks([ 0 vecAlphas(cs.order(1)) vecAlphas(cs.order(2)) vecAlphas(cs.order(3))])
            yticklabels({' ',['$\rho_{',infections{cs.order(1)},'}^\prime$'],...
                ['$\rho_{',infectionsMini{cs.order(2)},'}^\prime$'],'0.5',...
                ['$\rho_{',infections{cs.order(3)},'}^\prime$']})
        end
        if (N>=4)
            p_rho4   = plot(vecC,tab.one(4).rhohat,':','Linewidth',2,'Color',vert,'DisplayName',['$\hat\rho_{',infectionsMini{4},'}$']); %vert/ SEIIIS
            
            if N==4
                p_rhohat = plot(vecC(vecC<=cs.cs4dis),tab.rhohat.rhohat(vecC<=cs.cs4dis),'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$'); %bleu/SEIIS 1
                plot(vecC(vecC<=cs.cs3dis & vecC>cs.cs4dis),tab.rhohat.rhohat(vecC<=cs.cs3dis & vecC>cs.cs4dis),'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$'); %bleu/SEIIS 1
                plot(vecC(vecC<=cs.cs2dis & vecC>cs.cs3dis),tab.rhohat.rhohat(vecC<=cs.cs2dis & vecC>cs.cs3dis),'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$'); %bleu/SEIIS 1
                plot(vecC(vecC<=cs.cs1dis & vecC>cs.cs2dis),tab.rhohat.rhohat(vecC<=cs.cs1dis & vecC>cs.cs2dis),'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$'); %bleu/SEIIS 1
                plot(vecC(vecC>cs.cs1dis),tab.rhohat.rhohat(vecC>cs.cs1dis),'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$'); %bleu/SEIIS 1
                
                %p_rhohat = plot(vecC,tab.rhohat,'-','Linewidth',3,'Color','k','DisplayName','$\hat\rho$'); %bleu/SEIIS 1
                legend([p_rhohat,p_rho1,p_rho2,p_rho3,p_rho4],'Interpreter','latex')
                yticks([ 0 vecAlphas(cs.order(1)) vecAlphas(cs.order(2)) vecAlphas(cs.order(3)) vecAlphas(cs.order(4))])
                yticklabels({' ',['$\rho_{',infectionsMini{cs.order(1)},'}^\prime$'],...
                    ['$\rho_{',infectionsMini{cs.order(2)},'}^\prime$'],...
                    ['$\rho_{',infectionsMini{cs.order(3)},'}^\prime$'], '0.5'...
                    ['$\rho_{',infectionsMini{cs.order(4)},'}^\prime$']})
            end
        end
    end
    %lgd.FontSize = 10;
end

%% areas plot
%N = nSEIIS+nSICR+nSEIIIS;
gris=[78/255, 78/255, 78/255];
%
if N==2
    %lim_xinf = cs.cs2dis-(cs.c0-cs.cs2dis)/2;
    lim_xinf = cs.cs2dis-(cs.c0-cs.cs2dis)*0.8%/2;
    lim_xsup = cs.c0+(cs.c0-cs.cs2dis)/2;
    xlim([lim_xinf, lim_xsup])
    ex = (lim_xsup-lim_xinf)/15;
    ey = limy/20;
    %%
    % zones 2 diseases
    a2=area([lim_xinf lim_xinf cs.cs2dis cs.cs2dis],[0,limy,limy,0 ],...
        'LineStyle','None','DisplayName',[infections{1},' \& ', infections{2},' eliminated']); a2(1).FaceColor = gris;a2.FaceAlpha = 0.3;
    %text((lim_xinf+cs.cs2dis)/2,limy*0.9,'III','Interpreter','latex','FontSize',FontSizeGlobal,'HorizontalAlignment', 'center')
    text((lim_xinf+cs.cs2dis)/2,limy*0.85,'\fontsize{22}{0} \selectfont III','Interpreter','latex','HorizontalAlignment', 'center')

    a1=area([cs.cs2dis cs.cs2dis cs.cs1dis cs.cs1dis],[0,limy,limy,0 ],...
        'LineStyle','None','DisplayName',['only ', num2str(infections{cs.order(1)}),' eliminated']); a1(1).FaceColor = gris;a1.FaceAlpha = 0.1;  %a0(1).EdgeColor = gris*2;
    if ((cs.cs1dis-cs.cs2dis)>(cright-cleft)/100)
        %text((cs.cs2dis+cs.cs1dis)/2,limy*0.9,'IIb','Interpreter','latex','FontSize',FontSizeGlobal,'HorizontalAlignment', 'center')
        text((cs.cs2dis+cs.cs1dis)/2,limy*0.85,'\fontsize{22}{0} \selectfont IIb','Interpreter','latex','FontSize',FontSizeGlobal,'HorizontalAlignment', 'center')
    else
        if (cs.cs1dis~=cs.cs2dis)
            %anArrow = annotation('arrow') ;
            %anArrow.Parent = gca;  % or any other existing axes or figure
            %anArrow.Position = [(cs.cs1dis+cs.cs2dis)/2+ex, limy*0.9, -ex, -ey] ;
            
            text((cs.cs2dis+cs.cs1dis)/2,limy*0.85-ey,'$\nwarrow$','Interpreter','latex','FontSize',FontSizeGlobal,'HorizontalAlignment', 'left')
            %text((cs.cs2dis+cs.cs1dis)/2 + ex,limy*0.9-ey-ey,'IIb','Interpreter','latex','FontSize',FontSizeGlobal,'HorizontalAlignment', 'left')
            text((cs.cs2dis+cs.cs1dis)/2 + ex,limy*0.85-ey-ey,'\fontsize{22}{0} \selectfont IIb','Interpreter','latex','HorizontalAlignment', 'left')    
        end
    end
    a1bis=area([cs.cs2dis cs.cs2dis cs.c0 cs.c0],[0,limy,limy,0],...
        'LineStyle','-','DisplayName','no disease eliminated'); a1bis(1).FaceColor = gris;a1bis.FaceAlpha = 0.05; a1bis(1).EdgeColor = "none";%[0.95 0.95 0.95];
    %text((cs.cs1dis+cs.c0)/2,limy*0.9,'IIa','Interpreter','latex','FontSize',FontSizeGlobal,'HorizontalAlignment', 'center')
    text((cs.cs1dis+cs.c0)/2,limy*0.85,'\fontsize{22}{0} \selectfont IIa','Interpreter','latex','HorizontalAlignment', 'center')

    a0=area([cs.cs1dis cs.cs1dis lim_xsup lim_xsup],[0,limy,limy,0],...
        'LineStyle','-','DisplayName','no disease eliminated'); a0(1).FaceColor = gris;a0.FaceAlpha = 0.0; a0(1).EdgeColor = "none";% gris*3;
    %text((lim_xsup+cs.c0)/2,limy*0.9,'I','Interpreter','latex','FontSize',FontSizeGlobal,'HorizontalAlignment', 'center')
    text((lim_xsup+cs.c0)/2,limy*0.85,'\fontsize{22}{0} \selectfont I','Interpreter','latex','HorizontalAlignment', 'center')
    
    if(cs.cs2dis==cs.cs1dis)
        delete(a1)
    end
    lgd = legend([p_rhohat,p_rho1,p_rho2],'Interpreter','latex');
  %%  
    % print rhohat'
    %text(lim_xinf+ex,paramTab{1}.alpha+ey,['$\rho_{',paramTab{1}.disease,'}}^\prime$'],'Interpreter','latex','FontSize',FontSizeGlobal,'HorizontalAlignment', 'center')
    %text(lim_xinf+ex,paramTab{2}.alpha+ey,['$\rho_{',paramTab{2}.disease,'}^\prime$'],'Interpreter','latex','FontSize',FontSizeGlobal,'HorizontalAlignment', 'center')
    text('position',[lim_xinf+0.5*ex paramTab{2}.alpha+ey 0],'interpreter','latex','string',['\fontsize{22}{0}\selectfont$\rho$\fontsize{14}{0}\selectfont$_{',paramTab{2}.disease,'}^{\prime}$'],'HorizontalAlignment', 'left');
    text('position',[lim_xinf+0.5*ex paramTab{1}.alpha+ey 0],'interpreter','latex','string',['\fontsize{22}{0}\selectfont$\rho$\fontsize{14}{0}\selectfont$_{',paramTab{1}.disease,'}^{\prime}$'],'HorizontalAlignment', 'left');
    
end
if N==3
    lim_xinf = cs.cs3dis-(cs.c0-cs.cs3dis)/2;
    lim_xsup = cs.c0+(cs.c0-cs.cs3dis)/2;
    xlim([lim_xinf, lim_xsup])
    % zones 3 diseases
    a3=area([lim_xinf lim_xinf cs.cs3dis cs.cs3dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',[infections{cs.order(1)},', ',infections{cs.order(2)},' \& ', infections{cs.order(3)},' eliminated']); a3(1).FaceColor = gris;a3.FaceAlpha = 0.4;
    a2=area([cs.cs3dis cs.cs3dis cs.cs2dis cs.cs2dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',['only ',infections{cs.order(1)},' \& ', infections{cs.order(2)},' eliminated']); a2(1).FaceColor = gris;a2.FaceAlpha = 0.25;
    a1=area([cs.cs2dis cs.cs2dis cs.cs1dis cs.cs1dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',['only ', infections{cs.order(1)},' eliminated']); a1(1).FaceColor = gris;a1.FaceAlpha = 0.1;
    a0=area([cs.cs1dis cs.cs1dis lim_xsup lim_xsup],[0,limy,limy,0 ],'LineStyle','-','DisplayName','no disease eliminated'); a0(1).FaceColor = gris;a0.FaceAlpha = 0.0; a0(1).EdgeColor = gris*3;
    
    if(cs.cs2dis==cs.cs3dis)
        delete(a2)
    end
    if(cs.cs2dis==cs.cs1dis)
        delete(a1)
    end
    %title(['Testing rate in function of $c$ - ',['1:',modelTypes{1}],' 2:',modelTypes{2},' 3:',modelTypes{3}],'Interpreter','latex')
end
if N==4
    lim_xinf = cs.cs4dis-(cs.c0-cs.cs4dis)/2;
    lim_xsup = cs.c0+(cs.c0-cs.cs3dis)/2;
    xlim([lim_xinf, lim_xsup])
    % zones 3 diseases
    a4=area([lim_xinf lim_xinf cs.cs4dis cs.cs4dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName','HIV, syph., Ct \& Ng eliminated'); a4(1).FaceColor = gris;a4.FaceAlpha = 0.5;
    a3=area([cs.cs4dis cs.cs4dis cs.cs3dis cs.cs3dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',['only ',num2str(cs.order(1)),',',num2str(cs.order(2)),' \& ',num2str(cs.order(3)), ' eliminated']); a3(1).FaceColor = gris;a3.FaceAlpha = 0.4;
    a2=area([cs.cs3dis cs.cs3dis cs.cs2dis cs.cs2dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',['only ', num2str(cs.order(1)),' \& ',num2str(cs.order(2)), ' eliminated']); a2(1).FaceColor = gris;a2.FaceAlpha = 0.25;
    a1=area([cs.cs2dis cs.cs2dis cs.cs1dis cs.cs1dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',['only ', infections{1},' eliminated']); a1(1).FaceColor = gris;a1.FaceAlpha = 0.1;
    a0=area([cs.cs1dis cs.cs1dis lim_xsup lim_xsup],[0,limy,limy,0 ],'LineStyle','-','DisplayName','no disease eliminated'); a0(1).FaceColor = gris;a0.FaceAlpha = 0.0; a0(1).EdgeColor = gris*3;
    
    if(cs.cs3dis==cs.cs4dis)
        delete(a3)
    end
    if(cs.cs2dis==cs.cs3dis)
        delete(a2)
    end
    if(cs.cs2dis==cs.cs1dis)
        delete(a1)
    end
    
    %title(['Testing rate in function of $c$ - ',['1:',modelTypes{1}],' 2:',modelTypes{2},' 3:',modelTypes{3},' 4:',modelTypes{4}],'Interpreter','latex')
    %title(['Testing rate in function of $c$ - 1: Chlamydia, 2: Gonorrhea, 3:HIV, 4: Syphilis'],'Interpreter','latex')
end

%lmabels on the x axis
if N>2 %le faire manuellement
    %xticks([-1 tabcn.one(2) tabcn.one(1) 0 1])%order ?
    %xticklabels({'-1','$c_{',infectionsMini{2},'}^\prime$','$c_{',infectionsMini{1},'}^\prime$','0','1'})
    %xticks([-1 -0.5 -0.10 0 0.10 0.2 0.3 0.4 0.5 1])%order ?
    %xticklabels({'-1','-0.5','0','0.5','1'})
    ylim([0,limy])
end

%title(['Testing rate in function of $c$'],'Interpreter','latex')

%yticks([ 0 0.01 0.02 vecAlphas(cs.order(1)) 0.04 0.05 0.06 0.07 vecAlphas(cs.order(2)) 0.09 0.1 ])
%yticklabels({'0', '0.01', '0.02', ['$\rho_{',infectionsMini{cs.order(1)},'}^\prime$'],...
%     '0.04', '0.05', '0.06', '0.07', ['$\rho_{',infectionsMini{cs.order(2)},'}^\prime$'],'0.09', '0.1'})

%%
%xlabel("$c$","fontweight","bold",'FontSize',FontSizeGlobal,'Interpreter','latex')
xlabel('\fontsize{16}{0}\selectfont $c$',"fontweight","bold",'Interpreter','latex')
%ylabel("$\hat\rho$","fontweight","bold",'FontSize',FontSizeGlobal,'Interpreter','latex')
ylabel("\fontsize{16}{0}\selectfont$\hat\rho$","fontweight","bold",'Interpreter','latex')

ax.XAxis.FontSize = FontSizeGlobal;
ax.YAxis.FontSize = FontSizeGlobal;
ax.FontWeight = 'bold';
box on;
%legend boxoff;
%set(lgd, 'fontsize' , FontSizeGlobal , 'location' , 'east', 'color' , 'white' , 'box' , 'on','edgecolor','white' )
set(lgd, 'location' , 'east', 'color' , 'white' , 'box' , 'on','edgecolor','white' )

%lgd.FontSize = 15;
%set(findall(gcf,'-property','FontSize'),'FontSize',12)
