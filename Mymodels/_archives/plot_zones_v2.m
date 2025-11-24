%% plot
% functions plot
bleu = [0,0.4470,0.7410];jaune=[0.9290,0.6940,0.1250];rouge=[0.900,0.3250,0.0980];violet=[0.4940, 0.1840, 0.5560];vert=[0.4660, 0.6740, 0.1880];

figure()
p_rho1  = plot(vecC,tab.one(:,1),':','Linewidth',2,'Color',bleu,'DisplayName','$\hat\rho_1$'); %bleu/SEIIS 1
if (N>=2) %we plot rhohat for 1 disease
    hold on
    p_rho2   = plot(vecC,tab.one(:,2),':','Linewidth',2,'Color',jaune,'DisplayName','$\hat\rho_2$'); %jaune/SEIIS 1
    legend([p_rho1,p_rho2],'Interpreter','latex')
    if (N>=3)
        p_rho3   = plot(vecC,tab.one(:,3),':','Linewidth',2,'Color',rouge,'DisplayName','$\hat\rho_3$'); %rouge/SICR
        legend([p_rho1,p_rho2,p_rho3],'Interpreter','latex')
        if (N>=4)
            p_rho4   = plot(vecC,tab.one(:,4),':','Linewidth',2,'Color',vert,'DisplayName','$\hat\rho_4$'); %vert/ SEIIIS
            legend([p_rho1,p_rho2,p_rho3,p_rho4],'Interpreter','latex')
        end
    end

    %we plot combination of 2 diseases
    p_rho12  = plot(vecC,tab.two(:,1),'-.','Linewidth',1.5,'Color',vert,'DisplayName','$\hat\rho_{1\times 2}$');
    if N>=3
        p_rho13  = plot(vecC,tab.two(:,2),'-.','Linewidth',1.5,'Color',violet,'DisplayName','$\hat\rho_{1 \times 3}$');
        p_rho23  = plot(vecC,tab.two(:,4),'-.','Linewidth',1.5,'Color',(rouge+jaune)/2,'DisplayName','$\hat\rho_{2 \times 3}$');
        legend([p_rho1,p_rho2,p_rho3,p_rho4,p_rho12,p_rho13,p_rho23],'Interpreter','latex')
        if N>=4
            p_rho14  = plot(vecC,tab.two(:,3),'-.','Linewidth',1.5,'Color',(bleu+vert)/2,'DisplayName','$\hat\rho_{1 \times 4}$');
            p_rho24  = plot(vecC,tab.two(:,5),'-.','Linewidth',1.5,'Color',(vert+jaune)/2,'DisplayName','$\hat\rho_{2 \times 4}$');
            p_rho34  = plot(vecC,tab.two(:,6),'-.','Linewidth',1.5,'Color',(rouge+vert)/2,'DisplayName','$\hat\rho_{3 \times 4}$');
            legend([p_rho1,p_rho2,p_rho3,p_rho4,p_rho12,p_rho13,p_rho14,p_rho23,p_rho24,p_rho34],'Interpreter','latex')
        end
        %we plot combination of 3 diseases
        p_rho123  = plot(vecC,tab.three(:,1),'--','Linewidth',1.5,'Color',(violet+jaune)/3,'DisplayName','$\hat\rho_{1\times 2\times 3}$');
        if N>=4
            p_rho124  = plot(vecC,tab.three(:,2),'--','Linewidth',1.5,'Color',(bleu+jaune+vert)/3,'DisplayName','$\hat\rho_{1\times 2\times 4}$');
            p_rho134  = plot(vecC,tab.three(:,3),'--','Linewidth',1.5,'Color',(violet+vert)/3,'DisplayName','$\hat\rho_{1\times 3\times 4}$');
            p_rho234  = plot(vecC,tab.three(:,4),'--','Linewidth',1.5,'Color',(jaune+rouge+vert)/3,'DisplayName','$\hat\rho_{2\times 3\times 4}$');
            %we plot combination of 3 diseases
            p_rho1234  = plot(vecC,tab.four(:,1),'-','Linewidth',1.,'Color','k','DisplayName','$\hat\rho_{1\times 2\times 3 \times 4}$');
        end
    end
end
% 
% if(N==1)
%     p_rho1  = plot(vecC,tab.one,'-','Linewidth',2,'Color',bleu,'DisplayName','$\hat\rho$'); %bleu
% end
% if (N>=2)
%     p_rho1  = plot(vecC,tab.one(:,1),':','Linewidth',2,'Color',bleu,'DisplayName','$\hat\rho_1$'); %bleu
%     hold on
%     p_rho2   = plot(vecC,tab.one(:,2),':','Linewidth',2,'Color',jaune,'DisplayName','$\hat\rho_2$'); %jaune
%     legend([p_rho1,p_rho2],'Interpreter','latex')
%     p_rho12  = plot(vecC,tab.two(:,1),'--','Linewidth',2,'Color',vert,'DisplayName','$\hat\rho_{1\times 2}$');
% end
% if (N>=3)
%     p_rho3   = plot(vecC,tab.one(:,3),':','Linewidth',2,'Color',rouge,'DisplayName','$\hat\rho_3$'); %rouge
%     legend([p_rho1,p_rho2,p_rho3,p_rho12],'Interpreter','latex')
%     p_rho13  = plot(vecC,tab.two(:,2),'--','Linewidth',2,'Color',violet,'DisplayName','$\hat\rho_{1 \times 3}$');
%     p_rho23  = plot(vecC,tab.two(:,3),'--','Linewidth',2,'Color',(rouge+jaune)/2,'DisplayName','$\hat\rho_{2 \times 3}$');
%     p_rho_i  = plot(vecC((ics1+1):end),tab.rhohat((ics1+1):end),'-','Linewidth',1.2, 'Color', [0.2, 0.2, 0.2], 'DisplayName','$\hat\rho$');
%     p_rho_ii = plot(vecC((ics2+1):ics1),tab.rhohat((ics2+1):ics1),'-','Linewidth',1.2, 'Color', [0.2, 0.2, 0.2],'HandleVisibility','off');
%     p_rho_iii= plot(vecC((ics3+1):ics2),tab.rhohat((ics3+1):ics2),'-','Linewidth',1.2, 'Color', [0.2, 0.2, 0.2],'HandleVisibility','off');
%     p_rho_iv = plot(vecC(1:ics3),tab.rhohat(1:ics3),'-','Linewidth',1.2, 'Color', [0.2, 0.2, 0.2],'HandleVisibility','off');
%     %legend([p_rho12,p_rho13,p_rho23],'Interpreter','latex')
% end

modelTypes = repelem({'SEIIS', 'SICR', 'SEIIIS'}, [nSEIIS nSICR nSEIIIS]);

title(['Testing rate in function of $c$'],'Interpreter','latex')
xlabel("$c$","fontweight","bold",'Interpreter','latex')
ylabel("voluntary testing rate $\hat\rho$","fontweight","bold",'Interpreter','latex')
limy=1.2*maxalpha;
ylim([0,limy])

% areas plot
%N = nSEIIS+nSICR+nSEIIIS;
gris=[78/255, 78/255, 78/255];

if N==1
    figure(1)
    hold on
    lim_xinf = cs.cs1dis-(cs.c0-cs.cs1dis)/2;
    lim_xsup = cs.c0+(cs.c0-cs.cs1dis)/2;
    xlim([lim_xinf, lim_xsup])
    a2=area([lim_xinf lim_xinf cs.cs1dis cs.cs1dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName','eliminated'); a2(1).FaceColor = gris;a2.FaceAlpha = 0.4;  
    a1=area([cs.cs1dis cs.cs1dis cs.c0 cs.c0],[0,limy,limy,0 ],'LineStyle','none','DisplayName','controlled'); a1(1).FaceColor = gris;a1.FaceAlpha = 0.15;  
    a0=area([cs.c0 cs.c0 lim_xsup lim_xsup],[0,limy,limy,0 ],'LineStyle','none','DisplayName','no change'); a0(1).FaceColor = gris;a0.FaceAlpha = 0.05;  
end
if N==2
    lim_xinf = cs.cs2dis-(cs.c0-cs.cs2dis)/2;
    lim_xsup = cs.c0+(cs.c0-cs.cs2dis)/2;
    xlim([lim_xinf, lim_xsup])
    % zones 2 diseases
    a2=area([lim_xinf lim_xinf cs.cs2dis cs.cs2dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName','1 \& 2 eliminated'); a2(1).FaceColor = gris;a2.FaceAlpha = 0.4;  
    a1=area([cs.cs2dis cs.cs2dis cs.cs1dis cs.cs1dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',['only ', num2str(cs.order(1)),' eliminated']); a1(1).FaceColor = gris;a1.FaceAlpha = 0.15;  
    a0=area([cs.cs1dis cs.cs1dis lim_xsup lim_xsup],[0,limy,limy,0 ],'LineStyle','none','DisplayName','no disease eliminated'); a0(1).FaceColor = gris;a0.FaceAlpha = 0.05;  
    if(cs.cs2dis==cs.cs1dis)
        delete(a1)
    end
end
if N==3
    lim_xinf = cs.cs3dis-(cs.c0-cs.cs3dis)/2;
    lim_xsup = cs.c0+(cs.c0-cs.cs3dis)/2;
    xlim([lim_xinf, lim_xsup])
    % zones 3 diseases
    a3=area([lim_xinf lim_xinf cs.cs3dis cs.cs3dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName','1,2 \& 3 eliminated'); a3(1).FaceColor = gris;a3.FaceAlpha = 0.4;  
    a2=area([cs.cs3dis cs.cs3dis cs.cs2dis cs.cs2dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',['only ', num2str(cs.order(1)),' \& ',num2str(cs.order(2)), ' eliminated']); a2(1).FaceColor = gris;a2.FaceAlpha = 0.25;  
    a1=area([cs.cs2dis cs.cs2dis cs.cs1dis cs.cs1dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',['only ', num2str(cs.order(1)),' eliminated']); a1(1).FaceColor = gris;a1.FaceAlpha = 0.15;  
    a0=area([cs.cs1dis cs.cs1dis lim_xsup lim_xsup],[0,limy,limy,0 ],'LineStyle','none','DisplayName','no disease eliminated'); a0(1).FaceColor = gris;a0.FaceAlpha = 0.05;  
    
    if(cs.cs2dis==cs.cs3dis)
        delete(a2)
    end

    if(cs.cs2dis==cs.cs1dis)
        delete(a1)
    end
    title(['Testing rate in function of $c$ - ',['1:',modelTypes{1}],' 2:',modelTypes{2},' 3:',modelTypes{3}],'Interpreter','latex')

end
if N==4
    lim_xinf = cs.cs4dis-(cs.c0-cs.cs4dis)/2;
    lim_xsup = cs.c0+(cs.c0-cs.cs3dis)/2;
    xlim([lim_xinf, lim_xsup])
    % zones 3 diseases
    a4=area([lim_xinf lim_xinf cs.cs4dis cs.cs4dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName','1,2,3 \& 4 eliminated'); a4(1).FaceColor = gris;a4.FaceAlpha = 0.5;  
    a3=area([cs.cs4dis cs.cs4dis cs.cs3dis cs.cs3dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',['only ',num2str(cs.order(1)),',',num2str(cs.order(2)),' \& ',num2str(cs.order(3)), ' eliminated']); a3(1).FaceColor = gris;a3.FaceAlpha = 0.4;  
    a2=area([cs.cs3dis cs.cs3dis cs.cs2dis cs.cs2dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',['only ', num2str(cs.order(1)),' \& ',num2str(cs.order(2)), ' eliminated']); a2(1).FaceColor = gris;a2.FaceAlpha = 0.25;  
    a1=area([cs.cs2dis cs.cs2dis cs.cs1dis cs.cs1dis],[0,limy,limy,0 ],'LineStyle','none','DisplayName',['only ', num2str(cs.order(1)),' eliminated']); a1(1).FaceColor = gris;a1.FaceAlpha = 0.15;  
    a0=area([cs.cs1dis cs.cs1dis lim_xsup lim_xsup],[0,limy,limy,0 ],'LineStyle','none','DisplayName','no disease eliminated'); a0(1).FaceColor = gris;a0.FaceAlpha = 0.05;  
    
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
    title(['Testing rate in function of $c$ - 1: Chlamydia, 2: Gonorrhea, 3:HIV, 4: Syphilis'],'Interpreter','latex')

end