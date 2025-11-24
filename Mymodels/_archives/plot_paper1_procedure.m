
function [] = plot_paper1_procedure(vecRhomax1, c1, c2, alpha, vecC)
%procedure to plot same graphs for the paper 1
    close all;
    %figure()
    %upbound = (c2+(c2-c1)/2);
    limit_y = 1.3*alpha;
    a=area([c2 c2 (c2+(c2-c1)/2) (c2+(c2-c1)/2) ],[0 limit_y limit_y 0],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;
    hold on
    plot(c1*ones(1,100), 0:limit_y/99:limit_y ,'--','Color',[0.6,0.6,0.6],'LineWidth',3)
    hold on
    plot(c2*ones(1,100), 0:limit_y/99:limit_y ,'--','Color',[0.6,0.6,0.6],'LineWidth',3)
    %plot hatrho
    vecRhomax11 = vecRhomax1(vecC<=c2); vecC1 = vecC(vecC<=c2);
    plot(vecC1,alpha*ones(length(vecC1),1),'--','Color',[0.6,0.6,0.6],'LineWidth',3)
    %vecRhomax12 = vecRhomax1(vecC>=c2); vecC2 = vecC(vecC>=c2);
    plot(vecC1,vecRhomax11,'-','LineWidth',8,'Color',[44/255, 63/255, 81/255])
    %plot(vecC2,vecRhomax12,'--','LineWidth',4,'Color','blue')
    xlim([(c1-(c2-c1)/2),(c2+(c2-c1)/2)])
    ylim([0,limit_y])
    ey = 0.05*alpha; 
    ex = (vecC(end)-vecC(1))*0.05;
    text(c1,0-ey,'$c_1$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(c2,0-ey,'$c_2$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((c1-(c2-c1)/2)-ex,alpha,'$\rho\prime$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    set(gca,'YTickLabel',[],'XTickLabel',[]);
    text((c2+(c2-c1)/2),0-ey,'$c$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((c1-(c2-c1)/2)-ex,limit_y,'$\hat\rho$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(0,0-ey*1.5,'$0$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')

    %ZONES I II III
    text((c1+vecC1(1))/2,0.85*limit_y,'I','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((c1+c2)/2,0.85*limit_y,'II','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((c2+(c2+(c2-c1)/2))/2,0.85*limit_y,'III','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')

end