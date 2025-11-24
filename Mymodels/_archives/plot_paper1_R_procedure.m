
function [plt] = plot_paper1_R_procedure(vecR, c1, c2, vecC,R0)
%procedure to plot same graphs for the paper 1
    close all;
    %figure()
    %upbound = (c2+(c2-c1)/2);
    limit_y = 1.2*R0;
    a=area([c2 c2 (c2+(c2-c1)/2) (c2+(c2-c1)/2) ],[0 limit_y limit_y 0],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 1;
    hold on
    plt=plot(c1*ones(1,100), 0:limit_y/99:limit_y ,'--','Color',[0.6,0.6,0.6],'LineWidth',3)
    hold on
    plot(c2*ones(1,100), 0:limit_y/99:limit_y ,'--','Color',[0.6,0.6,0.6],'LineWidth',3)
    %plot hatrho
    vecR11 = vecR(vecC<=c2); vecC1 = vecC(vecC<=c2);
    plot(vecC1,R0*ones(length(vecC1),1),'--','Color',[0.6,0.6,0.6],'LineWidth',3)
    plot(vecC1,ones(length(vecC1),1),'--','Color',[0.6,0.6,0.6],'LineWidth',3)

    plot(vecC1,vecR11,'-','LineWidth',5,'Color',[44/255, 63/255, 81/255])
    xlim([(c1-(c2-c1)/2),(c2+(c2-c1)/2)])
    ylim([0.75,limit_y])
    ey = 0.05*R0; 
    ex = (vecC(end)-vecC(1))*0.075;
    text(c1,0.75-ey,'$c_1$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(c2,0.75-ey,'$c_2$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((c1-(c2-c1)/2)-ex,1,'$1$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((c1-(c2-c1)/2)-ex,R0,'$\mathtt R(0)$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')

    set(gca,'YTickLabel',[],'XTickLabel',[]);
    text((c2+(c2-c1)/2),0.75-ey,'$c$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((c1-(c2-c1)/2)-ex,limit_y,'$\mathtt R(\hat\rho)$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(0,0.78-ey*1.5,'$0$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')

    %ZONES I II III
    text((c1+vecC1(1))/2,0.9*limit_y,'I','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((c1+c2)/2,0.9*limit_y,'II','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((c2+(c2+(c2-c1)/2))/2,0.9*limit_y,'III','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')

end