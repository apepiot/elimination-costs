
function [pltrho] = plot_paper2_procedure1(vecRhomax1,c1,c2,c0,alpha1,alpha2,vecC,c11,c12,c21,c22,cmaxalpha,cswitch)
%procedure to plot same graphs for the paper 1
    %close all;
    %figure()
    maxalpha=max(alpha1,alpha2);
    minalpha=min(alpha1,alpha2);
    limit_y = 1.3*maxalpha;
    a1=area([(cmaxalpha-(c0-cmaxalpha)/2) (cmaxalpha-(c0-cmaxalpha)/2) cmaxalpha cmaxalpha],[0 limit_y limit_y 0],'LineStyle','none'); a1(1).FaceColor = [132/255, 151/255, 176/255];a1.FaceAlpha = 0.5;
    hold on
    a2=area([cmaxalpha cmaxalpha c2 c2],[0 limit_y limit_y 0],'LineStyle','none'); a2(1).FaceColor = [132/255, 151/255, 176/255];a2.FaceAlpha = 0.2;
    hold on
    
    %lignes verticales c1 et c2
    plot(cmaxalpha*ones(1,100), 0:limit_y/99:limit_y ,':','Color','k','LineWidth',1.5)
    hold on
    plot(c2*ones(1,100), 0:limit_y/99:limit_y ,':','Color','k','LineWidth',1.)
    plot(c0*ones(1,100), 0:limit_y/99:limit_y ,':','Color','k','LineWidth',1.5)
    
    %lignes verticales c11 c12 c21 c22
    %plot(c11*ones(1,100), 0:limit_y/99:limit_y ,':','Color','k','LineWidth',1)
    %plot(c12*ones(1,100), 0:limit_y/99:limit_y ,':','Color','k','LineWidth',1)
    %plot(c21*ones(1,100), 0:limit_y/99:limit_y ,':','Color','k','LineWidth',1)
    %plot(c22*ones(1,100), 0:limit_y/99:limit_y ,':','Color','k','LineWidth',1)
    plot(cswitch*ones(1,100), 0:limit_y/99:limit_y ,':','Color','k','LineWidth',1)
    
    %plot hatrho
    vecC1 = vecC(vecC<=cmaxalpha);
    vecC2 = vecC(vecC<=cmaxalpha);
    
    %ligne horizontale a rho'
    plot(vecC1,maxalpha*ones(length(vecC1),1),':','Color','k','LineWidth',1.5)
    plot(vecC2,minalpha*ones(length(vecC2),1),':','Color','k','LineWidth',1.5)
    vecRhomax12 = vecRhomax1(vecC>=cmaxalpha); vecC0 = vecC(vecC>=cmaxalpha);
    
    %plot rhohat
    pltrho = plot(vecC0,vecRhomax12,'-','LineWidth',4,'Color',[44/255, 63/255, 81/255],'DisplayName','$\hat\rho(c)$')
    %plot(vecC0,vecRhomax12,'--','LineWidth',4,'Color','blue')
    xlim([cmaxalpha-(c0-cmaxalpha)/2,(c0+(c0-cmaxalpha)/2)]) %(c1-(c0-c1)/2)
    ylim([0,limit_y])
    ey = 0.05*maxalpha; 
    ex = abs(maxalpha-vecC(1))*0.005;
    
    %abscisse
    text(c1,0-ey,'$c_2^{12}$','Interpreter','latex','FontSize',15, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(c0,0-ey,'$c_0^{12}$','Interpreter','latex','FontSize',15, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(c2,0-ey,'$c_1^{12}$','Interpreter','latex','FontSize',15, 'FontWeight','bold','HorizontalAlignment', 'center')

    text(c11,0-ey,'$c_1^1$','Interpreter','latex','FontSize',15, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(c12,0-ey,'$c_0^1$','Interpreter','latex','FontSize',15, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(c21,0-ey,'$c_{2}^2$','Interpreter','latex','FontSize',15, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(c22,0-ey,'$c_0^2$','Interpreter','latex','FontSize',15, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(cswitch,0-ey,'$c_s$','Interpreter','latex','FontSize',15, 'FontWeight','bold','HorizontalAlignment', 'center')
    
    %ordonnee
    text(vecC(1)-ex,alpha1,'$\rho_1\prime$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(vecC(1)-ex,alpha2,'$\rho_2\prime$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    
    set(gca,'YTickLabel',[],'XTickLabel',[]);
    text((c0+(c0-cmaxalpha)/2),0-ey,'$c$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(vecC(1)-ex,limit_y,'$\hat\rho$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')
    text(0,0-ey*1.5,'$0$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')

    %ZONES I II III
    %text((c1+min(vecC))/2,0.85*limit_y,'IV','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((cmaxalpha-(c0-cmaxalpha)/2+cmaxalpha)/2,0.85*limit_y,'III','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((cmaxalpha+c2)/2,0.85*limit_y,'II','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
    text((c2+vecC(end))/2,0.85*limit_y,'I','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')

end