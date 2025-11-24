function randomU_SISSIS_V1
fh = figure;
bh1 = uicontrol(fh,'Position',[700 200 100 50],...
                     'String','Aleatoire',...
                     'Callback',@button1_plot);
 

ah1 = axes('Parent',fh,'units','pixels',...
           'Position',[150 150 700 450]);

 %------------------------------------------------
     function button1_plot(hObject,eventdata)
        cond = 0;
        it = 0;
        while(~cond & it<1000)
            [beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');  
            cond =  plotU_SISSIS(beta1,beta2,gamma1,gamma2,s1,s2,b,mu);
            it=it+1
        end
     end
end

function cond = plotU_SISSIS (beta1,beta2,gamma1,gamma2,s1,s2,b,mu)
    R1 = beta1/(mu+gamma1);
    R2 = beta2/(gamma2+mu);
    alpha1 = beta1/s1*(1-1/R1);
    alpha2 = beta2/s2*(1-1/R2);
    rho1 = alpha1/2;
    rho2 = alpha2/2;
    
    vecDelta = 0:0.005:(max(alpha1, alpha2)+0.2);
    vecGamma1 = gamma1 + s1*vecDelta;
    vecGamma2 = gamma2 + s2*vecDelta;
    [U,vecPrev, vecR10,vecR20,P1,P2,P12] = utilityFunctionSIS2(beta1,beta2,vecGamma1,vecGamma2,vecDelta,b,mu,beta1,beta2);
    hold off
    
    [Umax,imax] = max(U);
    rhomax = vecDelta(imax);
    
    [U12max,i12max] = max(vecDelta.*P12);
    rho12max = vecDelta(i12max);
    
    if (alpha1<alpha2)
        rhol = rho1;
        rhok = rho2;
    else
        rhol = rho2;
        rhok = rho1;
    end
    
    ordreOK = rhol <= rho12max & rho12max<= rhok;
    ordre3OK = rhok <= rho12max & rho12max<= rhol;
    ordre2Ok = rhol - rho12max > 0.005;
    ordre4Ok = -rhok + rho12max > 0.005;
    
    cond = ordre4Ok
    
    if (cond )
    figure(1)
    plot(vecDelta, max(U,0), 'LineWidth',2)
    title([{'Utility of the SISxSIS model'},...
        {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
        num2str(gamma2), ' s1=', num2str(s1),' s2=', num2str(s2), '\mu=', num2str(mu), ' \pi=',  num2str(b)]}, ...
        {[' R_0^1=' num2str(round(beta1/(gamma1+mu),2)), ' R_0^2=' num2str(round(beta2/(gamma2+mu),2)),...
        ' \alpha_1=', num2str(round(beta1/s1*(1-1/R1),2)),...
        ' \alpha_2=', num2str(round(beta2/s2*(1-1/R2),2)), ...
        ' \rho_1=', num2str(round(rho1,2)), ' \rho_2=', num2str(round(rho2,2)), ' \rho_{12}=', num2str(round(rhomax,2)), ...
        ' Umax=', num2str(round(Umax,2)) ]},...
        {['\alpha_l=', num2str(round(min(alpha1, alpha2),2)), ' \rho_l=',num2str(round(rhol,3)), ...
         ' \rho_k=',num2str(round(rhok,3)), ' \rho_{12}^*=',num2str(round(rho12max,3)), ' \rho_{12}^* \in [\rho_l, \rho_k]:', num2str(ordreOK),' \rho_{12}^* \in [\rho_k, \rho_l]:', num2str(ordre3OK), ...
         '  \rho_{12}^* < \rho_l:', num2str(ordre2Ok), '  \rho_{12}^* > \rho_k:', num2str(ordre4Ok) ]} ])
    xlabel("Voluntary-testing rate \delta","fontweight","bold")
    ylabel("U(\delta)","fontweight","bold")
    hold on;
    
    plot(vecDelta, max(vecDelta.*max(P1,0),0),  '-.','LineWidth',3)
    plot(vecDelta, max(vecDelta.*max(P2,0),0),  '--','LineWidth',3)
    plot(vecDelta, max(vecDelta.*max(P12,0),0),  '.','MarkerSize',10)
%   plot(vecDelta, vecDelta.*(max(P1,0)+max(P2,0)))
%plot(vecDelta, vecDelta.*(P1+P2))

    legend('U', 'U_1', 'U_2', 'U_{12}', 'U_1+U_2')
    end
end