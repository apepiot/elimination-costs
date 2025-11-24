function randomU_SIRSIS
fh = figure;
bh1 = uicontrol(fh,'Position',[880 200 100 50],...
                     'String','Aleatoire',...
                     'Callback',@button1_plot);
 
%m1 = uicontrol(fh, 'Position',[700 350 50 50], 'Style','edit',...
%    'String', 'beta1','Callback',@manuel1)

ah1 = axes('Parent',fh,'units','pixels',...
           'Position',[150 150 700 450]);

 %------------------------------------------------
     function button1_plot(hObject,eventdata)
        cond = 0;%mettre a 0
        it=0;
        while(~cond & it<10000)%mettre a 1000
            [beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');  
            cond = plotU (beta1,beta2,gamma1,gamma2,s1,s2,b,mu);
            it=it+1
        end
     end

    function manuel1(beta1input, hObject,eventdata)
       
            [beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');  
            beta1 = str2double(beta1input.String);
            %beta2 = str2double(beta2input.String);

            cond = plotU (beta1,beta2,gamma1,gamma2,s1,s2,b,mu)
       
    end
end

function cond = plotU (beta1,beta2,gamma1,gamma2,s1,s2,b,mu)
    R1 = beta1/(mu+gamma1);
    R2 = beta2/(gamma2+mu);
    alpha1 = beta1/s1*(1-1/R1);
    alpha2 = beta2/s2*(1-1/R2);
    rho1 = beta1/s1*(1/sqrt(R1) - 1/R1);
    rho2 = alpha2/2;
    
    vecDelta =0:0.05:(max(alpha1, alpha2)+0.2);
    [U,vecPrev, vecR10,vecR20,P1,P2,P12,rho12max] = utilityFunctionSIRxSIS(beta1,beta2,gamma1,gamma2,vecDelta,vecDelta,s1,s2,b,mu);
    hold off
    
    [Umax,imax] = max(U);
    rhomax = vecDelta(imax);
    
    %[U12max,i12max] = max(vecDelta.*P12);
    %rho12max = vecDelta(i12max);
    
    if (alpha1<alpha2)
        rhol = rho1;
        rhok = rho2;
        alphal = alpha1;
        alphak = alpha2;
    else
        rhol = rho2;
        rhok = rho1;
        alphal = alpha2;
        alphak = alpha1;
    end
    
    ordreOK = rhol <= rho12max & rho12max<= rhok;
    ordre3OK = rhok <= rho12max & rho12max<= rhol;
    minrho = min(rhol,rhok);
    ordre2Ok = minrho - rho12max > 0.005;
    maxrho = max(rhol,rhok);
    ordre4Ok = -maxrho + rho12max > 0.005;
    
    %cond = (rhok<rhol) & (rho12max>rhol);
    cond = ( alphal>rhok & alphal<rho12max);
    
    if(cond)
        
    figure(1)
    plot(vecDelta, max(U,0), 'LineWidth',2)
    title([{'Utility of the SISxSIR model'},...
        {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
        num2str(gamma2), ' s1=', num2str(s1),' s2=', num2str(s2), '\mu=', num2str(mu), ' \pi=',  num2str(b)]}, ...
        {[' R_0^1=' num2str(round(beta1/(gamma1+mu),2)), ' R_0^2=' num2str(round(beta2/(gamma2+mu),2)),...
        ' \alpha_1=', num2str(round(beta1/s1*(1-1/R1),2)),...
        ' \alpha_2=', num2str(round(beta2/s2*(1-1/R2),2)), ...
        ' \rho_1=', num2str(round(rho1,2)), ' \rho_2=', num2str(round(rho2,2)), ' \rho_{12}=', num2str(round(rhomax,2)), ...
        ' Umax=', num2str(round(Umax,2)) ]},...
        {['\alpha_l=', num2str(round(min(alpha1, alpha2),2)), ' \rho_l=',num2str(round(rhol,3)), ...
         ' \rho_k=',num2str(round(rhok,3)), ' \rho_{12}^*=',num2str(round(rho12max,3))]},...
        {[' \rho_{12}^* \in [\rho_l, \rho_k]:', num2str(ordreOK),' \rho_{12}^* \in [\rho_k, \rho_l]:', num2str(ordre3OK) ...
         '  \rho_{12}^* < min(\rho_l, \rho_k):', num2str(ordre2Ok), '  \rho_{12}^* > max(\rho_l, \rho_k):', num2str(ordre4Ok) ]} ])
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