%% Plot Utility
function [outputArg1,outputArg2] = plot_U(paramTab)
if size(paramTab,2)==4
    vecAlpha = [paramTab{1}.alpha,paramTab{2}.alpha,paramTab{3}.alpha,paramTab{4}.alpha];
    [orderAlpha,elimOrder] = sort(vecAlpha);
    vecRho=0:0.005:0.2; vecRho=[vecRho,vecAlpha]; vecRho=sort(vecRho);
    
    U1234 = U_SEIIS2SICRSEIIIS_v4(paramTab{1},paramTab{2},paramTab{3},paramTab{4},mu,10,vecRho,c,f);
    U123 = U_SEIIS2SICR_v4(paramTab{1},paramTab{2},paramTab{3},mu,10,vecRho,c,f);
    U234 = U_SEIISSICRSEIIIS_v4(paramTab{2},paramTab{3},paramTab{4},mu,10,vecRho,c,f);
    U134 = U_SEIISSICRSEIIIS_v4(paramTab{1},paramTab{3},paramTab{4},mu,10,vecRho,c,f);
    U124 = U_SEIIS2SEIIIS_v4(paramTab{1},paramTab{2},paramTab{4},mu,10,vecRho,c,f);
    U13 = U_SEIISSICR_v4(paramTab{1},paramTab{3},mu,10,vecRho,c,f);
    U12 = U_SEIIS2_v4(paramTab{1},paramTab{2},mu,10,vecRho,c,f);
    U23 = U_SEIISSICR_v4(paramTab{2},paramTab{3},mu,10,vecRho,c,f);
    U14 = U_SEIISSEIIIS_v4(paramTab{1},paramTab{4},mu,10,vecRho,c,f);
    U24 = U_SEIISSEIIIS_v4(paramTab{2},paramTab{4},mu,10,vecRho,c,f);
    U34 = U_SICRSEIIIS_v4(paramTab{3},paramTab{4},mu,10,vecRho,c,f);
    U4 = U_SEIIIS_v4(paramTab{4},mu,10,vecRho,c,f);
    U3 = U_SICR_v4(paramTab{3},mu,10,vecRho,c,f);
    U2 = U_SEIIS_v4(paramTab{2},mu,10,vecRho,c,f);
    U1 = U_SEIIS_v4(paramTab{1},mu,10,vecRho,c,f);
    
    
    
    
    p     = plot(vecRho(1),U1234(1),'-','Color','black','DisplayName', 'U', 'Linewidth',2);
    hold on
    p1234 = plot(vecRho,U1234,'-','Color',[0.8500 0.3250 0.0980],'DisplayName', 'U_{hscg}', 'Linewidth',2);
    
    if elimOrder(1)==4 %syphilis eliminated first
        p123 = plot(vecRho(vecRho>=paramTab{4}.alpha),U123(vecRho>=paramTab{4}.alpha),'-','DisplayName', 'U_{cgh}', 'Linewidth',2,'Color',[0.9290 0.6940 0.1250]); %
        p123_part2 = plot(vecRho(vecRho<=paramTab{4}.alpha),U123(vecRho<=paramTab{4}.alpha),'--','DisplayName', 'U_{hcg}', 'Linewidth',2,'Color',[0.9290 0.6940 0.1250]); %
        if elimOrder(2)==3 %then HIV
            p12  = plot(vecRho(vecRho>=paramTab{3}.alpha),U12(vecRho>=paramTab{3}.alpha),'-','DisplayName', 'U_{cg}', 'Linewidth',2,'Color',[0.4660 0.6740 0.1880]); %CtxNg : vert
            p12_part2  = plot(vecRho(vecRho<=paramTab{3}.alpha),U12(vecRho<=paramTab{3}.alpha),'--','DisplayName', 'U_{cg}', 'Linewidth',2,'Color',[0.4660 0.6740 0.1880]); %CtxNg : vert
            if elimOrder(3)==2
                p1   = plot(vecRho,U1,'-','DisplayName', 'U_{c}', 'Linewidth',2); %Ct : jaune
            elseif elimOrder(3)==1
                p2   = plot(vecRho(vecRho>=paramTab{1}.alpha),U2(vecRho>=paramTab{1}.alpha),'-','DisplayName', 'U_{g}', 'Linewidth',2,'Color',[0 0.4470 0.7410]); %Ng : bleu
                p2_part2   = plot(vecRho(vecRho<=paramTab{1}.alpha),U2(vecRho<=paramTab{1}.alpha),'--','DisplayName', 'U_{g}', 'Linewidth',2,'Color',[0 0.4470 0.7410]); %Ng : bleu
            end
        elseif elimOrder(2)==2 %then Ng
            p13  = plot(vecRho,U13,'-','DisplayName', 'U_{ch}', 'Linewidth',2);
            if elimOrder(3)==3 %then HIV
                p1   = plot(vecRho,U1,'-','DisplayName', 'U_{c}', 'Linewidth',2); %Ct : jaune
            elseif elimOrder(3)==1 %then Ct
                p3   = plot(vecRho,U3,'-','DisplayName', 'U_{c}', 'Linewidth',2); %HIV : rouge
            end
        elseif elimOrder(2)==1 %then Ct
            p23  = plot(vecRho,U23,'-','DisplayName', 'U_{gh}', 'Linewidth',2);
            if elimOrder(3)==3 %then HIV
                p2   = plot(vecRho(vecRho>=paramTab{1}.alpha),U2(vecRho>=paramTab{1}.alpha),'-','DisplayName', 'U_{g}', 'Linewidth',2,'Color',[0 0.4470 0.7410]); %Ng : bleu
                p2_part2   = plot(vecRho(vecRho<=paramTab{1}.alpha),U2(vecRho<=paramTab{1}.alpha),'--','DisplayName', 'U_{g}', 'Linewidth',2,'Color',[0 0.4470 0.7410]); %Ng : bleu
            elseif elimOrder(3)==2 %then Ng
                p3   = plot(vecRho,U3,'-','DisplayName', 'U_{c}', 'Linewidth',2); %HIV : rouge
            end
        end
    elseif elimOrder(1)==3 %hiv eliminated first
        p124 = plot(vecRho,U124,'-','DisplayName', 'U_{cgs}', 'Linewidth',2);
        if elimOrder(2)==4
            %p12
            if elimOrder(3)==1
                %p2
            elseif elimOrder(3)==2
                %p1
            end
        elseif elimOrder(2)==2
            %p14
            if elimOrder(3)==4
                %p1
            elseif elimOrder(3)==1
                %p4
            end
        elseif elimOrder(2)==1
            p24  = plot(vecRho,U24,'-','DisplayName', 'U_{gs}', 'Linewidth',2);
            if elimOrder(3)==4
                %p2
            elseif elimOrder(3)==2
                %p4
            end
        end
    elseif elimOrder(1)==2
        p134 = plot(vecRho,U134,'-','DisplayName', 'U_{chs}', 'Linewidth',2);
        if elimOrder(2)==4
            %p13
            if elimOrder(3)==3
                %p1
            elseif elimOrder(3)==1
                %p3
            end
        elseif elimOrder(2)==3
            p14  = plot(vecRho,U14,'-','DisplayName', 'U_{cs}', 'Linewidth',2);
            if elimOrder(3)==4
                %p1
            elseif elimOrder(3)==1
                %p4
            end
        elseif elimOrder(2)==1
            p34  = plot(vecRho,U34,'-','DisplayName', 'U_{hs}', 'Linewidth',2);
            if elimOrder(3)==4
                %p3
            elseif elimOrder(3)==3
                %p4
            end
        end
    elseif elimOrder(1)==1
        p234 = plot(vecRho,U234,'-','DisplayName', 'U_{ghs}', 'Linewidth',2);
        if elimOrder(2)==4
            %p23
            if elimOrder(3)==2
                %p3
            elseif elimOrder(3)==3
                %p2
            end
        elseif elimOrder(2)==3
            p24  = plot(vecRho,U24,'-','DisplayName', 'U_{gs}', 'Linewidth',2);
            if elimOrder(3)==4
                %p2
            elseif elimOrder(3)==2
                %p4
            end
        elseif elimOrder(2)==2
           %p34
           if elimOrder(3)==4
               %p3
           elseif elimOrder(3)==3
               p4   = plot(vecRho,U4,'-','DisplayName', 'U_{s}', 'Linewidth',2);
           end
        end
    end

    ylim([0,max(U1234)*1.2])
    xlim([0,max(vecRho)])
    legend([p,p1234,p123,p12,p2])
    
    chH = get(gca,'Children');
    set(gca,'Children',[chH([1:2,4,6:8]);chH([3,5])])
end

end

