function [bestStrategy,worstStrategy,costStratTable,cleft,cright,costStratTable_sorted,tab_refined,tab,costElimTable] = findBestStrategy_v2(cleft,cright,N,paramTab,mu,b,biasFactor,optFindRhohat)
%% This section gives rhohat=argmaxU for a given model for several values of c
vecC = linspace(cleft,cright,optFindRhohat.sampleSize);
[tab,~,tabcn,~] = findRhohat_v5(N,paramTab,mu,b,vecC,biasFactor);
tab_refined = tab;
%vecC = tab.c;
vecAlphas = [paramTab{1}.alpha, paramTab{2}.alpha,paramTab{3}.alpha,paramTab{4}.alpha];
[costElimTable] = createCostElimTable(tab,tabcn,vecAlphas);
%figure(1)
%plot_4dis(vecC,tab)

precision = optFindRhohat.Tol*0.1;
%% TO DO :
%1. on calcule pour tous les modeles
%si on n'a pas trouve la borne pour eliminer toutes les diseases a
%gauche et la premiere (a droite), il faut elargir l'intervalle de
%recherche

%2. refining vecC if some strategies are almost the same
%(c_elim_mod1-c_elim_mod2)<sampleStep

%% thresholds of cost for disease elimination
if N==4
    %     if (vecC(1)>min(tabcn.one)) %on verifie qu'on calule bien sur les bornes min(cnn) et max(cnn)
    %         newC = min(tabcn.one);
    %         [newTab,~,~,~] = findRhohat_v5(N,paramTab,mu,b,newC,biasFactor);
    %         tab_refined = concatenateTabRhohat(tab_refined,newTab);
    %         [costElimTable] = createCostElimTable(tab_refined,tabcn,vecAlphas);
    %         vecC = tab_refined.c;
    %         figure(1)
    %         plot_4dis(vecC,tab_refined)
    %     end
    %     if (vecC(end)<max(tabcn.one)) %on verifie qu'on calule bien sur les bornes min(cnn) et max(cnn)
    %         newC = max(tabcn.one);
    %         [newTab,~,~,~] = findRhohat_v5(N,paramTab,mu,b,newC,biasFactor);
    %         tab_refined = concatenateTabRhohat(tab_refined,newTab);
    %         [costElimTable] = createCostElimTable(tab_refined,tabcn,vecAlphas);
    %         vecC = tab_refined.c;
    %         figure(1)
    %         plot_4dis(vecC,tab_refined)
    %     end
    
    %% On elargit l'intervalle
    %on ajoute les valeurs de cnn si elles ne sont pas deja presentes
    newC = tabcn.one;
    [newTab,~,~,~] = findRhohat_v5(N,paramTab,mu,b,newC,biasFactor);
    tab_refined = concatenateTabRhohat(tab_refined,newTab);
    [costElimTable,decreaseC_reason2] = createCostElimTable(tab_refined,tabcn,vecAlphas);
    vecC = tab_refined.c;
    %figure(3)
    %plot_4dis(vecC,tab_refined)
    
    increaseC = abs(max(max(table2array(costElimTable(:,2:5))))-vecC(end))<precision; t=0;
    while(increaseC && t<1000)
        disp('increasing boundary on the right')
        newC = vecC(end)+(max(vecC(2:end) - vecC(1:(end-1))));
        [newTab,~,~,~] = findRhohat_v5(N,paramTab,mu,b,newC,biasFactor);
        tab_refined = concatenateTabRhohat(tab_refined,newTab);
        [costElimTable,decreaseC_reason2] = createCostElimTable(tab_refined,tabcn,vecAlphas);
        vecC = tab_refined.c;
        increaseC = abs(max(max(table2array(costElimTable(:,2:5))))-vecC(end))<precision;
        t=t+1;
        %figure(2)
        %plot_4dis(vecC,tab_refined)
    end
    t=0;
    decreaseC_reason1 = abs(min(min(table2array(costElimTable(:,2:5))))-vecC(1))<precision;
    while((decreaseC_reason1 || decreaseC_reason2) && t<1000)
        disp('decreasing boundary on the left')
        newC = vecC(1)-(max(vecC(2:end) - vecC(1:(end-1))));
        [newTab,~,~,~] = findRhohat_v5(N,paramTab,mu,b,newC,biasFactor);
        tab_refined = concatenateTabRhohat(tab_refined,newTab);
        [costElimTable,decreaseC_reason2] = createCostElimTable(tab_refined,tabcn,vecAlphas);
        vecC = tab_refined.c;
        decreaseC_reason1 = abs(min(min(table2array(costElimTable(:,2:5))))-vecC(1))<precision;
        t=t+1;
    end
    
    %[alphaSorted, elimOrder] = sort(vecAlphas,'ascend');
    %elimOrder : diseases elimination order when c is decreasing
    
    %% Precision autour des valeurs de c qui peuvent poser probleme
    
    bestStrategyHasNotBeenFound = 1;
    while (bestStrategyHasNotBeenFound)
        [~,~,costStratTable_sorted,~] = createCostElimStrat(costElimTable);
        
        %si on ne peut pas distinguer les 2 premières strategies, alors
        % 1. problème d'échantillonnage
        % 2. égalité dans les modèles
        % 3. hasard
        k=0; err = optFindRhohat.Tol+1; iterMax=10; HIVEqu=0;
        while (costStratTable_sorted{1,'HIV'}==costStratTable_sorted{2,'HIV'} && k<iterMax && err>optFindRhohat.Tol)
            if 0
                %here : if model is a submodel of the other          
            else
                disp('affinage HIV')
                cost = costStratTable_sorted{1,'HIV'};
                cSuivant = tab_refined.c(find(abs(tab_refined.c-cost)<precision)+1);
                newC = (cost+cSuivant)/2;
                [newTab,~,~,~] = findRhohat_v5(N,paramTab,mu,b,newC,biasFactor);
                tab_refined = concatenateTabRhohat(tab_refined,newTab);
                [costElimTable] = createCostElimTable(tab_refined,tabcn,vecAlphas);
                [tempBest,tempWorst,costStratTable_sorted,costStratTable] = createCostElimStrat(costElimTable);
                k=k+1;
                %tempBest
                err = (cSuivant-cost);
            end
            HIVEqu = (k==iterMax); %we consider that the two models are equals
        end
        
        k=0;err=optFindRhohat.Tol+1;SyphEqu=0;
        while (HIVEqu && costStratTable_sorted{1,'Syphilis'}==costStratTable_sorted{2,'Syphilis'} && k<iterMax && err>optFindRhohat.Tol)
            if 0
                %here : if model is a submodel of the other
            else
                disp('affinage syph.')
                cost = costStratTable_sorted{1,'Syphilis'};
                cSuivant = tab_refined.c(find(abs(tab_refined.c-cost)<precision)+1);
                newC = (cost+cSuivant)/2;
                [newTab,~,~,~] = findRhohat_v5(N,paramTab,mu,b,newC,biasFactor);
                tab_refined = concatenateTabRhohat(tab_refined,newTab);
                [costElimTable] = createCostElimTable(tab_refined,tabcn,vecAlphas);
                [tempBest,tempWorst,costStratTable_sorted,costStratTable] = createCostElimStrat(costElimTable);
                
                k=k+1;
                %tempBest
                err = (cSuivant-cost);
            end
            SyphEqu = (k==iterMax);
        end
        
        %tri selon Ct
        k=0;err=optFindRhohat.Tol+1; CtEqu=0;
        while (HIVEqu && SyphEqu && costStratTable_sorted{1,'Chlam'}==costStratTable_sorted{2,'Chlam'} && k<iterMax && err>optFindRhohat.Tol)
            if 0
                %here : if model is a submodel of the other
            else
                disp('affinage Ct')
                cost = costStratTable_sorted{1,'Chlam'};
                cSuivant = tab_refined.c(find(abs(tab_refined.c-cost)<precision)+1);
                newC = (cost+cSuivant)/2;
                [newTab,~,~,~] = findRhohat_v5(N,paramTab,mu,b,newC,biasFactor);
                tab_refined = concatenateTabRhohat(tab_refined,newTab);
                [costElimTable] = createCostElimTable(tab_refined,tabcn,vecAlphas);
                [tempBest,tempWorst,costStratTable_sorted,costStratTable] = createCostElimStrat(costElimTable);
                k=k+1;
                %tempBest
                err = (cSuivant-cost);
            end
            CtEqu = (k==iterMax);
        end
        
        %tri selon Ng
        k=0;err=optFindRhohat.Tol+1;
        while (HIVEqu && SyphEqu && CtEqu && costStratTable_sorted{1,'Gono'}==costStratTable_sorted{2,'Gono'} && k<iterMax && err>optFindRhohat.Tol)
            if 0
                %here : if model is a submodel of the other
            else
                disp('affinage Ng.')
                cost = costStratTable_sorted{1,'Gono'};
                cSuivant = tab_refined.c(find(abs(tab_refined.c-cost)<precision)+1);
                newC = (cost+cSuivant)/2;
                [newTab,~,~,~] = findRhohat_v5(N,paramTab,mu,b,newC,biasFactor);
                tab_refined = concatenateTabRhohat(tab_refined,newTab);
                [costElimTable] = createCostElimTable(tab_refined,tabcn,vecAlphas);
                [tempBest,tempWorst,costStratTable_sorted,costStratTable] = createCostElimStrat(costElimTable);
                k=k+1;
                %tempBest
                err = (cSuivant-cost);
            end
            totalEqu = (k==iterMax); %egalite de 2 modeles
        end
        bestStrategyHasNotBeenFound = 0;
    end %end while
    bestStrategy = tempBest;
    worstStrategy = tempWorst;
    cleft = tab_refined.c(1);
    cright = tab_refined.c(end);
else
    bestStrategy='NA';
    worstStrategy='NA';
    costStratTable='NA';
end

end



