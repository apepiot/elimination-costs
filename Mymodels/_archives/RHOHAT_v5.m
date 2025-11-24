% New version : here we sample over each strategy %%%%PC%%%%%
% i.e. 1500 new set of parameters for each strategy
% strategy are independent
% we look for the cost of elimination of each disease 
clear all;
%% README
%attention en relançant a ne pas ecraser les fichiers
nbSim = 300;%1500
path2 = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\';
%path2 = '_resultats/';

%%
%Liste des strategies
%1:Ct, 2:Ng, 3:HIV, 4:syph
strategies{1} = {1,2,3,4};
strategies{2} = {[1,2],3,4};
strategies{3} = {[1,3],2,4};
strategies{4} = {[1,4],2,3};
strategies{5} = {1,[2,3],4};
strategies{6} = {1,[2,4],3};
strategies{7} = {1,2,[3,4]};
strategies{8} = {[1,2],[3,4]};
strategies{9} = {[1,3],[2,4]};
strategies{10} = {[1,4],[2,3]};
strategies{11} = {[1,2,3],4};
strategies{12} = {[1,2,4],3};
strategies{13} = {[1,3,4],2};
strategies{14} = {1,[2,3,4]};
strategies{15} = {[1,2,3,4]};
b=2;
for f=1
    elimCosts = zeros(nbSim,4,15);
    tStart = tic;
    t1=tic;
    for i=[3,9,13]
        strat = strategies{i}; %e.g. {[1,2,4],[3]}
        nbModels = length(strat); %nombre de modeles dans la strategies e.g. 3
        for sim=1:nbSim %1->1500
            [paramTab,mu,vecAlphas] = sampleParameters(true,true,true,true,b);        
            for m=1:nbModels %e.g. m=1->3
                %strat{m} %e.g. [1,2,4] puis [3]
                N = length(strat{m}); %nombre de maladies dans le modele e.g. N=3 puis 1   
                fprintf(['debut mod [',num2str(strat{m}),'] \n'])
                if N==1 %si c'est un modele a une maladie, on connait la valeur exacte de cnn
                    n = strat{m}; %num of the disease
                    [~,~,celim_n] = findRhohat1d_v4(paramTab{n}.modelType,paramTab{n},mu,b,0,f);
                    elimCosts(sim,n,i) = celim_n;                  
                    fprintf(['strat:',num2str(i), ' sim:',num2str(sim), ' mod:[', num2str(strat{m}),'] c_',num2str(n),'=', num2str(round(celim_n,3)),'\n'])
                else 
                    vecC = linspace(-0.1,0.0,5); % apriori (TO DO : reajuster l'a priori avec les valeurs trouvees precedemment)
                    %strat{m} % e.g. [1,2,4], N=3;
                    [tabRhohat,tabC0,tabCnn,tabTimes,vecC] = findRhohat_v5(N,paramTab(strat{m}),mu,b,vecC,f);
                    tab_refined = tabRhohat;
                    t=0;
                    while(tab_refined.rhohat.rhohat(1)<max(vecAlphas) && t<100) %si on ne calcule pas suffisamment à gauche
                        %disp('on augmente a gauche')
                        newC = vecC(1)-(max(vecC(2:end) - vecC(1:(end-1))));
                        [newTab,~,~,~] = findRhohat_v5(N,paramTab(strat{m}),mu,b,newC,f);
                        tab_refined = concatenateTabRhohat(tab_refined,newTab);
                        vecC = tab_refined.c;
                        t=t+1;
                    end
                    t=0;
                    while (tab_refined.rhohat.rhohat(end)>min(vecAlphas) && t<100)
                        %disp('on augmente a droite')
                        newC = vecC(end)+(max(vecC(2:end) - vecC(1:(end-1))));
                        [newTab,~,~,~] = findRhohat_v5(N,paramTab(strat{m}),mu,b,newC,f);
                        tab_refined = concatenateTabRhohat(tab_refined,newTab);
                        vecC = tab_refined.c;
                        t=t+1;
                    end
                    %on cherche le cout d'elimination de chaque maladie et on affine si besoin
                    for n=strat{m} %e.g. n=1 puis n=2 puis n=4
                        alphan = paramTab{n}.alpha;
                        affiner=true;
                        while affiner && N>1 %if N
                            %disp('on affine')
                            temp1 = vecC(tab_refined.rhohat.rhohat>=alphan);
                            borneinf = temp1(end);
                            temp2 = vecC(tab_refined.rhohat.rhohat<alphan);
                            bornesup = temp2(1);
                            precision = bornesup-borneinf;
                            t=0;
                            if (precision>=1e-3)
                                t=t+1;
                                newC = (bornesup+borneinf)/2;
                                [newTab,~,~,~] = findRhohat_v5(N,paramTab(strat{m}),mu,b,newC,f);
                                tab_refined = concatenateTabRhohat(tab_refined,newTab);
                                vecC = tab_refined.c;
                                t=t+1;
                            else
                                %fprintf(['precision=',num2str(precision),'\n'])
                                affiner=false;
                                celim_n = newC; 
                            end
                        end
                        elimCosts(sim,n,i) = celim_n;
                        fprintf(['strat:',num2str(i), ' sim:',num2str(sim), ' mod:[', num2str(strat{m}),'] c_',num2str(n),'=', num2str(round(celim_n,3)), '\n'])
                    end
                end
                fprintf(['strat:',num2str(i), ' sim:',num2str(sim), ' mod:[', num2str(strat{m}),'] time:',num2str(round(toc(t1),1)),'s','\n'])
                t1=tic;
            end        
        end
        fileName = ['Strat_',num2str(i),'_b_',num2str(f),'_tot_',num2str(nbSim),'.txt'];
        tabElim = array2table(elimCosts(:,:,i),"VariableNames",["1","2","3","4"]);
        writetable(tabElim,[path2,fileName],'WriteVariableNames',true)
    end
end
tEnd=toc(tStart);