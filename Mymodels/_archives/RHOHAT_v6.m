% New version : here we sample over each strategy %%%%PC%%%%%
% i.e. 1500 new set of parameters for each strategy
% strategy are independent
% we look for the cost of elimination of each disease 
% difference with RHOHAT_v5: new model with sictp, but forget to include 
% mandatory routine testing rate of stis under prep

clear all;
%% README
%---------------------------------------------------------------%
nbSimPerRound = 10;%1500
nbRounds = 3;
vecPHIV=[0.5]; %ptage of people on PrEP
vecF = [1];
vecStratNos=1:15;
backUpFolder='Runtest'; %mettre un nom de dossier qui n'existe pas
%---------------------------------------------------------------%
path2 = ['C:\Users\Moi\Documents\IPLESP\These\Codes\Simulations_StrategiesPC\',backUpFolder,'\'];

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
b=2000; %input parameters

[status, msg, msgID] = mkdir([path2]);
if isempty(msg)
    for roundNo=1:nbRounds
        recapSimID = NaN(length(vecF)*length(vecStratNos)*nbSimPerRound,6);
        currentTotSimSinceRun=0;
        for f=vecF
            for pHIV=vecPHIV
                elimCosts = NaN(nbSimPerRound,5,15,'double');
                tStart = tic;
                t1=tic;
                for i=vecStratNos %strategies
                    strat = strategies{i};
                    nbModels = length(strat); %nombre de modeles dans la strategies e.g. 3
                    %tabParametersByStrat
                    for sim=1:nbSimPerRound %1->1500
                        unIDsim = str2double([char(datetime('now','TimeZone','local','Format','yyMMddHHmmss')),num2str(sim)]); %unique ID based on the current time 
                        currentTotSimSinceRun = currentTotSimSinceRun+1;
                        recapSimID(currentTotSimSinceRun,1:6) = [unIDsim,roundNo,i,f,sim,pHIV]; %common to Ct,Ng,HIV, syph

                        [paramTab,mu,vecAlphas] = sampleParameters_v3(true,true,true,true,b);   %Ct,Ng,HIV,syph  
                        
                        % R0= Rp_SICTP(...,rho=0) or rho=rhob? shouldn't be lower than 1
                        paramTab{3}.modelType='SICTP';
                        paramTab{3}.p = pHIV;
                        paramTab{3}.eta = 4;
                        paramTab{3}.mu = mu;
                        paramTab{3}.zeta=randPERT(46,60,71,1)/100;
                        paramTab{3}.alpha0 = paramTab{3}.alpha;
                        [paramTab{3}.RSICTP,~,paramTab{3}.alpha] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
                            paramTab{3}.theta,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
                            paramTab{3}.p,paramTab{3}.mu,b,paramTab{3}.rhob);
                        vecAlphas(3) = paramTab{3}.alpha;
                        
                        %Store parameters in a file
                        allParametersTabCt(currentTotSimSinceRun,:) = [struct2table(paramTab{1}),table(unIDsim)];
                        allParametersTabNg(currentTotSimSinceRun,:) = [struct2table(paramTab{2}),table(unIDsim)];
                        allParametersTabHIV(currentTotSimSinceRun,:) = [struct2table(paramTab{3}),table(unIDsim)];
                        allParametersTabS(currentTotSimSinceRun,:) = [struct2table(paramTab{4}),table(unIDsim)];

                        if paramTab{3}.RSICTP<=1
                            %add the parameter values in a table
                            elimCosts(sim,1:5,i) = [NaN,NaN,NaN,NaN,unIDsim];
                        else
                            for m=1:nbModels %e.g. m=1->3 %strat{11} has 2 models: [1,2,4] and [3]
                                N = length(strat{m}); %nombre de maladies dans le modele e.g. N=3 puis 1   
                                fprintf(['debut mod [',num2str(strat{m}),'] \n'])
                                if N==1 %si c'est un modele a une maladie, on connait la valeur exacte de cnn
                                    n = strat{m}; %num of the disease
                                    [~,~,celim_n] = findRhohat1d_v5(paramTab{n}.modelType,paramTab{n},mu,b,0,f);
                                    elimCosts(sim,n,i) = celim_n;                  
                                    fprintf(['strat:',num2str(i), ' sim:',num2str(sim), ' mod:[', num2str(strat{m}),'] c_',num2str(n),'=', num2str(round(celim_n,3)),'\n'])
                                else 
                                    vecC = linspace(-0.5,0.02,4);  %a priori %% TO DO : reajuster l'a priori avec les valeurs trouvees precedemment
                                    %encadrer avec les valeurs c max et min obtenus
                                    %avec le single disease model
                                    [tabRhohat,tabC0,tabCnn,tabTimes,vecC] = findRhohat_v6(N,paramTab(strat{m}),mu,b,vecC,f);
                                    tab_refined = tabRhohat;
                                    t=0;
                                    while(tab_refined.rhohat.rhohat(1)<max(vecAlphas(strat{m})) && t<100) %si on ne calcule pas suffisamment à gauche
                                    %while(tab_refined.rhohat.rhohat(1)<max(vecAlphas) && t<100)    
                                        disp('on augmente vecC a gauche')
                                        newC = vecC(1)-(max(vecC(2:end) - vecC(1:(end-1))));
                                        [newTab,~,~,~] = findRhohat_v6(N,paramTab(strat{m}),mu,b,newC,f);
                                        tab_refined = concatenateTabRhohat(tab_refined,newTab);
                                        vecC = tab_refined.c;
                                        t=t+1;
                                    end
                                    t=0;
                                    while (tab_refined.rhohat.rhohat(end)>min(vecAlphas(strat{m})) && t<100)
                                    %while (tab_refined.rhohat.rhohat(end)>min(vecAlphas) && t<100)
                                        disp('on augmente vecC a droite')
                                        newC = vecC(end)+(max(vecC(2:end) - vecC(1:(end-1))));
                                        [newTab,~,~,~] = findRhohat_v6(N,paramTab(strat{m}),mu,b,newC,f);
                                        tab_refined = concatenateTabRhohat(tab_refined,newTab);
                                        vecC = tab_refined.c;
                                        t=t+1;
                                    end
                                    %on cherche le cout d'elimination de chaque maladie en affinant par dichotomie
                                    for n=strat{m} %e.g. n=1 puis n=2 puis n=4
                                        alphan = paramTab{n}.alpha;
                                        affiner=true;
                                        while affiner && N>1
                                            disp('on affine')
                                            temp1 = vecC(tab_refined.rhohat.rhohat>=alphan);
                                            borneinf = temp1(end);
                                            temp2 = vecC(tab_refined.rhohat.rhohat<alphan);
                                            bornesup = temp2(1);
                                            precision = bornesup-borneinf;
                                            t=0;
                                            if (precision>=1e-3)
                                                t=t+1;
                                                newC = (bornesup+borneinf)/2;
                                                [newTab,~,~,~] = findRhohat_v6(N,paramTab(strat{m}),mu,b,newC,f);
                                                tab_refined = concatenateTabRhohat(tab_refined,newTab);
                                                vecC = tab_refined.c;
                                                t=t+1;
                                            else
                                                fprintf(['precision=',num2str(precision),'\n'])
                                                affiner=false;
                                                celim_n = newC; 
                                            end
                                        end
                                        elimCosts(sim,n,i) = celim_n;
                                        fprintf(['strat:',num2str(i), ' sim:',num2str(sim), ' mod:[', num2str(strat{m}),'] c_',num2str(n),'=', num2str(round(celim_n,3)), '\n'])
                                    end                            
                                end
                                elimCosts(sim,5,i) = unIDsim;
                                fprintf(['strat:',num2str(i), ' sim:',num2str(sim), ' mod:[', num2str(strat{m}),'] time:',num2str(round(toc(t1),1)),'s','\n'])
                                t1=tic;
                            end 
                        end
                    end
                    tabElim = array2table(elimCosts(:,:,i),"VariableNames",["Ct","Ng","HIV","syph","IDsim"]);
                    [status, msg, msgID] = mkdir([path2,'_round_',num2str(roundNo)]);
                    %if isempty(msg) %folder round_i does not exist and has been created
                        fileName = ['_round_',num2str(roundNo),'\strat_',num2str(i),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'_',num2str("1"),'-',num2str(nbSimPerRound),'_elim_cost.txt'];
                        writetable(tabElim,[path2,fileName],'WriteVariableNames',true)
                    %else
                    %end     
                    writetable(array2table(recapSimID,"VariableNames",["IDsim","roundNo","stratNo","b","simNo","pHIV"]),...
                        [path2,'_round_',num2str(roundNo),'\recapSimID'],'WriteVariableNames',true)   

                    %Parameters backup
                    writetable(allParametersTabCt,[path2,'_round_',num2str(roundNo),'\allParametersTabCt'],'WriteVariableNames',true)
                    writetable(allParametersTabNg,[path2,'_round_',num2str(roundNo),'\allParametersTabNg'],'WriteVariableNames',true)
                    writetable(allParametersTabHIV,[path2,'_round_',num2str(roundNo),'\allParametersTabHIV'],'WriteVariableNames',true)
                    writetable(allParametersTabS,[path2,'_round_',num2str(roundNo),'\allParametersTabS'],'WriteVariableNames',true)
                end
            end
        end
    end
    tEnd=toc(tStart);
else
    error([backUpFolder, ' already exists. Change the name of the back-up folder'])
end


