clear all;
%Création de jeux de paramètres
pathBackup = 'C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/ParameterAnalysis/paramSets_1';
mkdir(pathBackup)
b=100;
for roundNo=1:100
    simRndNo=0;
    pathRound = [pathBackup,'/round_',num2str(roundNo)];
    [err,msg] = mkdir(pathRound);
    if isempty(msg)
        tic
        for nbEch=1:10
            [paramTab,~,~] = sampleParameters_v3_extent(1,1,1,1,b,0);
            IDset_id = [char(datetime('now','TimeZone','local','Format','MMddHHmmss')),num2str(nbEch,'%04.f');];  %unique ID based on the current time
            IDech_id = str2double(IDset_id);
            paramTab{1} = orderfields(paramTab{1}); paramTab{2} = orderfields(paramTab{2}); paramTab{4} = orderfields(paramTab{4});
            paramSetsCt(nbEch,:) = [struct2table(paramTab{1}),table(IDech_id,nbEch,roundNo)];
            paramSetsNg(nbEch,:) = [struct2table(paramTab{2}),table(IDech_id,nbEch,roundNo)];
            paramSetsS(nbEch,:)  = [struct2table(paramTab{4}),table(IDech_id,nbEch,roundNo)];
            for pHIV=0:0.1:1
                simRndNo=simRndNo+1;
                %-- Update pHIV and the set of parameters
                paramTab{3}.p = pHIV;
                [paramTab{3}.R_prep_base,~,~,paramTab{3}.Ptot_prep_base,paramTab{3}.Pun_prep_base] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
                    paramTab{3}.theta0,paramTab{3}.gamma0,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
                    paramTab{3}.p,paramTab{3}.mu,b,paramTab{3}.rhob);
                [paramTab{3}.R_prep_0,~,paramTab{3}.alpha_prep,paramTab{3}.Ptot_prep_0,paramTab{3}.Pun_prep_0] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
                    paramTab{3}.theta0,paramTab{3}.gamma0,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
                    paramTab{3}.p,paramTab{3}.mu,b,0);
                paramTab{3} = orderfields(paramTab{3});
                %-- end of update pHIV and the set of parameters
                
                paramSetsHIV(simRndNo,:)= [struct2table(paramTab{3}),table(IDech_id,simRndNo,nbEch,roundNo)];
                
            end
        end
        
        writetable(paramSetsCt,[pathRound,'/allParametersSets_Ct'],'WriteVariableNames',true)
        writetable(paramSetsNg,[pathRound,'/allParametersSets_Ng'],'WriteVariableNames',true)
        writetable(paramSetsHIV,[pathRound,'/allParametersSets_HIV'],'WriteVariableNames',true)
        writetable(paramSetsS,[pathRound,'/allParametersSets_syphilis'],'WriteVariableNames',true)
        
        clear paramSetsCt
        clear paramSetsNg
        clear paramSetsHIV
        clear paramSetsS
        
        tpsRound = toc
        if (tpsRound<1)
            pause(5)
        end
    else
        disp(msg)
    end
end



