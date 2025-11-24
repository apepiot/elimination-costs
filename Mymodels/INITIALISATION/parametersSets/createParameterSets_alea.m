clear all;
%Création de jeux de paramètres
parametrizationNo = '11';
vecP=0.6:0.1:1; %0:0.1:0.5;
pathBackup = ['C:/Users/Moi/Documents/IPLESP/These/Codes/multi-voluntary-testing/Mymodels/ParameterAnalysis/paramSets_',parametrizationNo,'/'];
mkdir(pathBackup)
b=10000;
for roundNo=1:10 %1:160 
    nbEch=0;
    simRndNo=0;
    pathRound = [pathBackup,'round_',num2str(roundNo)];
    [err,msg] = mkdir(pathRound);
    if isempty(msg)
        tic
        for i=1:50  %simu par p par modele par dossier (id_ech jusqu'a 9999, 50(i)*6(p)=300 :ok)
            for pHIV=vecP
                nbEch=nbEch+1;
                [paramTab,~,~] = sampleParameters_v3_extent(1,1,1,1,b,0);
                IDset_id = [char(datetime('now','TimeZone','local','Format','MMddHHmmss')),num2str(nbEch,'%04.f');];  %unique ID based on the current time
                IDech_id = str2double(IDset_id);
                paramTab{1} = orderfields(paramTab{1}); paramTab{2} = orderfields(paramTab{2}); paramTab{4} = orderfields(paramTab{4});
                paramSetsCt(nbEch,:) = [struct2table(paramTab{1}),table(IDech_id,nbEch,roundNo)];
                paramSetsNg(nbEch,:) = [struct2table(paramTab{2}),table(IDech_id,nbEch,roundNo)];
                paramSetsS(nbEch,:)  = [struct2table(paramTab{4}),table(IDech_id,nbEch,roundNo)];

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

                %disp(IDech_id)
                pause(0)
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



