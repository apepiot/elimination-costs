function [] = createStatTests_v2(f,pHIV,path,test)
[status, msg, msgID] = mkdir([path,'\_tests']);

if test=="Student"
    myTest = @myStudentTest;
elseif test == "StudentMatlab"
    myTest = @ttest2;
elseif test == "MWW"
    myTest = @myMWWTest;
elseif test == "KS" %Kolmogorov-Smirnov
    myTest = @kstest2;
else
    error('the name of the test is unknown')
end

recapAllTest = zeros(15,15,4); %strat vs strat vs infection %p_values
for stratA = 1:15
    
    stratAfile = readtable([path,'\strat_',num2str(stratA),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),...
                '_elim_cost','_concatenated.txt']);
    recapTest = table('Size',[14,10],'VariableTypes',["string",repmat("double",1,9)],...
        'VariableNames',[test,"z_Ct","p_Ct","z_Ng","p_Ng","z_HIV","p_HIV","z_s","p_s","N_i"]);
    
    for stratB=1:15
        tot = readtable([path,'\strat_',num2str(stratB),'_b_',num2str(f),'_p_',num2str(round(pHIV*100)),...
                '_elim_cost','_concatenated.txt']);
        recapTest(stratB,test) = {['strat_',num2str(stratA),'_vs_strat_',num2str(stratB)]};
        [z1,p1] = myTest(stratAfile.Ct,tot.Ct);
        recapTest(stratB,"z_Ct") = {z1};
        recapTest(stratB,"p_Ct") = {p1};
        recapAllTest(stratA,stratB,1) = p1;
        
        %[H2,P2,CI2,STATS2] = ttest2(strat1.Ng,tot.Ng,'Vartype','unequal');
        [z2,p2] = myTest(stratAfile.Ng,tot.Ng);
        recapTest(stratB,"z_Ng") = {z2};
        recapTest(stratB,"p_Ng") = {p2};
        recapAllTest(stratA,stratB,2) = p2;
        
        %[H3,P3,CI3,STATS3] = ttest2(strat1.HIV,tot.HIV,'Vartype','unequal');
        [z3,p3] = myTest(stratAfile.HIV,tot.HIV);
        recapTest(stratB,"z_HIV") = {z3};
        recapTest(stratB,"p_HIV") = {p3};
        recapAllTest(stratA,stratB,3) = p3;
        
        %[H4,P4,CI4,STATS4] = ttest2(strat1.syph,tot.syph,'Vartype','unequal');
        [z4,p4] = myTest(stratAfile.syph,tot.syph);
        recapTest(stratB,"z_s") = {z4};
        recapTest(stratB,"p_s") = {p4};
        recapTest(stratB,"N_i") = {length(tot.syph)};
        recapAllTest(stratA,stratB,4) = p4;
        
        recapTest{:,2:10} = round(recapTest{:,2:10},3);

    end
    % diff significative si p<0.05
    writetable(recapTest,[path,'_tests\_',char(test),'_versus_strategy_',num2str(stratA),'.txt'])
end
recapAllTest = round(recapAllTest,3);

%Ct
recapAllTestCt = array2table(recapAllTest(:,:,1),'VariableNames',"strat_" + (1:15));
recapAllTestCt{:,"Strat"} = ["strat_" + (1:15)]';
writetable(recapAllTestCt, [path,'_tests\_',char(test),'_p_values_Ct_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt'])

%Ng
recapAllTestNg = array2table(recapAllTest(:,:,2),'VariableNames',"strat_" + (1:15));
recapAllTestNg{:,"Strat"} = ["strat_" + (1:15)]';
writetable(recapAllTestNg, [path,'_tests\_',char(test),'_p_values_Ng_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt'])

%HIV
recapAllTestHIV = array2table(recapAllTest(:,:,3),'VariableNames',"strat_" + (1:15));
recapAllTestHIV{:,"Strat"} = ["strat_" + (1:15)]';
writetable(recapAllTestHIV, [path,'_tests\_',char(test),'_p_values_HIV_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt'])

%S
recapAllTestS = array2table(recapAllTest(:,:,4),'VariableNames',"strat_" + (1:15));
recapAllTestS{:,"Strat"} = ["strat_" + (1:15)]';
writetable(recapAllTestS, [path,'_tests\_',char(test),'_p_values_S_b_',num2str(f),'_p_',num2str(round(pHIV*100)),'.txt'])
end

