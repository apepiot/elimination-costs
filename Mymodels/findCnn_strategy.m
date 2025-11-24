function [elimCosts] = findCnn_strategy(paramTab,paramRho,mu,b,f,strat,aprioriBndsC,vecAlphas1d,afficherOutput)
nbModels = length(strat);

elimCosts = table('Size',[1,5],'VariableTypes',{'double','double','double','double','string'},'VariableNames',{'HIV','syphilis','Ct','Ng','code_err'});
elimCosts(1,1:5) = array2table(NaN(1,5));

%elimCosts(1,1:5) = array2table(rand(1,5));
%elimCosts(:,:).code_err = ['test'];

for m=1:nbModels %e.g. m=1->3 %strat{11} has 2 models: [1,2,4] and [3]
    %N = length(strat{m}); %nombre de maladies dans le modele e.g. N=3 puis 1
    %fprintf(['debut mod [',num2str(strat{m}),'] \n'])
    kit = strat{m}; k=[];
    
    aprioriBndsC.inf = -1;  aprioriBndsC.sup = 1;
    [alphaMax,~,msg_alpha] = findAlpha(paramTab,paramRho,mu,b,kit,afficherOutput); %il faudrait le calculer avant, dans le code RHOHAT_v7
    [cnn, ~, msg_kit] = find_Cnn_kit(paramTab,mu,b,paramRho,f,kit,aprioriBndsC,alphaMax,afficherOutput);
    
    if ismember('HIV',kit)
        elimCosts(:,:).HIV = cnn.HIV;
        k = [k,'h'];
    end
    if ismember('syphilis',kit)
        elimCosts(:,:).syphilis = cnn.syphilis;
        k = [k,'s'];
    end
    if ismember('Ct',kit)
        elimCosts(:,:).Ct = cnn.Ct;
        k = [k,'c'];
    end
    if ismember('Ng',kit)
        elimCosts(:,:).Ng = cnn.Ng;
        k = [k,'g'];
    end
    
    
    %Error msg
    msg_strat = ['a_',k,':',msg_alpha, '-c_',k,':',msg_kit];
    
end
elimCosts(:,:).code_err = msg_strat;
end


