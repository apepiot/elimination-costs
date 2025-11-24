function [bestStrategy,worstStrategy,costStratTable,cnn,c0,costStratTable_sorted,tab,tabco,tabcn] = findBestStrategy(computeC0Cnn,optFindCnnAndC0,optFindRhohat,N,paramTab,mu,b,biasFactor,plotRhohat)
    
    vecAlphas=zeros(1,size(paramTab,2));
    for i=1:size(paramTab,2)
    vecAlphas(i) = paramTab{i}.alpha;
    end
    if(computeC0Cnn)       
        errMax = optFindCnnAndC0.errMax;
        Tmax   = optFindCnnAndC0.Tmax;
        aPrioriCnn = optFindCnnAndC0.aPrioriCnn;
        aPrioriC0 = optFindCnnAndC0.aPrioriC0;
        [c0,cnn,t0,t1] = findCnnAndC0(aPrioriCnn,aPrioriC0,Tmax,errMax,N,paramTab,mu,b,biasFactor);
        cleft  = cnn - (c0-cnn)*5/9 ;
        cright = c0 + (c0-cnn)*5/9; %c0 + abs(c0)*errMax;
        disp('c0 and cnn trouves')
    else
        %error('Define boundaries of vecC : c_min and c_max to compute argmax : code to implement')
        cleft = optFindCnnAndC0.aPrioriCnn;
        cright = optFindCnnAndC0.aPrioriC0;
        cnn = optFindCnnAndC0.aPrioriCnn;
        c0 = optFindCnnAndC0.aPrioriC0;
    end
    
    
    %% This section gives rhohat=argmaxU for a given model
    vecC = [cleft,linspace(cnn-(c0-cnn)/3,c0+(c0-cnn)/3,optFindRhohat.sampleSize),cright];
    tic;
    [tab,tabco,tabcn,tabTimes] = findRhohat_v4(N,paramTab,mu,b,vecC,biasFactor);
    tps = toc;

    %% cswitch for the n disease model (swith from argmax Uijk to Uij) and PLOT 
    if (plotRhohat)
            nSEIIS=0;nSICR=0;nSEIIIS=0;
            for i=1:size(paramTab,2)
                nSEIIS = nSEIIS + strcmp(paramTab{i}.modelType,'SEIIS');
                nSICR = nSICR + strcmp(paramTab{i}.modelType,'SICR');
                nSEIIIS = nSEIIIS + strcmp(paramTab{i}.modelType,'SEIIIS');
            end
        if N==1
            disp('plot : a verifier, notamment les arguments dans les fonctions pour trouver les cs')
            if nSICR==1
                [~,cs.c0] = U1_SICR_v4(paramTab{1},mu,b,0,0,1); %c2 = c0 
                [~,cs.cs1dis] = U1_SICR_v4(paramTab{1},mu,b,paramTab{1}.alpha,0,1); %c1=c1
            end
            if nSEIIS==1 
                [~,cs.c0] = U1_SEIISv4(paramTab{1},mu,b,0,0,1); %c2 = c0
                [~,cs.cs1dis] = U1_SEIISv4(paramTab{1},mu,b,paramTab{1}.alpha,0,1);
            end
            if nSEIIIS==1
                [~,cs.c0] = U1_SEIIIS_v4(paramTab{1},mu,b,0,0,1);
                [~,cs.cs1dis] = U1_SEIIIS_v4(paramTab{1},mu,b,paramTab{1}.alpha,0,1); %c1=c1
            end
        elseif N>=2
            cs = findThresholds_v4(nSEIIS,nSICR,nSEIIIS, tab, vecAlphas, vecC);
        end

        if (N>=2)
            ics1 = find(cs.cs1dis==vecC);
            ics2 = find(cs.cs2dis==vecC);
            if (N>=3)
                ics3 = find(cs.cs3dis==vecC);
                if (N>=4)
                    ics4 = find(cs.cs4dis==vecC);
                end
            end
        end
        plot_zones_v3bis;
    end


    %% thresholds of cost for disease elimination
    if N==4
        %rowNames = split(num2str(tab.single),' ');
        %costElimTable = zeros(4^2-1,5);
        costElimTable = table('Size',[4^2-1,5], 'VariableNames',{'model', 'Chlam', 'Gono', 'HIV', 'Syphilis'},...
            'VariableTypes', {'string','double','double','double','double'});

        costElimTable = standardizeMissing(costElimTable,0);
        %1 disease model
        for i=1:4
            costElimTable{i,1} = {num2str(i)};
            costElimTable{i,i+1} = round(tabcn.one(i),3,"significant");
        end

        %2 disease model
        rhohat2d_reshape = reshape([tab.two(:).rhohat],[],6);
        for k=1:6
            nameModel = num2str(tab.duos(k,:));
            costElimTable{k+4,1} = {nameModel(find(~isspace(nameModel)))};

            for dis=tab.duos(k,:)
                j=max(find(rhohat2d_reshape(:,k)>=vecAlphas(dis)));
                if ~isempty(j)
                    costElimTable{k+4,dis+1} = round(vecC(j),3,"significant");
                end
            end
        end

        %3 disease model
        rhohat3d_reshape = reshape([tab.three(:).rhohat],[],4);

        for k=1:4
            nameModel = num2str(tab.trios(k,:));
            costElimTable{k+4+6,1} = {nameModel(find(~isspace(nameModel)))};

            for dis=tab.trios(k,:)
                j=max(find(rhohat3d_reshape(:,k)>=vecAlphas(dis)));
                if ~isempty(j)
                    costElimTable{k+4+6,dis+1} = round(vecC(j),3,"significant");
                end      
            end
        end

        %4 disease model
        nameModel = num2str('1234');
        costElimTable{4+6+4+1,1} = {nameModel};
        for dis=1:4
            j=max(find(tab.rhohat.rhohat>=vecAlphas(dis)));
            if ~isempty(j)
                costElimTable{1+4+6+4,dis+1} = round(vecC(j),3,"significant");
            end      
        end

        %% Odering strategies by cost
        %e.g. 1x2 + 3 + 4 ou 1x2 + 3x4

        costStratTable = table('Size',[15,5], 'VariableNames',{'strategies', 'Chlam', 'Gono', 'HIV', 'Syphilis'},...
        'VariableTypes', {'string','double','double','double','double'});
        costStratTable.strategies = {{'1','2','3','4'}, {'1', '2', '34'}, {'1', '234'},...
                     {'1', '23', '4'}, {'1', '24', '3'},...
                     {'12', '3', '4'}, {'12', '34'},...
                     {'13', '2', '4'}, {'13', '24'},...
                     {'14', '2', '3'}, {'14', '23'},...
                     {'123','4'},{'124','3'},{'134','2'},...
                     {'1234'}}';
        strategies = costStratTable.strategies;

         %costs associated to each strategy
         nbStrat = size(strategies,1);
         for i=1:nbStrat
             nbKits = size(strategies{i},2);
             for j=1:nbKits
                kitCost = costElimTable{costElimTable.model == strategies{i}{j},2:5};
                costStratTable{i,2:5} = sum([kitCost; costStratTable{i,2:5}],1,'omitnan') ;
             end
         end

         %detailing strategies
          for i=1:nbStrat
            costStratTable.strategiesdetailed{i} = strjoin(costStratTable.strategies{i},'+');
            %costStratTable.strategiesdetailed{i};
          end

        %% sorting by costs (HIV, syphilis, chlamydia, gonorrhea)
        costStratTable_sorted = sortrows(costStratTable,[{'HIV'},{'Syphilis'},{'Chlam'},{'Gono'}],'descend');
        costStratTable_sorted = costStratTable_sorted(:,[{'strategiesdetailed','HIV'},{'Syphilis'},{'Chlam'},{'Gono'}]);

        %writetable(costStratTable_sorted,'myData.txt','Delimiter',' ')  

        bestStrategy = costStratTable_sorted{1,1};
        
        if(sum(table2array(costStratTable_sorted(1,2:5))==table2array(costStratTable_sorted(2,2:5)))==4)
            bestStrategy = {'equality '};
        end
        
        worstStrategy = costStratTable_sorted{nbStrat,1};
        
    else
       bestStrategy='NA';
       worstStrategy='NA';
       costStratTable='NA';
    end
 
    tab.c = vecC;
end



