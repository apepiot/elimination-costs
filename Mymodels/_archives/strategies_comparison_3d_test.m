%% Strategy comparison
%% code for 2 or 3 diseases et for 1 ?
clear all;
s1=1;s2=1;s3=1;
mu = 1/35;
nSIS=2;nSIR=1; N=nSIS+nSIR;

len = 5; %sampling R(0)'s
R10 = linspace(1.1,5,len);
R20 = linspace(1.1,5,len);
R30 = linspace(1.1,5,len);

d1 = 1/12:1/12:6/12;
nbItertot = length(d1)^3*len^3;
elimCostsCVT = zeros(nbItertot,9); 
elimCostsTVT = zeros(nbItertot,9); iterCurrent = 0;
timesTable   = zeros(nbItertot,1);

for gamma1 = 1./d1
    if(nSIS==2)
        d2 = 1/12:1/12:1/gamma1; %the problem is symmetrical (we just need to run half the values)
    else
        d2 = d1;
    end
    for gamma2 = 1./d2
        for gamma3 = 1./d1
            
            disp(['gamma1=',num2str(gamma1),' gamma2=',num2str(gamma2),' gamma3=',num2str(gamma3)])
            %[R1,R2,R3] = meshgrid(R10,R20,R30);
            BETA1 = R10.*(gamma1+mu);
            BETA2 = R20.*(gamma2+mu);
            BETA3 = R30.*(gamma3+mu);

            ALPHA1 = BETA1./s1.*(1-1./R10);
            ALPHA2 = BETA2./s2.*(1-1./R20);
            ALPHA3 = BETA3./s3.*(1-1./R30);
            C11 = 1./R10-1;
            C22 = (1./R20-1)*(nSIS==2)+(-mu./BETA2.*(R20-1))*(nSIR==1 & nSIS==1);           
            C33 = (1./R30-1)*(nSIS==3)+(-mu./BETA3.*(R30-1))*(nSIR==1 & nSIS==2);

            CELI=zeros(4,len); CELJ=zeros(4,len); CELK=zeros(4,len); %store cs/the second/first disease is eliminated
            CEL1=zeros(4,len); CEL2=zeros(4,len); CEL3=zeros(4,len); %store cs/ disease 1 (resp. 2) is eliminated (R1,R2,R3,c...)
            for i=1:len
                for j=1:len
                    for k=1:len
                        iterCurrent = iterCurrent+1;
                        disp(['Iter ', num2str(iterCurrent), '/', num2str(nbItertot)])
                        tic
                        %disp(['i=',num2str(i),' j=',num2str(k),' k=',num2str(k)])
                        c11 = C11(i); c01 = -c11; %SIS
                        c22 = C22(j); c02 = -c22*(nSIS==2) + mu/BETA2(j)*(1-1/R20(j))*(nSIR==1 & nSIS==1); %SIS^2 ou SISxSIR
                        c33 = C33(k); c03 = -c33*(nSIS==3) + mu/BETA3(k)*(1-1/R30(k))*(nSIR==1 & nSIS==2); %SIS^2,SIR ou SIS^3
                        %j : first disease eliminated

                        mincnn = min([c11,c22,c33]); maxc0n = max([c01,c02,c03]);
                        maxcnn = max([c11,c22,c33]);

                        %interval of c
                        vecC = linspace(1.5*mincnn,maxc0n,50);
                        diff = 0;
                        step_init = vecC(2)-vecC(1);
                        vecAlpha = [ALPHA1(i),ALPHA2(j),ALPHA3(k)]; %alphaj = min(vecAlpha);

                        % we look for the cost of each disease elimination
                        paramSIS = [BETA1(i),gamma1,1;BETA2(j),gamma2,s2];
                        paramSIR = [BETA2(j),gamma2,s2]*(nSIS==1 & nSIR==1)+[BETA3(k),gamma3,s3]*(nSIS==2 & nSIR==1);

                        findRhohat_ij = @(c) findRhohat(nSIS,nSIR,0,paramSIS(1:nSIS,:),paramSIR(1:nSIR,:),[],mu,5,c);
                        findCs_ij = @(c,tabRho) findThresholds(nSIS,nSIR,0,tabRho,vecAlpha,[BETA1(i),BETA2(j),BETA3(k)],[gamma1,gamma2,gamma3],[s1,s2,s3],5,mu,c);
                        [tab,tabco,tabcn] = findRhohat_ij(vecC);
                        cs_init = findCs_ij(vecC,tab);
                        tab_init = tab; 
                        vecC_init = vecC;
                        elimOrder = cs_init.order;

                        for n = 1:N %n = disease number
                            disp(['nDis:',num2str(n)])
                            cnijnotfound=1; %cnij : cost at which CVT eliminates n (n in {i,j})
                            %below : to determine either cnn<cnij or cnn>cnij,
                            %cnij needs to be sufficiently accurate
                            cnn = c11*(n==1) + c22*(n==2)+c33*(n==3);
                            cs = cs_init; %thresholds found during the first run of findThresholds (reinitializing for each disease to avoid crushing of the following cnij)
                            tab = tab_init;
                            alphan = vecAlpha(n);
                            vecC = vecC_init; step = step_init;

                            while (cnijnotfound)
                                %disp('cnijnotfound')
                                if n==elimOrder(1) %if n is the first disease eliminated
                                    cnij = cs.cs1dis;
                                elseif n==elimOrder(2) %if n is the second eliminated
                                    cnij = cs.cs2dis;
                                elseif n==elimOrder(3)
                                    cnij = cs.cs3dis;
                                end
                                if(isempty(cnij)) %vecC does not include cnij
                                    %disp(['cnij empty, i=',num2str(i),', j=',num2str(j)]);
                                    if(tab.rhohat(1)<alphan) %vecC did not go far enough on the left side
                                        vecC = linspace(vecC(1)-10*step_init, vecC(1),10);
                                        [tab,tabco,tabcn] = findRhohat_ij(vecC);
                                        cs = findCs_ij(vecC,tab);
                                        %disp('if tab.rhohat(1)<alphan');
                                        continue; %next iteration of the while loop
                                    elseif (tab.rhohat(end)>=alphan)
                                        %disp('if tab.rhohat(end)>alphan')
                                        vecC = linspace(vecC(end),vecC(end)+10*step_init,10);
                                        [tab,tabco,tabcn] = findRhohat_ij(vecC);
                                        cs = findCs_ij(vecC,tab);
                                        continue; %next iteration of the while loop
                                    end
                                end

                                %if cnij has been found (not empty)
                                diff = abs(cnij-cnn);
                                f1 = find(cnn>vecC); f2 = find(cnn<=vecC); %if empty means that cnn is not in vecC
                                if(isempty(f2))
                                    cnn_sup = cnn; cnn_inf = vecC(f1(end)); %the two values surrounding cjj %see notes from 27/04, 28/04
                                elseif(isempty(f1))
                                    cnn_sup = vecC(f2(1)); cnn_inf = cnn;
                                else
                                    cnn_sup = vecC(f2(1)); cnn_inf = vecC(f1(end));
                                end
                                step = cnn_sup-cnn_inf;

                                if(diff<step && abs(diff-step)>1e-10)%if diff=0 then cnij=cnn (not in the loop)
                                    % the difference between cnij and cnn is not big enough to tell if one is bigger than the other,
                                    % here we refine vecC
                                    vecC = linspace(cnn_inf,cnn_sup,10);
                                    %disp(['diff<step:', num2str(abs(diff-step)),' i=',num2str(i),', j=',num2str(j)]);
                                    [tab,tabco,tabcn] = findRhohat_ij(vecC);
                                    cs = findCs_ij(vecC,tab);
                                    continue;
                                end
                                if((diff>=step || abs(diff-step)<=1e-10) && ~isempty(cnij))
                                    cnijnotfound=0;
                                    if n==1
                                        c1123 = cnij;
                                        %CEL1(i,j) = c112;
                                    elseif n==2
                                        c2123 = cnij;
                                        %CEL2(i,j) = c212;
                                    elseif n==3
                                        c3123 = cnij;
                                    end
                                end
                            end %end while(cnijnotfound)
                        end
                        %CELJ(i,j) = c112*(cs.order(1)==1) + c212*(cs.order(1)==2); %store costs at which the first disease is eliminated (not the exact value but sufficiently to compare with cjj)
                        %CELI(i,j) = c112*(cs.order(2)==1) + c212*(cs.order(2)==2);
                        elimCostsCVT(iterCurrent,:) = [gamma1,gamma2,gamma3,R10(i),R20(j),R30(k),c1123,c2123,c3123];
                        elimCostsTVT(iterCurrent,:) = [gamma1,gamma2,gamma3,R10(i),R20(j),R30(k),c11,c22,c33];
                        timesTable(iterCurrent) = toc; %store time spent on the current iteration
                    end%end R30
                end %end R2(0)
            end %end R1(0)
        end%end gamma3
    end %end gamma2
end %end gamma1

%%
elimCostsCVT = elimCostsCVT(1:iterCurrent,:);
elimCostsTVT = elimCostsTVT(1:iterCurrent,:);
timesTable = timesTable(1:iterCurrent,:);
histogram(timesTable)

%adding the symmetrical part due to SISxSIS
if nSIS==2 && nSIR==1
    symPart = elimCostsCVT(:,1)~=elimCostsCVT(:,2);
    elimCostsCVT_total = [elimCostsCVT; elimCostsCVT(symPart,[2,1,3,5,4,6,8,7,9])]; %inverting gamma1-gamma2, R1-R2, c1123-c2123
    elimCostsTVT_total = [elimCostsTVT; elimCostsTVT(symPart,[2,1,3,5,4,6,8,7,9])];
end

%filename='C:\Users\Moi\Documents\IPLESP\These\Graphes\strategies_comparison\SIS2SIR\TVT.txt'; %PC
filename1 = '/Users/amandine/Desktop/These/Graphes/strategies_comparison/SIS2SIR/CVT.txt';%mac
filename2 = '/Users/amandine/Desktop/These/Graphes/strategies_comparison/SIS2SIR/TVT.txt';%mac
%writematrix(elimCostsCVT_total, filename1)
%writematrix(elimCostsTVT_total, filename2)

%%
%filename1 = '/Users/amandine/Desktop/These/Graphes/strategies_comparison/SIS2SIR/CVT.txt';%mac
%filename2 = '/Users/amandine/Desktop/These/Graphes/strategies_comparison/SIS2SIR/TVT.txt';%mac

filename1 = 'C:\Users\Moi\Documents\IPLESP\These\Graphes\strategies_comparison\SIS2SIR\CVT.txt'; %PC
filename2 = 'C:\Users\Moi\Documents\IPLESP\These\Graphes\strategies_comparison\SIS2SIR\TVT.txt'; %PC

elimCostsCVT_total = readmatrix(filename1);
elimCostsTVT_total = readmatrix(filename2);

%CVT vs TVT
CVTvsTVT1 = [elimCostsCVT_total(:,[1 2 3 4 5 6]),elimCostsCVT_total(:,7)>elimCostsTVT_total(:,7)];
CVTvsTVT2 = [elimCostsCVT_total(:,[1 2 3 4 5 6]),elimCostsCVT_total(:,8)>elimCostsTVT_total(:,8)];
CVTvsTVT3 = [elimCostsCVT_total(:,[1 2 3 4 5 6]),elimCostsCVT_total(:,9)>elimCostsTVT_total(:,9)];

% aggregate all costs by (R1(0),R2(0),R3(0))
[unR1R2R3, ~, indexUnique] = unique(elimCostsCVT_total(:,[4,5,6]) , 'rows'); %selection by R1,R2,R3
firstColumn  = accumarray(indexUnique, CVTvsTVT1(:,7) , size(unR1R2R3(:,1)) , @(x) sum(x));
secondColumn = accumarray(indexUnique, CVTvsTVT2(:,7) , size(unR1R2R3(:,1)) , @(x) sum(x));
thirdColumn  = accumarray(indexUnique, CVTvsTVT3(:,7) , size(unR1R2R3(:,1)) , @(x) sum(x));
nColumn1     = accumarray(indexUnique, CVTvsTVT1(:,7) , size(unR1R2R3(:,1)) , @(x) length(x));
nColumn2     = accumarray(indexUnique, CVTvsTVT2(:,7) , size(unR1R2R3(:,1)) , @(x) length(x));
nColumn3     = accumarray(indexUnique, CVTvsTVT3(:,7) , size(unR1R2R3(:,1)) , @(x) length(x));

outputArray  = cat(2 , unR1R2R3, firstColumn,secondColumn,thirdColumn,nColumn1,nColumn2,nColumn3);

% colArray = outputArray(:,4)./outputArray(:,7)*[0,1,0]
% plot3(outputArray(:,1),outputArray(:,2),outputArray(:,3), 'o', 'Color',colArray')

%disease 1
figure(1)
scatter3(outputArray(:,1),outputArray(:,2),outputArray(:,3), 200, outputArray(:,4)./outputArray(:,7),'filled')
xlabel('R_1(0)')
ylabel('R_2(0)')
zlabel('R_3(0)')
cb = colorbar();
title({['color: proportion of cases where the cost to eliminate disease 1'],...
    ['is lower with CVT than with TVT'],...
    ['SIS(1)xSIS(2)xSIR(3) :  \gamma_1(0),\gamma_2(0),\gamma_3(0)=1/D such that D=1/12:1/12:6/12']})
%disease 2
figure(2)
scatter3(outputArray(:,1),outputArray(:,2),outputArray(:,3), 200, outputArray(:,5)./outputArray(:,8),'filled')
xlabel('R_1(0)')
ylabel('R_2(0)')
zlabel('R_3(0)')
cb = colorbar();
title({['color: proportion of cases where the cost to eliminate disease 2'],...
    ['is lower with CVT than with TVT'],...
    ['SIS(1)xSIS(2)xSIR(3) :  \gamma_1(0),\gamma_2(0),\gamma_3(0)=1/D such that D=1/12:1/12:6/12']})

%disease 3
figure(3)
scatter3(outputArray(:,1),outputArray(:,2),outputArray(:,3), 200, outputArray(:,6)./outputArray(:,9),'filled')
xlabel('R_1(0)')
ylabel('R_2(0)')
zlabel('R_3(0)')
cb = colorbar();
title({['color: proportion of cases where the cost to eliminate disease 3'],...
    ['is lower with CVT than with TVT'],...
    ['SIS(1)xSIS(2)xSIR(3) :  \gamma_1(0),\gamma_2(0),\gamma_3(0)=1/D such that D=1/12:1/12:6/12']})
