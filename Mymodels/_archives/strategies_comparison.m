%% Strategy comparison
%% code for 2 diseases
clear all;
elimCosts =[]; %prop = []; 
s1=1;s2=1;
mu = 1/35;
nSIS = 1;nSIR=1; N=nSIS+nSIR;

len = 15; %sampling R(0)'s
R10 = linspace(1.1,5,len);
R20 = linspace(1.1,5,len);
%gamma2=365/21; gamma1=12/1; 
d1 = 1/12:1/12:6/12;
for gamma1 = 1./d1
    if(nSIS==2)
        d2 = 1/12:1/12:1/gamma1; %the problem is symmetrical
    else
        d2 = d1;
    end
    for gamma2 = 1./d2
        disp(['gamma1=',num2str(gamma1),' gamma2=',num2str(gamma2)])
   
        [R1,R2] = meshgrid(R10,R20);
        BETA1 = R1.*(gamma1+mu);
        BETA2 = R2.*(gamma2+mu);

        ALPHA1 = BETA1./s1.*(1-1./R1);
        ALPHA2 = BETA2./s2.*(1-1./R2);
        C11 = 1./R1-1; 
        C22 = (1./R2-1)*(nSIS==2)+(-mu./BETA2.*(R2-1))*(nSIR==1 & nSIS==1);
        
        CELI=zeros(len); CELJ=zeros(len); %store cs/the second/first disease is eliminated
        CEL1=zeros(len);CEL2=zeros(len); %store cs/ disease 1 (resp. 2) is eliminated
        for i=1:len
            %disp(['i=',num2str(i)])
            for j=1:len
                c11 = C11(i,j); c01 = -c11;
                c22 = C22(i,j); c02 = -c22*(nSIS==2) + mu/BETA2(i,j)*(1-1/R2(i,j))*(nSIR==1 & nSIS==1);
                %j : first disease eliminated
                
                mincnn = min(c11,c22);maxc0n = max(c01,c02);
                maxcnn = max(c11,c22);

                %interval of c
                vecC = linspace(1.5*mincnn,maxc0n,50); 
                diff = 0;
                step = vecC(2)-vecC(1);
                vecAlpha = [ALPHA1(i,j),ALPHA2(i,j)]; %alphaj = min(vecAlpha);
                %cjj = c11*(ALPHA1(i,j)<ALPHA2(i,j)) + c22*(ALPHA1(i,j)>=ALPHA2(i,j));

                % we look for the cost of each disease elimination
                paramSIS = [BETA1(i,j),gamma1,1;BETA2(i,j),gamma2,s2];
                paramSIR = [BETA2(i,j),gamma2,s2]*(nSIS==2 & nSIR==1);
                findRhohat_ij = @(c) findRhohat(nSIS,nSIR,0,paramSIS(1:nSIS,:),paramSIR(1:nSIR,:),[],mu,5,c); 
                findCs_ij = @(c,tabRho) findThresholds(nSIS,nSIR,0,tabRho,vecAlpha,[BETA1(i,j),BETA2(i,j)],[gamma1,gamma2],[s1,s2],5,mu,c);
                [tab,tabco,tabcn] = findRhohat_ij(vecC);        
                cs_init = findCs_ij(vecC,tab);
                tab_init = tab;
                vecC_init = vecC;
                elimOrder = cs_init.order;

                for n = 1:N %n = disease number
                    %disp(['nDis:',num2str(n)])
                    cnijnotfound=1; %cnij : cost at which CVT eliminates n (n in {i,j})
                    %below : to determine either cnn<cnij or cnn>cnij,
                    %cnij needs to be sufficiently accurate
                    cnn = c11*(n==1) + c22*(n==2);
                    cs = cs_init; %thresholds found during the first run of findThresholds (reinitializing for each disease to avoid crushing of the following cnij)
                    tab = tab_init;
                    alphan = vecAlpha(n);

                    while (cnijnotfound)
                        %disp('cnijnotfound')
                        if n==elimOrder(1) %if n is the first disease eliminated
                            cnij = cs.cs1dis;
                        elseif n==elimOrder(2) %if n is the second eliminated
                            cnij = cs.cs2dis;
                        end

                        %cnij = (cs.cs1dis)*(cs.order(1)==n) + (cs.cs2dis)*(cs.order(2)==n); %sum(A) returns 0 if A is empty
                        if(isempty(cnij)) %vecC does not include cnij
                            %disp(['cnij empty, i=',num2str(i),', j=',num2str(j)]);
                            if(tab.rhohat(1)<alphan) %vecC did not go far enough on the left side
                                vecC = linspace(vecC(1)-10*step, vecC(1),10);
                                [tab,tabco,tabcn] = findRhohat_ij(vecC);        
                                cs = findCs_ij(vecC,tab);
                                %disp('if tab.rhohat(1)<alphan');
                                continue; %next iteration of the while loop
                            elseif (tab.rhohat(end)>=alphan)
                                %disp('if tab.rhohat(end)>alphan')
                                vecC = linspace(vecC(end),vecC(end)+10*step,10);
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
                                c112 = cnij;
                                %CEL1(i,j) = c112;
                            elseif n==2
                                c212 =cnij;
                                %CEL2(i,j) = c212;
                            end
                        end
                    end %end while(cnijnotfound)
                end
                %CELJ(i,j) = c112*(cs.order(1)==1) + c212*(cs.order(1)==2); %store costs at which the first disease is eliminated (not the exact value but sufficiently to compare with cjj)          
                %CELI(i,j) = c112*(cs.order(2)==1) + c212*(cs.order(2)==2);
                elimCosts = [elimCosts;gamma1,gamma2,R1(i,j),R2(i,j),c112,c212,c11,c22];
            end %end R2(0)
        end %end R1(0)
        
        %matrices that stores the relative difference of gain obtained with CVT
        
%          dis2drives1 = (ALPHA1<=ALPHA2).*(CELJ-C11)./abs(C11); %
%          dis1drives2 = (ALPHA1>=ALPHA2).*(CELJ-C22)./abs(C22);
%          tot = dis2drives1+dis1drives2; minv = min(min(tot));maxv=max(max(tot));
%          dis2drives1((ALPHA1>ALPHA2) & dis2drives1==0)=NaN;
%          dis1drives2((ALPHA1<ALPHA2) & dis1drives2==0)=NaN;
        
%         sizeM = size(tot,1)*size(tot,2);
%         prop = [prop; gamma1,gamma2,sum(sum(dis2drives1>0,'omitnan')),...
%             sum(sum(dis2drives1<=0,'omitnan')),...
%             sum(sum(dis1drives2>0,'omitnan')),...
%             sum(sum(dis1drives2<=0,'omitnan'))];

        %type = 'SIS2';
        %saveFiles = false;
        %plot_3Dgraph; %one disease drives the other (one graph for each couple (gamma1,gamma2)  
    end %end gamma2
end %end gamma1

%% adding the other "side" of gamma's triangle %SISxSIS
% such that d2>1/gamma1; %see notes from 28/04
% gamma1<-gamma2
% gamma2<-gamma1
% dis1drives2<-dis2drives1
% dis2drivesdis1
% attention: do not include cases where gamma1=gamma2

%prop2 = prop(prop(:,1)~=prop(:,2),[2,1,5,6,3,4]);
%proptot = [prop;prop2];
% filename='C:\Users\Moi\Documents\IPLESP\These\Graphes\strategies_comparison\SIS2\SIS2.txt'; %PC
% writematrix(proptot, filename)

if nSIS==2 && nSIR==0
    elimCosts_symPart = elimCosts(elimCosts(:,1)~=elimCosts(:,2),:);
    elimCosts_total = [elimCosts; elimCosts_symPart(:,[2,1,4,3,6,5,8,7])];
end

%filename='C:\Users\Moi\Documents\IPLESP\These\Graphes\strategies_comparison\SISSIR\SISSIR_comp_cnij_cnn.txt'; %PC
%writematrix(elimCosts_total, filename)

%%
%filename = '/Users/amandine/Desktop/These/Graphes/strategies_comparison/SISSIR/SISSIR_comp_cnij_cnn.txt';%mac
filename = 'C:\Users\Moi\Documents\IPLESP\These\Graphes\strategies_comparison\SIS2\SIS2_comp_cnij_cnn.txt'; %PC
elimCosts_total = readmatrix(filename);

% c112-c11... >0 : 1, <0 : -1
% aggregate by R1, R2, gamma1, gamma2

elimCosts_total(:,9)  = elimCosts_total(:,5)>elimCosts_total(:,7); %1 : CVT better for disease 1 than TVT
elimCosts_total(:,10) = elimCosts_total(:,6)>elimCosts_total(:,8); %1 : CVT better for disease 2 than TVT

% aggregate all costs by (R1(0),R2(0))
[unR1R2, ~, indexUnique] = unique(elimCosts_total(:,[3,4]) , 'rows'); %selection by R1,R2
firstColumn  = accumarray(indexUnique, elimCosts_total(:,9) , size(unR1R2(:,1)) , @(x) sum(x));
secondColumn = accumarray(indexUnique, elimCosts_total(:,10) , size(unR1R2(:,1)) , @(x) sum(x));
%sumColumn    = accumarray(indexUnique, elimCosts_total(:,3) , size(unR1R2(:,1)) , @(x) sum(x))
nColumn1    = accumarray(indexUnique, elimCosts_total(:,9) , size(unR1R2(:,1)) , @(x) length(x));
nColumn2    = accumarray(indexUnique, elimCosts_total(:,10) , size(unR1R2(:,1)) , @(x) length(x));
outputArray  = cat(2 , unR1R2, firstColumn,secondColumn,nColumn1,nColumn2)

%x = outputArray(:,1);
%y = outputArray(:,2);
%s = scatter(x,y,'filled');
%distfromzero = outputArray(:,3);
%s.AlphaData = distfromzero;
%s.MarkerFaceAlpha = 'flat';

 len =10;
R10 = linspace(1.1,5,len);
R20 = linspace(1.1,5,len);
[R1,R2] = meshgrid(R10,R20);


figure(1)
Z1 = reshape(outputArray(:,3)./nColumn1,length(unique(outputArray(:,1))),length(unique(outputArray(:,2)))); %bycolumn (R1 first)
s1 = surf(R1,R2,Z1,'FaceColor','interp');
s1.EdgeColor = 'none';
xlabel('R_1(0)')
ylabel('R_2(0)')
view(2)
title({['color: proportion of cases where the cost to eliminate disease 1'],...
    ['is lower with CVT than with TVT'],...
    ['SIS(1)xSIS(2) : \gamma_1(0),\gamma_2(0)=1/D such that D=1/12:1/12:6/12']})
colorbar;

figure(2)
Z2 = reshape(outputArray(:,4)./nColumn2,length(unique(outputArray(:,1))),length(unique(outputArray(:,2)))); %bycolumn
s2 = surf(R1,R2,Z2,'FaceColor','interp');
s2.EdgeColor = 'none';
xlabel('R_1(0)')
ylabel('R_2(0)')
view(2)
title({['color: proportion of cases where the cost to eliminate disease 2'],...
    ['is lower with CVT than with TVT'],...
    ['SIS(1)xSIS(2) :  \gamma_1(0),\gamma_2(0)=1/D such that D=1/12:1/12:6/12']})
colorbar
