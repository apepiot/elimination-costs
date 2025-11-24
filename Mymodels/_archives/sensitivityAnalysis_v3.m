%november 2023
clear all; close all;

newSet   = 0;
oldSetNo = 8;
vecRho=0.:0.05:0.3;
%-------------------%

%Evaluating the impact of the different voluntary testing rates on the
%prevalences of each infection with the four-disease model (version v5)

resultsPath = 'C:/Users/Moi/Documents/IPLESP/These/Codes/Analyse/modele_v5';

if newSet
    %lastParameterSet number
    fileID = fopen([resultsPath,'/parameters/lastParameterSet.txt'],'r');
    lastParameterSetNo = fscanf(fileID,'%f');
    fclose(fileID);

    %Adding a new parameter set
    newParameterSetNo = lastParameterSetNo+1;
    b=100;

    [paramTab,mu,~] = sampleParameters_v3(true,true,true,true,b);
    paramTab{3}.eta = 4; paramTab{3}.zeta=randPERT(46,60,71,1)/100; paramTab{3}.p = 0.;
else
    %open the old set
    Recuppath = [resultsPath,'/parameters/set_',num2str(oldSetNo)];
    paramTab{1} = table2struct(readtable([Recuppath,'_Ct.txt']));
    paramTab{2} = table2struct(readtable([Recuppath,'_Ng.txt']));
    paramTab{3} = table2struct(readtable([Recuppath,'_HIV.txt']));
    paramTab{4} = table2struct(readtable([Recuppath,'_syphilis.txt']));
    mu=1/35; b=100;
    %newParameterSetNo = 
end


%% Evaluating the impact of rho_i on the prevalences 
% Initialize the vector of rho's
createParamRho; %all the rhos at 0
paramRho.eta_c_prep=4;
paramRho.eta_g_prep=4;
paramRho.eta_s_prep=1;
paramRho.eta_h_prep=4;

paramTab{4}.rhob = 0; %paramRho.rho_s=paramTab{4}.rhob;
paramTab{1}.rhob = 0.0; %paramRho.rho_c=paramTab{1}.rhob;
paramTab{2}.rhob = 0; %paramRho.rho_g=paramTab{2}.rhob;
paramTab{3}.rhob = 0; %paramRho.rho_h=paramTab{3}.rhob;

% loop of rho_i
vecRho_i = vecRho;
%vecRho_i=0.08
timesiter = [];
P = zeros(length(vecRho_i),4,4);
j=0; %number of the subplot
for k=[3,4,1,2]
    i=1; j=j+1;
    for rho=vecRho_i
        disp(['rho=',num2str(rho)])
        paramTab{k}.rhob = rho;
        paramRho.rho_s = paramTab{4}.rhob; %paramTab{4}.rhob=paramRho.rho_s;
        paramRho.rho_c = paramTab{1}.rhob; %paramTab{1}.rhob=paramRho.rho_c;
        paramRho.rho_g = paramTab{2}.rhob; %paramTab{2}.rhob=paramRho.rho_g;
        paramRho.rho_h = paramTab{3}.rhob; %paramTab{3}.rhob=paramRho.rho_h;
    
        tic
        %paramRho.rho = rhoh; paramTab{3}.rhob=rhoh;
        [ES,~,~,nbIter,~,changeSolver] = P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,1,...
            'knitroampl',15,{paramTab{k}.disease});
        [HIV,syphilis,Ct,Ng,tot,N] = assigningPrevalence(ES);
        P(i,1,k) = sum(ES(HIV.asymptomaticHIV))/N;
        P(i,2,k) = sum(ES(syphilis.asymptomatic))/N;
        P(i,3,k) = sum(ES(Ct.asymptomatic))/N;
        P(i,4,k) = sum(ES(Ng.asymptomatic))/N;
        
        [ES_sictp,ES_seiiis,ES_ct,ES_ng] = singleDisES(paramTab,mu,b);
        sum(ES_sictp([2:3,5:6]))/N
        sum(ES_seiiis(2:5))/N
        sum(ES_ct([2,3]))/N
        sum(ES_ng([2,3]))/N
        
        timesiter(i,1,k) = toc;
        timesiter(i,2,k) = changeSolver;
        i=i+1;
    end
    figure(1)
    subplot(2,2,j)
    plot(vecRho_i,P(:,1,k),'r','LineWidth',2);
    hold on
    plot(vecRho_i,P(:,2,k),'y','LineWidth',2);
    plot(vecRho_i,P(:,3,k),'b','LineWidth',2);
    plot(vecRho_i,P(:,4,k),'g','LineWidth',2);
    legend('$\Pi_{h}$','$\Pi_{s}$','$\Pi_{c}$','$\Pi_{g}$','Interpreter','latex')
    xlabel(['$\rho_{',paramTab{k}.mini_d,'}$'],'Interpreter','latex')
    paramTab{k}.rhob = 0;
end

%% Evaluating the impact of rho_ij on the prevalences (+add rho_i)
paramRho.rho_s = paramTab{4}.rhob; %paramTab{4}.rhob=paramRho.rho_s;
paramRho.rho_c = paramTab{1}.rhob; %paramTab{1}.rhob=paramRho.rho_c;
paramRho.rho_g = paramTab{2}.rhob; %paramTab{2}.rhob=paramRho.rho_g;
paramRho.rho_h = paramTab{3}.rhob; %paramTab{3}.rhob=paramRho.rho_h;
vecRho_i = vecRho;
infectionsNo = [3,4,1,2];
for i=1:4
    k=infectionsNo(i);
    
    %P_j(rho_k)
    figure(2)
    subplot(4,4,4*(i-1)+i)
    plot(vecRho_i,P(:,1,k),'r','LineWidth',2);
    hold on
    plot(vecRho_i,P(:,2,k),'y','LineWidth',2);
    plot(vecRho_i,P(:,3,k),'b','LineWidth',2);
    plot(vecRho_i,P(:,4,k),'g','LineWidth',2);
    %legend('$\Pi_{h}$','$\Pi_{s}$','$\Pi_{c}$','$\Pi_{g}$','Interpreter','latex')
    xlabel(['$\rho_{',paramTab{k}.mini_d,'}$'],'Interpreter','latex')
    
    %P_k(rho_ij)
    for j=(i+1):4
        ell=infectionsNo(j);
        Pij=[];m=0;
        
        for rho=vecRho
            m=m+1;
            %str2sym(['paramRho.rho_',paramTab{ell}.mini_d,paramTab{k}.mini_d]) = rho;
            paramRho.(['rho_',paramTab{k}.mini_d,paramTab{ell}.mini_d]) = rho;
            [ES,~,~,~,~,changeSolver] = P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,1,...
                'knitroampl',25);
            [HIV,syphilis,Ct,Ng,tot,N] = assigningPrevalence(ES);
            Pij(m,1) = sum(ES(HIV.asymptomaticHIV))/N;
            Pij(m,2) = sum(ES(syphilis.asymptomatic))/N;
            Pij(m,3) = sum(ES(Ct.asymptomatic))/N;
            Pij(m,4) = sum(ES(Ng.asymptomatic))/N;
        end
        paramRho.(['rho_',paramTab{k}.mini_d,paramTab{ell}.mini_d]) = 0;
        
        figure(2)
        subplot(4,4,4*(i-1)+j)
        plot(vecRho,Pij(:,1),'r','LineWidth',2);
        hold on
        plot(vecRho,Pij(:,2),'y','LineWidth',2);
        plot(vecRho,Pij(:,3),'b','LineWidth',2);
        plot(vecRho,Pij(:,4),'g','LineWidth',2);
        xlabel(['$\rho_{',paramTab{k}.mini_d,paramTab{ell}.mini_d,'}$'],'Interpreter','latex')
    end
end


%% Evaluating the impact of rho_ijk
kits3 = {'hsc','hsg','hcg','scg'};
i=0;
for k = kits3
    i=i+1;
    P=[];
    m=0;
    for rho=vecRho
        m=m+1;
        paramRho.(['rho_',k{:}]) = rho;
        [ES,~,~,~,~,changeSolver] = P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,1,...
            'knitroampl',25);
        [HIV,syphilis,Ct,Ng,tot,N] = assigningPrevalence(ES);
        P(m,1) = sum(ES(HIV.asymptomaticHIV))/N;
        P(m,2) = sum(ES(syphilis.asymptomatic))/N;
        P(m,3) = sum(ES(Ct.asymptomatic))/N;
        P(m,4) = sum(ES(Ng.asymptomatic))/N;
    end
    paramRho.(['rho_',k{:}]) = 0;

    figure(3)
    subplot(2,2,i)
    plot(vecRho,P(:,1),'r','LineWidth',2);
    hold on
    plot(vecRho,P(:,2),'y','LineWidth',2);
    plot(vecRho,P(:,3),'b','LineWidth',2);
    plot(vecRho,P(:,4),'g','LineWidth',2);
    xlabel(['$\rho_{',k{:},'}$'],'Interpreter','latex')
end


%% Evaluating the impact of rho_ijkl on the prevalences
P=[];
vecRho_ijkl = vecRho;
m=0;
for rho=vecRho_ijkl
    m=m+1;
    paramRho.rho_hscg = rho;
    [ES,~,~,~,~,changeSolver] = P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,1,...
        'knitroampl',25);
    [HIV,syphilis,Ct,Ng,tot,N] = assigningPrevalence(ES);
    P(m,1) = sum(ES(HIV.asymptomaticHIV))/N;
    P(m,2) = sum(ES(syphilis.asymptomatic))/N;
    P(m,3) = sum(ES(Ct.asymptomatic))/N;
    P(m,4) = sum(ES(Ng.asymptomatic))/N;
end
paramRho.rho_hscg = 0;

figure(4)
plot(vecRho_ijkl,P(:,1),'r','LineWidth',2);
hold on
plot(vecRho_ijkl,P(:,2),'y','LineWidth',2);
plot(vecRho_ijkl,P(:,3),'b','LineWidth',2);
plot(vecRho_ijkl,P(:,4),'g','LineWidth',2);
xlabel(['$\rho_{hscg}$'],'Interpreter','latex')

oldP = paramTab{3}.p;
vecP = 0:0.1:1;
m=0;
P=[];
for p=vecP
    m=m+1;
    paramTab{3}.p = p;
    [ES,~,~,~,~,changeSolver] = P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,1,...
        'knitroampl',25);
    [HIV,syphilis,Ct,Ng,tot,N] = assigningPrevalence(ES);
    P(m,1) = sum(ES(HIV.asymptomaticHIV))/N;
    P(m,2) = sum(ES(syphilis.asymptomatic))/N;
    P(m,3) = sum(ES(Ct.asymptomatic))/N;
    P(m,4) = sum(ES(Ng.asymptomatic))/N;
end
paramTab{3}.p = oldP;

figure(5)
plot(vecP,P(:,1),'r','LineWidth',2);
hold on
plot(vecP,P(:,2),'y','LineWidth',2);
plot(vecP,P(:,3),'b','LineWidth',2);
plot(vecP,P(:,4),'g','LineWidth',2);
xlabel(['$\rho_{hscg}$'],'Interpreter','latex')

%%

if newSet
    figure(1)
    saveas(gcf,[resultsPath,'/figures/set_',num2str(newParameterSetNo),'_P(rho_i).png'])
    figure(2)
    saveas(gcf,[resultsPath,'/figures/set_',num2str(newParameterSetNo),'_P(rho_ij).png'])
    figure(3)
    saveas(gcf,[resultsPath,'/figures/set_',num2str(newParameterSetNo),'_P(rho_ijk).png'])
    figure(4)
    saveas(gcf,[resultsPath,'/figures/set_',num2str(newParameterSetNo),'_P(rho_ijkl).png'])
    figure(5)
    saveas(gcf,[resultsPath,'/figures/set_',num2str(newParameterSetNo),'_P(p).png'])

    writetable(struct2table(paramTab{1}), [resultsPath,'/parameters/set_',num2str(newParameterSetNo),'_Ct.txt'])
    writetable(struct2table(paramTab{2}), [resultsPath,'/parameters/set_',num2str(newParameterSetNo),'_Ng.txt'])
    writetable(struct2table(paramTab{3}), [resultsPath,'/parameters/set_',num2str(newParameterSetNo),'_HIV.txt'])
    writetable(struct2table(paramTab{4}), [resultsPath,'/parameters/set_',num2str(newParameterSetNo),'_syphilis.txt'])
    writetable(struct2table(paramRho), [resultsPath,'/parameters/set_',num2str(newParameterSetNo),'_rhos.txt'])

    %Updating the number of the last parameter set
    writematrix(newParameterSetNo, [resultsPath,'/parameters/lastParameterSet.txt'])
end