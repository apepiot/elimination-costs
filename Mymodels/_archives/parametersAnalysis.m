close all; clear all;
b='NaN'; N=10000;
alphaCt = zeros(1,N);alphaNg=alphaCt;alphaHIV=alphaCt;alphaS=alphaCt;

for k=1:N
    [paramTab,mu,vecAlphas] = sampleParameters_v2(true,true,true,true,b);
    alphaCt(k) = vecAlphas(1);
    alphaNg(k) = vecAlphas(2);
    alphaHIV(k) = vecAlphas(3);
    alphaS(k) = vecAlphas(4);
end
%%
close all
alphaTot(1,:) = alphaCt; alphaTot(2,:)=alphaNg;alphaTot(3,:)=alphaHIV;alphaTot(4,:)=alphaS;

fig = figure(1)
BinWidth=0.01%(max(alphaTot(:))-min(alphaTot(:)))/60;
histogram(alphaHIV,'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV

hold on

%BinWidth=min([10,(max(alphaS)-min(alphaS))./10]);
histogram(alphaS,'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis

hold on

%BinWidth=min((max(alphaCt)-min(alphaCt))./10);
histogram(alphaCt,'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct

hold on

%BinWidth=min(max(alphaNg)-min(alphaNg))./10);
histogram(alphaNg,'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng

%set(gca,'ytick',[])
%set(gca,'fontsize',10)
xlim([-0.01, 1])
legend('$\rho_{h}^\prime$',...
    '$\rho_{s}^\prime$',...
    '$\rho_{c}^\prime$',...
    '$\rho_{g}^\prime$',...
    'Interpreter','latex','FontSize',12,'Box','off')
%%
path = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\graphes\';
saveas(fig,[path,'hist_rhos.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[path,'hist_rhos.pdf'])
%%

% For LaTex / table with distributions of rho
pathW = 'C:\Users\Moi\Documents\IPLESP\These\Rapports\Archives\tables\';
fileID = fopen([pathW,'rhos.txt'],'w');

%Initialization
begin = {'\begin{table}[h]','\centering', ...
    '\begin{tabular}{|l|c|}','\hline'};
for k=1:length(begin)
    fprintf(fileID,'%12s\r\n',begin{k});
end
fprintf(fileID,'%12s\r\n',' & mean [95\% CI] \\ \hline \hline');

fprintf(fileID,'%15s',['$\rho_{h}^\prime$ & ', num2str(round(mean(alphaHIV),3)),...
                                            '~[',num2str(round(prctile(alphaHIV,2.5),3)),',',...
                                            num2str(round(prctile(alphaHIV,97.5),3)),']',...
                                            ' \\ \hline']);
fprintf(fileID,'%15s',['$\rho_{s}^\prime$ & ', num2str(round(mean(alphaS),3)),...
                                            '~[',num2str(round(prctile(alphaS,2.5),3)),',',...
                                            num2str(round(prctile(alphaS,97.5),3)),']',...
                                            ' \\ \hline']);
fprintf(fileID,'%15s',['$\rho_{c}^\prime$ & ', num2str(round(mean(alphaCt),3)),...
                                            '~[',num2str(round(prctile(alphaCt,2.5),3)),',',...
                                            num2str(round(prctile(alphaCt,97.5),3)),']',...
                                            ' \\ \hline']);
fprintf(fileID,'%15s',['$\rho_{g}^\prime$ & ', num2str(round(mean(alphaNg),3)),...
                                            '~[',num2str(round(prctile(alphaNg,2.5),3)),',',...
                                            num2str(round(prctile(alphaNg,97.5),3)),']',...
                                            ' \\ \hline']);


fprintf(fileID,'%12s\r\n','\end{tabular}');
fprintf(fileID,'%12s\r\n','\end{table}');

fclose(fileID);

%%
clear all; close all;
b=2;
N=1000;
R=zeros(N,4);alpha=zeros(N,4);
for i=1:N
    [paramTab,mu,vecAlphas] = sampleParameters_v3(true,true,true,true,b);
    paramCt = paramTab{1};
    [R(i,1),~,alpha(i,1)] = Rp_SEIIS_v4(paramCt.beta,paramCt.nu,paramCt.eps,paramCt.sigma,paramCt.gamma,mu,b,paramCt.rhob);
    paramNg = paramTab{2};
    [R(i,2),~,alpha(i,2)] = Rp_SEIIS_v4(paramNg.beta,paramNg.nu,paramNg.eps,paramNg.sigma,paramNg.gamma,mu,b,paramNg.rhob);
    paramHIV= paramTab{3};
    [R(i,3),~,alpha(i,3)] = Rp_SICR_v4(paramHIV.betaI,paramHIV.betaC,paramHIV.theta,paramHIV.sigma,paramHIV.gamma,mu,b,paramHIV.rhob);
    paramS  = paramTab{4};
    [R(i,4),~,alpha(i,4)] = Rp_SEIIIS_v4(paramS.beta,paramS.sigma,paramS.tau,paramS.nu,paramS.gamma1,paramS.theta,paramS.gamma3,mu,b,paramS.rhob);
end

figure(1)
histogram(R(:,3))
title('R')

figure(2)
histogram(1-1./R(:,3))
title('prevalence of HIV')

%%
clear all; close all;
b = 2;
N = 1000;
vecP = [0,0.25,0.5,0.75];
R = zeros(N,length(vecP)); alpha = zeros(N,length(vecP));
k = 0;
zeta = 0.6; eta = 4;
for p=vecP
    k=k+1;
    for i=1:N
        [paramTab,mu,~] = sampleParameters_v3(false,false,true,false,b);
        paramHIV = paramTab{1};
        [R(i,k),~,alpha(i,k)] = Rp_SICTP(paramHIV.betaI,paramHIV.betaC,...
            paramHIV.theta,paramHIV.sigma,zeta,eta,p,mu,b,paramHIV.rhob);
    end    
end

figure(1)
histogram(R(:,1))
hold on
histogram(R(:,2))
histogram(R(:,3))
histogram(R(:,4))
plot([1 1],[0,N/5],'r-')
legend(['p=',num2str(vecP(1))],...
    ['p=',num2str(vecP(2))],...
    ['p=',num2str(vecP(3))],...
    ['p=',num2str(vecP(4))])
xlabel('$R(\rho_0)$','Interpreter','latex')

figure(2)
histogram(alpha(:,1))
hold on
histogram(alpha(:,2))
histogram(alpha(:,3))
histogram(alpha(:,4))

legend(['p=',num2str(vecP(1))],...
    ['p=',num2str(vecP(2))],...
    ['p=',num2str(vecP(3))],...
    ['p=',num2str(vecP(4))])
xlabel('$\rho^\prime$','Interpreter','latex')

