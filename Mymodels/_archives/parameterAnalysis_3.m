
%% 6. Impact of P_prev_base on P_und_prep_base
close all
nbSim =1500;
mu=1/30.6; b=100;
vecvalPHIV = [10/100,15/100,20/100,25/100,30/100,35/100];
BinWidth=0.001;
fig = figure(7)
for valPHIV=vecvalPHIV
    Pund_base=[];
    for i = 1:nbSim
        PHIV_prevagay_VT = valPHIV;%randPERT(12.0,14.3,16.9,1)/100;
        R_sict_VT      = 1/(1-PHIV_prevagay_VT);
        sigma_sict     = 52/randPERT(6.7,8.2,9.8,1);
        rhobh          = 1/randPERT(1.6,1.9,2.2,1);
        theta_sict0     = 1/randPERT(4,4.4,10,1);
        theta_sict      = theta_sict0+rhobh;
        gamma_sict0     = 0;
        gamma_sict      = rhobh+gamma_sict0;
        ratioBeta       = randPERT(8.4,9.1,9.6,1);
        betaC_sict      = R_sict_VT*(sigma_sict+gamma_sict+mu)*(theta_sict+mu)/(ratioBeta*(theta_sict+mu)+sigma_sict); %check
        betaI_sict = ratioBeta*betaC_sict; %check
        [R0_sict,~,alpha_sict] = Rp_SICR_v4(betaI_sict,betaC_sict,theta_sict0,sigma_sict,gamma_sict0,mu,b,0); %check
        %P0_sict = 1-1/R0_sict;
        %undiagnosed prevalence :
        Pund_base(i) = mu*(theta_sict+sigma_sict+mu)/(betaI_sict*(theta_sict+mu)+betaC_sict*sigma_sict)*(R_sict_VT-1);
        %Pund_0    = mu*(theta_sict0+sigma_sict+mu)/(betaI_sict*(theta_sict0+mu)+betaC_sict*sigma_sict)*(R0_sict-1);
    end
    histogram(Pund_base,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
    hold on
end
xlim([0,0.03])
legend({'$\Pi(\rho^0_h)=10\%$','$\Pi(\rho^0_h)=15\%$','$\Pi(\rho^0_h)=20\%$',...
    '$\Pi(\rho^0_h)=25\%$','$\Pi(\rho^0_h)=30\%$','$\Pi(\rho^0_h)=35\%$' },'Interpreter','latex')
xlabel('$\Pi_{h,p=0}^{und}(\rho^0_h)$','Interpreter','latex')
saveas(fig,[backupPath,'hist_Ph_und_base_acc_to_Ph_base.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_Ph_und_base_acc_to_Ph_base.pdf'])


%% 7. Impact of P_prev_base on P_und_prep_base
close all
nbSim =1500;
mu=1/30.6; b=100;
vecRhobh = 1./[1,2,3,4,5];
BinWidth=0.001;
fig = figure(8)
for valRhobh=vecRhobh
    Pund_base=[];
    for i = 1:nbSim
        PHIV_prevagay_VT = randPERT(12.0,14.3,16.9,1)/100;
        R_sict_VT      = 1/(1-PHIV_prevagay_VT);
        sigma_sict     = 52/randPERT(6.7,8.2,9.8,1);
        rhobh          = valRhobh;%1/randPERT(1.6,1.9,2.2,1);
        theta_sict0     = 1/randPERT(4,4.4,10,1);
        theta_sict      = theta_sict0+rhobh;
        gamma_sict0     = 0;
        gamma_sict      = rhobh+gamma_sict0;
        ratioBeta       = randPERT(8.4,9.1,9.6,1);
        betaC_sict      = R_sict_VT*(sigma_sict+gamma_sict+mu)*(theta_sict+mu)/(ratioBeta*(theta_sict+mu)+sigma_sict); %check
        betaI_sict = ratioBeta*betaC_sict; %check
        [R0_sict,~,alpha_sict] = Rp_SICR_v4(betaI_sict,betaC_sict,theta_sict0,sigma_sict,gamma_sict0,mu,b,0); %check
        %P0_sict = 1-1/R0_sict;
        %undiagnosed prevalence :
        Pund_base(i) = mu*(theta_sict+sigma_sict+mu)/(betaI_sict*(theta_sict+mu)+betaC_sict*sigma_sict)*(R_sict_VT-1);
        %Pund_0    = mu*(theta_sict0+sigma_sict+mu)/(betaI_sict*(theta_sict0+mu)+betaC_sict*sigma_sict)*(R0_sict-1);
    end
    histogram(Pund_base,'EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
    hold on
end
xlim([0,0.03])
legend({'$\rho^0_h=1$','$\rho^0_h=1/2$','$\rho^0_h=1/3$',...
    '$\rho^0_h=1/4$','$\rho^0_h=1/5$' },'Interpreter','latex')
xlabel('$\Pi_{h,p=0}^{und}(\rho^0_h)$','Interpreter','latex')
saveas(fig,[backupPath,'hist_Ph_und_base_acc_to_rho_base.png'])
set(fig,'PaperOrientation','landscape');
saveas(fig,[backupPath,'hist_Ph_und_base_acc_to_rho_base.pdf'])


%%prevalence sans les traités ? utiliser la fonction scitp qui rend les
%%prevalences, en mettant p a 0
[Rp,Lambdap,alpha,Ptot,Pun] = Rp_SICTP(betaI_sict,betaC_sict,theta_sict0,gamma_sict0,sigma_sict,0,0,0,mu,b,rhobh)
