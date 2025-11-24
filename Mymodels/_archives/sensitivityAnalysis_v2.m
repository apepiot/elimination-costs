%% Sensitivity analysis on the model
% Finding the proportions of the prevalences that go to 0, according to p_h
%tabRecap=[];
%clearvars -except tabRecap; 
clear all; close all;
tabRecap=[];
n=100; b=2;
P_inf_tot = zeros(n,4);
vecpHIV=[0.3,0.4,0.5,0.6,0.7];
k=1;
for pHIV = vecpHIV
for i=1:n
    disp([i,pHIV])
    tic
    [paramTab,mu,vecAlphas] = sampleParameters_v3(true,true,true,true,b);   %Ct,Ng,HIV,syph
    %
    paramTab{3}.modelType='SICTP';
    paramTab{3}.p = pHIV; paramTab{3}.eta = 4; paramTab{3}.zeta=randPERT(46,60,71,1)/100;
    paramTab{3}.alpha0 = paramTab{3}.alpha;
    paramTab{3}.mu = mu;
    [paramTab{3}.RSICTP,~,paramTab{3}.alpha] = Rp_SICTP(paramTab{3}.betaI,paramTab{3}.betaC,...
        paramTab{3}.theta,paramTab{3}.sigma,paramTab{3}.zeta,paramTab{3}.eta,...
        paramTab{3}.p,paramTab{3}.mu,b,paramTab{3}.rhob);
    disp(paramTab{3}.RSICTP)
    
    paramRho.rho_h = paramTab{3}.rhob;
    paramRho.rho_s = paramTab{4}.rhob;
    paramRho.rho_c = paramTab{1}.rhob;
    paramRho.rho_g = paramTab{2}.rhob;
    
    paramRho.rho_hs = 0;
    paramRho.rho_hc = 0;
    paramRho.rho_hg = 0;
    paramRho.rho_sc = 0;
    paramRho.rho_sg = 0;
    paramRho.rho_cg = 0;
    paramRho.rho_hsc = 0;
    paramRho.rho_hsg = 0;
    paramRho.rho_hcg = 0;
    paramRho.rho_scg = 0;
    paramRho.rho_hscg = 0;
    paramRho.eta_s_prep = 1;
    paramRho.eta_c_prep = 4;
    paramRho.eta_g_prep = 4;
    paramRho.eta_s_art = 0;
    paramRho.eta_c_art = 0;
    paramRho.eta_g_art = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    
    c=0;f=1;
    %paramTab{2}.beta= 13.035779548202
    [~,~,~,P_inf] = U1234_SICTPSEIIISSEIIS2_v3(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,c,f,'ode45');
    
    P_inf_tot(i,1)=P_inf(1);
    P_inf_tot(i,2)=P_inf(2);
    P_inf_tot(i,3)=P_inf(3);
    P_inf_tot(i,4)=P_inf(4);
    toc    
end

tabRecap(k,1) = pHIV;
tabRecap(k,2) = paramRho.eta_s_art;
tabRecap(k,3) = paramRho.eta_c_art;
tabRecap(k,4) = paramRho.eta_g_art;
tabRecap(k,5) = sum(P_inf_tot(:,1)>0.005)/N;
tabRecap(k,6) = sum(P_inf_tot(:,2)>0.005)/N;
tabRecap(k,7) = sum(P_inf_tot(:,3)>0.005)/N;
tabRecap(k,8) = sum(P_inf_tot(:,4)>0.005)/N;
k=k+1;
end
load handel
sound(y,Fs)
%pour resoudre le probleme d'optim selon une variable rho particuliere :...
%faire une fonction auxiliaire qui prend en entree cette variable et qui
%appelle le system d'ODE
tabRecap





%%
clear all;close all;
vecPCt=[]; vecPNg=[]; vecPHIV=[]; vecPS=[]; vecPprevagay=[];
b=2; 
for i=1:500
    [paramTab,mu,vecAlphas] = sampleParameters_v3(true,true,true,true,b);
    vecPCt(i)= paramTab{1}.P;
    vecPNg(i)= paramTab{2}.P;
    vecPHIV(i)= paramTab{3}.PSICT;
    vecPS(i)= paramTab{4}.P;
    vecPprevagay(i)= paramTab{3}.Pprev;
end


close all;
fig = figure()
BinWidth=0.01%(max(alphaTot(:))-min(alphaTot(:)))/60;
histogram(vecPHIV,'FaceColor','red','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %HIV
hold on
histogram(vecPS,'FaceColor','yellow','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %syphilis
histogram(vecPCt,'FaceColor','blue','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ct
histogram(vecPNg,'FaceColor','green','EdgeColor','none','Normalization','probability','BinWidth',BinWidth); %Ng
legend('$\Pi_h(\rho_0^h))$ (prevenir)',...
    '$\Pi_s(\rho_0^s)$ (adapt. prevenir)',...
    '$\Pi_c(\rho_0^c)$ (adapt. prevenir)',...
    '$\Pi_g(\rho_0^g)$ (adapt. prevenir)',...
    'Interpreter','latex','FontSize',12,'Box','off')


mean(vecPHIV)/mean(vecPprevagay)
%--------------------------------------------------------------------------%
