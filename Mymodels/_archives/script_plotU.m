%This script plots the utility/cost function for a n-disease model (n=1,..,4)
clear all; close all;

%Utility or Cost ?
sgnU = -1;
c       = 0;
%Which diseases ?
hiv=true;syph=true;ct=true;ng=true;

% Parameters
mu=0.0328;b=2;

%Ct
paramCt.beta  = 1.1344%0.9698;
paramCt.sigma = 24.5225%34.0030;
paramCt.gamma = 27.6773%28.4550;
paramCt.nu    = 0.9822%0.8282;
paramCt.eps   = 0.0344%0.0419;

%Ng
paramNg.beta  = 7.8105%'NaN';
paramNg.sigma = 26.0714%'NaN';
paramNg.gamma = 29.8610%'NaN';
paramNg.nu    = 1.5000%'NaN';
paramNg.eps   = 0.8303%'NaN';

%HIV
paramHIV.betaI = 1.4201%1.4286; 
paramHIV.betaC = 0.1561%0.1575; 
paramHIV.gamma = 0 ;
paramHIV.theta = 0.1270%0.1425;
paramHIV.sigma = 6.9468%5.7693; 

%Syphilis
paramS.beta   = 0.0881%'NaN'; 
paramS.sigma  = 4.3483%'NaN'; 
paramS.tau    = 9.2706%'NaN'; 
paramS.gamma1 = 0%'NaN'; 
paramS.theta  = 3.0841%'NaN'; 
paramS.gamma3 = 0.0510%'NaN';
paramS.nu=0;


%%
N=hiv+ct+ng+syph;
paramCt.R='NaN';paramCt.alpha='NaN';
paramNg.R='NaN';paramNg.alpha='NaN';
paramHIV.R='NaN';paramHIV.alpha='NaN';
paramS.R='NaN';paramS.alpha='NaN';

if ct
    [paramCt.R,~,paramCt.alpha]   = Rp_SEIIS_v4(paramCt.beta,paramCt.nu,...
        paramCt.eps,paramCt.sigma,paramCt.gamma,mu,b,0);
end
if ng
    [paramNg.R,~,paramNg.alpha]   = Rp_SEIIS_v4(paramNg.beta,paramNg.nu,...
        paramNg.eps,paramNg.sigma,paramNg.gamma,mu,b,0);
end
if hiv 
    [paramHIV.R,~,paramHIV.alpha] = Rp_SICR_v4(paramHIV.betaI,paramHIV.betaC,...
        paramHIV.theta,paramHIV.sigma,paramHIV.gamma,mu,b,0);
end
if syph
    [paramS.R,~,paramS.alpha] = Rp_SEIIIS_v4(paramS.beta,paramS.sigma,...
        paramS.tau,paramS.nu,paramS.gamma1,paramS.theta,paramS.gamma3,mu,b,0);
end

paramCt.disease  = 'Ct'; paramCt.modelType = 'SEIIS';
paramNg.disease  = 'Ng'; paramNg.modelType = 'SEIIS';
paramHIV.disease = 'HIV'; paramHIV.modelType = 'SICR';
paramS.disease   = 'Syph.'; paramS.modelType = 'SEIIIS';


%paramTabAll{1} = paramCt;
%paramTabAll{2} = paramNg;
%paramTabAll{3} = paramHIV;
%paramTabAll{4} = paramS;
%paramTab = paramTabAll([Ct,Ng,HIV,Syph]);


%% Manually: SICTxSEIIS (HIVxCt)
close all;
vecAlphas = [paramCt.alpha,paramHIV.alpha];
maxalpha = max(vecAlphas);
minalpha = min(vecAlphas);

vecRHO = 0:maxalpha/1001:(maxalpha*1.5);
Uc     = U_SEIIS_v4(paramCt,mu,b,vecRHO,c,1);
Uh     = U_SICR_v4(paramHIV,mu,b,vecRHO,c,1);
Uhc    = U_SEIISSICR_v4(paramCt,paramHIV,mu,b,vecRHO,c,1);

plot(0,0,'Color','black','LineStyle','-','LineWidth',2,'DisplayName','C')
hold on
plot(vecRHO,sgnU*Uhc,'Color',[0 0.4470 0.7410],'LineStyle','-','DisplayName','C_{hc}','LineWidth',2)
%hold on
plot(vecRHO,sgnU*Uh,'Color',[0.8500 0.3250 0.0980],'LineStyle','--',...
    'LineWidth',2,'DisplayName','C_h')
plot(vecRHO(vecRHO<=minalpha),sgnU*Uc(vecRHO<=minalpha),...
    'Color',[0.9290 0.6940 0.1250],'LineStyle','--','LineWidth',2,'HandleVisibility','off')
plot(vecRHO(vecRHO>=minalpha),sgnU*Uc(vecRHO>=minalpha),...
    'Color',[0.9290 0.6940 0.1250],'LineStyle','-','LineWidth',2,...
    'DisplayName','C_c')


legend('Location','southeast')
xlim([0,vecRHO(end)])
xlabel('\rho')
ylabel('C(\rho)')


%% Manually: SICTxSEIIISxSEIIS^2
close all;
vecAlphas = [paramCt.alpha,paramNg.alpha,paramHIV.alpha,paramS.alpha];
minalpha = min(vecAlphas);
maxalpha = max(vecAlphas);
vecRHO = 0:maxalpha/1001:(maxalpha*1.35);

Uc     = U_SEIIS_v4(paramCt,mu,b,vecRHO,c,1);
Uh     = U_SICR_v4(paramHIV,mu,b,vecRHO,c,1);
Ug     = U_SEIIS_v4(paramNg,mu,b,vecRHO,c,1);
Us     = U_SEIIIS_v4(paramS,mu,b,vecRHO,c,1);
Uhc    = U_SEIISSICR_v4(paramCt,paramHIV,mu,b,vecRHO,c,1);
Uhg    = U_SEIISSICR_v4(paramNg,paramHIV,mu,b,vecRHO,c,1);
Uhs    = U_SICRSEIIIS_v4(paramHIV,paramS,mu,b,vecRHO,c,1);
Usc    = U_SEIISSEIIIS_v4(paramCt,paramS,mu,b,vecRHO,c,1);
Ucg    = U_SEIIS2_v4(paramCt,paramNg,mu,b,vecRHO,c,1);
Usg    = U_SEIISSEIIIS_v4(paramNg,paramS,mu,b,vecRHO,c,1);
Uhcg   = U_SEIIS2SICR_v4(paramCt,paramNg,paramHIV,mu,b,vecRHO,c,1);
Uhsc   = U_SEIISSICRSEIIIS_v4(paramCt,paramHIV,paramS,mu,b,vecRHO,c,1);
Uhsg   = U_SEIISSICRSEIIIS_v4(paramNg,paramHIV,paramS,mu,b,vecRHO,c,1);
Uscg   = U_SEIIS2SEIIIS_v4(paramCt,paramNg,paramS,mu,b,vecRHO,c,1);
Uhscg  = U_SEIIS2SICRSEIIIS_v4(paramCt,paramNg,paramHIV,paramS,...
                                    mu,b,vecRHO,c,1);
%% 
close all

plot(0,0,'Color','black','LineStyle','-','LineWidth',2,'DisplayName','C')
hold on
plot(vecRHO,sgnU*Uhscg,'Color',[0 0.4470 0.7410],'LineStyle','-','DisplayName','C_{hscg}','LineWidth',2)
plot(vecRHO(vecRHO>=vecAlphas(4)),sgnU*Uhcg(vecRHO>=vecAlphas(4)),'Color',[0.8500 0.3250 0.0980],'LineStyle','-',...
    'LineWidth',2,'DisplayName','C_{hcg}')
plot(vecRHO(vecRHO<vecAlphas(4)),sgnU*Uhcg(vecRHO<vecAlphas(4)),'Color',[0.8500 0.3250 0.0980],'LineStyle','--',...
    'LineWidth',2,'DisplayName','C_{hcg}','HandleVisibility','off')
plot(vecRHO(vecRHO>=vecAlphas(3)),sgnU*Ucg(vecRHO>=vecAlphas(3)),'Color',[0.4660 0.6740 0.1880],'LineStyle','-',...
    'LineWidth',2,'DisplayName','C_{cg}')
plot(vecRHO(vecRHO<vecAlphas(3)),sgnU*Ucg(vecRHO<vecAlphas(3)),'Color',[0.4660 0.6740 0.1880],'LineStyle','--',...
    'LineWidth',2,'DisplayName','C_{cg}','HandleVisibility','off')
plot(vecRHO(vecRHO>=vecAlphas(1)),sgnU*Ug(vecRHO>=vecAlphas(1)),'Color',[0.9290 0.6940 0.1250],'LineStyle','-',...
    'LineWidth',2,'DisplayName','C_{g}')
plot(vecRHO(vecRHO<vecAlphas(1)),sgnU*Ug(vecRHO<vecAlphas(1)),'Color',[0.9290 0.6940 0.1250],'LineStyle','--',...
    'LineWidth',2,'DisplayName','C_{g}','HandleVisibility','off')

legend('Location','southeast')
xlim([0,vecRHO(end)])
%ylim([sgnU*1.2*max(Uhscg),0])
xlabel('\rho')
ylabel('C(\rho)')