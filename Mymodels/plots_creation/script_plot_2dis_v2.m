% Script to plot the solution of the two-disease game
% essentially for the second paper, version sept. 2023
% Difference with script_plot_2dis: SICTP model instead of SICT
clear all; close all
%select the infections of the model
Ct=true; Ng=false; HIV=true; Syph=false; plotFigure=true;
biasFactor = 1;
newsetofParameters = true; 
%if false, choose a set:
fileSet = "set4.m";


%nbSim = 1;
nSEIIS=Ct+Ng; nSICR=HIV; nSEIIIS=Syph;
N = nSICR+nSEIIS+nSEIIIS;
if newsetofParameters
    [paramTab,mu,vecAlphas] = sampleParameters_v3(Ct,Ng,HIV,Syph,1);   %Ct,Ng,HIV,syph  
    
    iHIV = Ct+Ng+HIV;
    paramTab{iHIV}.RSICTP=0;
    while HIV && (paramTab{iHIV}.RSICTP<=1 || paramTab{iHIV}.alpha>1.5)
        [paramHIV,~,~] = sampleParameters_v3(false,false,HIV,false,1);
        paramTab{iHIV} = paramHIV{:};
        paramTab{iHIV}.modelType='SICTP';

        paramTab{iHIV}.p = 0.5;
        paramTab{iHIV}.eta = 4;
        paramTab{iHIV}.mu = mu;
        paramTab{iHIV}.zeta=randPERT(46,60,71,1)/100;

        [paramTab{iHIV}.RSICTP,~,paramTab{iHIV}.alpha] = Rp_SICTP(paramTab{iHIV}.betaI,paramTab{iHIV}.betaC,...
                                paramTab{iHIV}.theta,paramTab{iHIV}.sigma,paramTab{iHIV}.zeta,paramTab{iHIV}.eta,...
                                paramTab{iHIV}.p,paramTab{iHIV}.mu,1,paramTab{iHIV}.rhob);
        vecAlphas(Ct+Ng+HIV) = paramTab{iHIV}.alpha;

    end 
else 
    run(fileSet)
end
%%
%paramCt = paramTab{Ct}; paramNg = paramTab{(Ct+Ng)*Ng}; 
%paramHIV= paramTab{iHIV*HIV}; paramS  = paramTab{(Ct+Ng+HIV+Syph)*Syph};

cleft=-0.1; cright=0.4;
vecC = [-0.3,-0.23,-0.22,-0.2,-0.18:0.001:-0.005,-0.0044:0.0001:-0.003,-0.002:0.001:0.26,0.27:0.01:0.3,0.5]; 
[tab,tabco,tabcn,tabTimes] = findRhohat_v6(N,paramTab,mu,1,vecC,biasFactor);
cs = findThresholds_v5(Ct+Ng,HIV,Syph,tab,vecAlphas,vecC); 
%%
plot_zones_v4

yticks([0,0.5,1,1.5,2])
yticklabels({'0','0.5','1','',''})

%yticks([ 0 0.01 0.02 vecAlphas(cs.order(1)) 0.04 0.05 0.06 0.07 vecAlphas(cs.order(2)) 0.09 0.1 ])
%yticklabels({'0', '0.01', '0.02', ['$\rho_{',infMini{cs.order(1)},'}^\prime$'],...
%     '0.04', '0.05', '0.06', '0.07', ['$\rho_{',infMini{cs.order(2)},'}^\prime$'],'0.09', '0.1'})

