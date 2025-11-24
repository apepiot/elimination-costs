
%1:'Chlam', 2:'Gono', 3:'HIV', 4:'Syphilis'


%%
clear all;
% P1     = 7.4/100;[6.4,8.4]/100;%randPERT(6.4,7.4,8.4,1)/100;
% vecR1     = 1./(1-P1);
% vecsigma1 = 365/11;%365./[7,21];%365./randPERT(7,11,21,1);
% vecnu1    = 12./12%12./[6,36];%randPERT(6,12,36,1);
% vecgamma1 = 365./14%365./[5,20];%randPERT(5,14,20,1);%365/(31.5-1./sigma1); 
% veceps1   = 0.11%[0,0.5];%randPERT(0,0.11,0.5,1);
% vecMu = 1./30.6%1./[27.2,33.7];%randPERT(27.2,30.6,33.7,1);
%%
P1     = [6,8]./100%randPERT(6,7,8,1)/100;
vecR1     = 1./(1-P1);
vecsigma1 = 365./[1,14];%365/randPERT(1,5,14,1);
vecnu1    = 12./[4,12]%12/randPERT(4,6,12,1);
veceps1   = [0.5,1]%randPERT(0.5,0.8,1,1);
vecgamma1 = 365./[5,20]%365/randPERT(5,14,20,1);%365/3; 
vecMu = 1./30.6

vecBeta=[];
for R1=vecR1
    for sigma1=vecsigma1
        for nu1=vecnu1
            for gamma1=vecgamma1
                for eps1=veceps1
                    for mu=vecMu
                        beta1  = R1/((sigma1*(gamma1*(1-eps1) + mu + nu1))./((mu+sigma1).*(gamma1 + mu + nu1).*(mu + nu1)));
                        vecBeta = [vecBeta,beta1];
                    end
                end
            end
        end
    end
end
%% 
PHIV = randPERT(12.0,14.3,16.9,1)/100;
RSICR     = 1/(1-PHIV);%2.2; %P=1-1/R so R=1/(1-P);
sigmaSICR = 52/randPERT(6.7,8.2,9.8,1); %(1-0.39)*365/(8.2*7)%0.6*(); %8.2 semaines
thetaSICR = 1/randPERT(4,4.4,10,1);%1/9.8;
gammaSICR = 0; %0.39*365/(8.2*7); %0.4*(); %proportion de diag au stade precoce : 40
ratioBeta = randPERT(8.4,9.1,9.6,1);
betaCSICR = RSICR*(sigmaSICR+gammaSICR+mu)*(thetaSICR+mu)/(ratioBeta*(thetaSICR+mu)+sigmaSICR);
betaISICR = ratioBeta*betaCSICR;

%%
clear all
PHIV = [12,16.9]./100;%randPERT(12.0,14.3,16.9,1)/100; %14.3/100;%
vecRSICR     = 1./(1-PHIV);%2.2; %P=1-1/R so R=1/(1-P);
vecsigmaSICR = 52./[6.7,9.8];%randPERT(6.7,8.2,9.8,1); %(1-0.39)*365/(8.2*7)%0.6*(); %8.2 semaines %52/8.2;%
vecthetaSICR = 1./[4,10];%randPERT(4,4.4,10,1);%1/9.8;
vecgammaSICR = 0; %0.39*365/(8.2*7); %0.4*(); %proportion de diag au stade precoce : 40
vecratioBeta = [8.4,9.6];%randPERT(8.4,9.1,9.6,1);%9.1;%
vecMu = 1./[27.2,33.7];%1./30.6;%
vecBetaI=[];vecBetaC=[];
for RSICR=vecRSICR
    for sigmaSICR=vecsigmaSICR
        for thetaSICR=vecthetaSICR
            for gammaSICR=vecgammaSICR
                for mu=vecMu
                    for ratioBeta=vecratioBeta  
                        betaCSICR = RSICR*(sigmaSICR+gammaSICR+mu)*(thetaSICR+mu)/(ratioBeta*(thetaSICR+mu)+sigmaSICR);
                        vecBetaC = [vecBetaC,betaCSICR];
                        vecBetaI = [vecBetaI,ratioBeta*betaCSICR];
                    end
                end
            end
        end
    end
end

min(vecBetaC) 
max(vecBetaC)
min(vecBetaI) 
max(vecBetaI)
%%
clear all

PS      = 7.7/100%[6.7,8.7]/100%randPERT(6.7,7.7,8.7,1)/100;
vecRS      = 1./(1-PS);
vecsigmaS  = 365/25%365./[10,95]%randPERT(10,25,95,1);
vectauS    = 365/46%365./[10,130]%randPERT(10,46,130,1);%(1-0.55*0.31)*365/45; %0.6*
vecthetaS  = 12/3.6%12./[1,12]%randPERT(1,3.6,12,1);
gamma1S = 0%(0.55*0.31)*365/45; %0.2*
vecgamma3S = 1./20%1./[10,30]%randPERT(10,20,30,1);%1/5;
nuS     = 0;
vecBeta=[];
vecMu = 1./30.6% 1./[27.2,33.7]

for RS=vecRS
    for sigmaS=vecsigmaS
        for tauS=vectauS
            for thetaS=vecthetaS
                for gamma3S=vecgamma3S

                    for mu=vecMu
                    vecBeta= [vecBeta, RS./(sigmaS.*((gamma3S+mu+nuS).*(mu+tauS+thetaS)+tauS*thetaS)./((mu + sigmaS).*(mu + thetaS).*(gamma3S + mu + nuS).*(gamma1S + mu + tauS)))];
                    end
                end
            end
        end
    end
end

min(vecBeta)
max(vecBeta)