%% Utility septembre 2024
clear all; close all;
%--------------------------
vecC = -1:0.1:1;
b = 1000; mu = 1/30.6;

gamma1 = 1/2.9;
s1     = 1;
rho1b  = 0.12;

gamma2 = 365/14*0.11;
s2 = 0.85; %ct
rho2b = 0.12;

P1b = 16.1/100; R1b = 1/(1-P1b);
beta1 = R1b*(gamma1+s1*rho1b+mu);

P2b = 5.3/100; R2b = 1/(1-P2b);
beta2 = R2b*(gamma2+s2*rho2b+mu);

R10 = beta1/(gamma1 + mu);
R20 = beta2/(gamma2 + mu);
alpha1 = beta1*(1-1/R10)/s1;
alpha2 = beta2*(1-1/R20)/s2;
vecRho = sort([alpha1,alpha2,linspace(0,max(alpha1,alpha2)*1.5,200)]); %1:hiv, 2:ct
%--------------------------

i=1;P12=[];P1=[];P2=[];
Nequ = b/mu;
for rho=vecRho
    gamma2p = gamma2 + s2*rho;
	gamma1p = gamma1 + s1*rho;
    gamma12 = gamma1+gamma2+(1-(1-s1)*(1-s2))*rho;
    gamma1t = gamma1p-s1*s2*rho;
    gamma2t = gamma2p-s1*s2*rho;
    
    R1p = beta1/(gamma1p + mu);
	R2p = beta2/(gamma2p + mu);
    
    Lambda1p = mu*(R1p-1)/beta1;
    Lambda2p = beta2*(1-1/R2p);
    
    
    %SES12   = b/mu*(mu*(R1p-1) + gamma2p+mu)/(mu*(R1p-1)+beta2)/R1p;%ok
    SES12   =  b/(mu*R1p)*(Lambda1p+gamma2p+mu)/(Lambda1p+beta2)
    I2ES12  = b/(mu*R1p)-SES12;
    I12ES12 = b/mu*Lambda1p/R1p*((Lambda2p+gamma1p+mu)/(gamma1p+mu) - (Lambda1p+mu+gamma2p)/(Lambda1p+beta2))/(beta2+gamma1p-s1*s2*rho);
    I1ES12  = b/mu*Lambda1p/R1p*(gamma2t/(gamma1p+mu)+(Lambda1p+gamma2p+mu)/(Lambda1p+beta2))/(Lambda2p+gamma1p+mu+gamma2t);
    IR2ES12 = b/mu*(1-1/R2p)-I2ES12-I12ES12;
    R1ES12  = b/(mu*R2p)-SES12 -I1ES12;
    N       = b/mu;
    NVT    = (N-R1ES12);
    P12(i)  = max((I2ES12+I12ES12+I1ES12+IR2ES12+R1ES12)/N,0);
    
    SES1 = b/(mu*R1p); SES12+I2ES12;
    IES1 = b/beta1*(R1p-1); I12ES12+I1ES12; %pas ok
    R1ES = gamma1p*b/(mu*beta1)*(R1p-1); IR2ES12+R1ES12 %pas ok
    
    (SES1+IES1+R1ES)-b/mu
    
    %P1(i)   = max((I12ES12+I1ES12)/(NVT-IR2ES12),0);
    P1(i)  = max(IES1/(SES1+IES1),0);
    P2(i)   = max(1-1/R2p,0);

    i=i+1;
end

P12 = P12.*(vecRho>alpha1 & vecRho>alpha2) + P1.*(vecRho<=alpha2 & vecRho>alpha1) + P2.*(vecRho<=alpha1 & vecRho>alpha2); 

for c=vecC
    % utility
    U12 = vecRho.*(P12-c);
    U1 = vecRho.*(P1-c);
    U2 = vecRho.*(P2-c);
    figure(1)
    hold on
    %plot(vecRho, U12,'DisplayName',['c=',num2str(c)])
    %plot(vecRho, U12,'DisplayName','U_{12}')
    %plot(vecRho, U2, 'DisplayName','U_{2}')
    plot(vecRho,U1, 'DisplayName','U_{1}')
end
legend()
%legend('U_1','U_2','U_{12}')

[aa,bb] = max(U12);
vecRho(bb)

tab = [vecRho',P1',P2',P12']

%%
close all

addpath('C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\_archives\')
vecC=sort([-0.041,-0.0347,-0.032,-0.0315,-0.031,linspace(-0.15,0.4,100)]); %sort([0.135,0.132,0.133,0.137,-0.051,-0.052,-0.055,-0.25:0.01:0,0.001:0.0001:0.003,0.004:0.001:0.005,0.01:0.01:0.4]);
rho12=[]; rho2=[]; i=1; rhohat =[]; rho1=[]; rho12th=[]; up_bnd=[];
for c=vecC
    CSIS2 = @(rho) -U12_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    rho12(i) = min([max(fmincon(CSIS2,0),0),max(alpha1,alpha2)]);

    rho2(i)  = min(max(beta2/(2*s2)*(1-1/R20-c),0),alpha2);
    rhohat1_th = beta1*(sqrt(R10*mu/(mu+beta1*c))-1)/(R10*s1);
    rho1(i)  = min(max(rhohat1_th,0),alpha1);
    
    U_12 = U_SIRSIS7(rho12(i),beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    U_2  = U_SIRSIS7(rho2(i),beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    
    if U_12<=U_2
        rhohat(i) = rho2(i);
    else 
        rhohat(i) = rho12(i);
    end
    
    i=i+1;
end

figure()
plot(vecC,min(rhohat,alpha2))
hold on
plot(vecC,rho1)
plot(vecC,min(rho2,alpha2))
legend('\rho_{12}','\rho_1','\rho_2')

tab = [vecC',rho12',rho1',rho2',rhohat',rho12th',up_bnd'];
tab=round(tab,5)
format short


