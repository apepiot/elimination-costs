
%%  Analogy with a SIS model and comparaison of the 2 models %%%%

% In this script, we compare the prevalence of a SISxSIS model and its
% SIS model analogy where beta = beta1+beta2-mu et gamma = gamma1/R02 + gamma2/R01

clear all; close all;

%% Parameters
beta1 = 0.5;
beta2 = 0.5;
gamma1_c = 0.2; s1 = 1; gamma1_v = 0;
gamma2_c = 0.2; s2 = 1; gamma2_v = 0;

gamma1 = gamma1_c + s1*gamma1_v;
gamma2 = gamma2_c + s2*gamma2_v;
e1 = beta1;
e2 = beta2;
p = 0.06;
mu = 0.02;

% Initial populations
N0 = 15;
S0 = 9;
I10 = 2;
I20 = 3;
II0 = N0-S0-I10-I20;

%% Computation of $R_0*
R10 = beta1/(gamma1+mu)
R20 = beta2/(gamma2+mu)


%% *Computing the ODE's system* (SISxSIS model)

MaxTime = 300;
[t, pop] = ode45(@(t,y) SISxSIS(t,y,[beta1 beta2 gamma1 gamma2 p mu e1 e2]),[0 :0.1: MaxTime],[S0 I10 I20 II0]);
S=pop(:,1); I1=pop(:,2); I2=pop(:,3); II=pop(:,4);
popSIS2 = pop;

popSIS2_end = pop(end,:)
tSIS2 = t;

%% Theoretical equilibrium
Nequ = p/mu;
E0 = [p/mu,0,0,0];
E1 = [Nequ/R10, Nequ - Nequ/R10,0,0];
E2 = [Nequ/R20, 0, Nequ - Nequ/R20,0];

S12 = (gamma1/R20 + gamma2/R10 + mu)/(beta1+beta2-mu)*p/mu;
E12 = [S12, Nequ/R20-S12, Nequ/R10 - S12, S12 + Nequ*(1-1/R10 - 1/R20)];


%% Plot the population dynamics
figure(1)
hold on;
plot(t,S,"-r",t,(I1 + I2 + II), 'b')
xlim([0 MaxTime])
xlabel("Time","fontweight","bold")
ylabel("Number","fontweight","bold")
%h = legend("S_{SIS^2}","I_{SIS^2}");
%title("SISxSIS and SIS models")

%title([{'Susceptibles and prevalence of SISxSIS and the SIS associated'},...
 %   {['\beta_1 =  ',num2str(beta1), ' \beta_2 =  ',num2str(beta2),' \gamma_1 =  ',num2str(gamma1),' \gamma_2 =  ',num2str(gamma2), '\mu =', num2str(mu), ' \pi =',  num2str(p), ', R_{10} = ', num2str(round(R10,2)),  ' R_{20} = ', num2str(round(R20,2))]}])


PopEqu_SIS2 = popSIS2(end,:);
sum(popSIS2(end,2:4)) %somme de tous les infectes




%% Seen as a SIS with : 
gamma = gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1;
beta = beta1 + beta2 - mu;

I0 = I10 + I20 + II0;

%% *Computation of $R_0$*
R_0 = beta/(gamma+mu);

%% Computing the ODE's system
%MaxTime = 1000;
[t, pop]= ode45(@(t,y) SIS(t,y,[beta gamma p mu], 'frequency'),[0:0.1:MaxTime],[S0 I0]);
S = pop(:,1); I = pop(:,2);

popSIS = pop;
tSIS = t;
popSIS_end = popSIS(end,:)

%% Plot
%figure(2)

plot(t,S,"--r",t,I, '--b')
xlim([0 MaxTime])
%xlabel("Time","fontweight","bold")
%ylabel("Number","fontweight","bold")
h = legend("S_{SIS^2}","I_{SIS^2}","S_{SIS}","I_{SIS}");
title([{'Susceptibles and prevalence of SISxSIS and the SIS associated'},...
    {['\beta_1 =  ',num2str(beta1), ' \beta_2 =  ',num2str(beta2),' \gamma_1 =  ',num2str(gamma1),' \gamma_2 =  ',num2str(gamma2), '\mu =', num2str(mu), ' \pi =',  num2str(p)]},...
    {[' \gamma = ', num2str(round(gamma,2)), ' \beta = ', num2str(round(beta,2)), ' R_{0}^1 = ', num2str(round(R10,2)),  ' R_{0}^2 = ', num2str(round(R20,2))]}])
    
%title("SIS model")
hold off;

figure(3)
plot(t, popSIS(:,2) - sum(popSIS2(:,2:4)')')
diff = (abs(popSIS(:,2) - sum(popSIS2(:,2:4)')'))./sum(popSIS2(:,2:4)')';

figure(4)
plot(t, diff*100)
title("Relative difference between the prevalences of the SISxSIS model and its SIS approximation")
xlabel('Time')
ylabel('(I_1+I_2+II_{12} - I)/(I_1+I_2+II) * 100')


%% Numerical utility function of the SISxSIS model

N = S+I;
diff(end)

S(end)/N(end)*(beta1*I2(end) + beta2* I1(end) - mu*(N(end) - S(end))) - ( gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1)*(N(end) - S(end))
for i =1:300
    S(i)/N(i)*(beta1*I2(i) + beta2* I1(i) - mu*(N(i)-S(i))) - ( gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1)*(N(i) - S(i))
end



i=10000



%% Evolution de la prevalence a l'equilibie en fonction de gamma_1' pour les modèles SIS et SIS2
close all;

vecGamma1_v = 0:0.01:1.5;
vecGamma1 = gamma1_c + s1*vecGamma1_v;
gamma2 = gamma2_c + s2*gamma2_v;

vecR10 = beta1./(vecGamma1+mu);
R20 = beta2/(gamma2+mu);


% Modele SIS2
[U,PSIS2,vecR10,R20] = utilityFunctionSIS2(beta1,beta2,vecGamma1,gamma2,p,mu,beta1,beta2);

USIS2 = vecGamma1_v.*PSIS2*mu/p;

figure(5)
hold on
plot(vecGamma1_v, PSIS2, 'r')

figure(6) 
hold on
plot(vecGamma1_v, USIS2, 'r')



% Modele SIS
vecGamma = vecGamma1./R20 + gamma2./vecR10;
beta = beta1 + beta2 - mu;
R_0 = beta./(vecGamma+mu);

[USIS,PSIS,Umax] = utilityFunctionSIS(beta,vecGamma, 0, 1,p,mu, 'frequency'); %gamma_cl=0 car deja inclus dans veGamma_v

USIS = vecGamma1_v.*PSIS*mu/p;

figure(5)
hold on
plot(vecGamma1_v, PSIS, '--b')
h = legend("I_{SIS^2}^a","I_{SIS^2}^b","I_{SIS}");
%h = legend("I_{SIS^2}","I_{SIS}");
title([{'Prevalences of SISxSIS and the SIS associated'},...
    {['\beta_1 =  ',num2str(beta1), ' \beta_2 =  ',num2str(beta2),' \gamma_2 =  ',num2str(gamma2), ' \mu =', num2str(mu), ' \pi =',  num2str(p)]},...
    {[' \gamma = ', num2str(round(gamma,2)), ' \beta = ', num2str(round(beta,2)), ' R_{0}^1 = [', num2str(round(min(vecR10),2)),',', num2str(round(max(vecR10),2)), '] R_{0}^2 = ', num2str(round(R20,2))]}])

xlabel("\gamma_1'")
ylabel("I(\gamma_1')")


figure(6)
hold on
plot(vecGamma1_v, USIS, '--b' )
h = legend("U_{SIS^2}^a","U_{SIS^2}^b","U_{SIS}");
%h = legend("U_{SIS^2}","U_{SIS}");
title([{'Utility function of SISxSIS and the SIS associated'},...
    {['\beta_1 =  ',num2str(beta1), ' \beta_2 =  ',num2str(beta2),' \gamma_2 =  ',num2str(gamma2), ' \mu =', num2str(mu), ' \pi =',  num2str(p)]},...
    {[' \gamma = ', num2str(round(gamma,2)), ' \beta = ', num2str(round(beta,2)), ' R_{0}^1 = [', num2str(round(min(vecR10),2)),',', num2str(round(max(vecR10),2)), '] R_{0}^2 = ', num2str(round(R20,2))]}])

xlabel("\gamma_1'")
ylabel("U(\gamma_1')")

%title("SIS model")


Umax



