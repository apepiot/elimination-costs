% SIRxSIS^2 model - combined testing strategy

clear all

%%
condition = 0;
while ~condition
    [beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);
    gamma3=gamma1+gamma2/2;
    beta3=(beta1+beta2)/3;
    s3=0.9;
    
    
    R1=beta1/(gamma1+mu);
    R2=beta2/(gamma2+mu);
    R3=beta3/(gamma3+mu);
    
    condition = (R1>1 & R2>1 & R3>1);
end



%% Verification de l'endemic equilibria

% Populations initiales
S0      = b/mu*1.1;
I10     = b/mu*0.1;
I20     = b/mu*0.1;
I30     = b/mu*0.01;    %I_12 (coinfection)
I120    = 0; I230 = 0; I130 = 0;I1230 = 0.;
R10 = 0;IR2=0;IR3=0;IR23=0;

% Parametres du systeme d'ODE
tspan = 0:0.1:10000;
Y0 = [S0; I10; I20; I30; I120; I230; I130; I1230;R10;IR2;IR3;IR23];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) SIRSIS2(t,Y,b,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,rho,mu,'frequency'),tspan,Y0, options);
T = Ys(end,:)
Nnum = sum(T);

%prevalences
P1num = (T(2)+T(5)+T(7)+T(8))/Nnum;
P2num = (T(3)+T(5)+T(6)+T(8)+T(10)+T(12))/Nnum;
P3num = (T(4)+T(6)+T(7)+T(8)+T(11)+T(12))/Nnum;
P12num = 1-(T(1)+T(4)+sum(T(9:12)))/Nnum;
P23num = 1-(T(1)+T(2)-T(9))/Nnum;
P13num = 1-(T(1)+T(3)+T(9)+T(10))/Nnum;
P123num = 1 - (T(1)+T(9))/Nnum;


gamma1p=gamma1+s1*rho;gamma2p=gamma2+s2*rho;gamma3p=gamma3+s3*rho;
R1p=beta1/(gamma1p+mu);
R2p=beta2/(gamma2p+mu);
R3p=beta3/(gamma3p+mu);

P1th = mu/beta1*(R1p-1);%ok
P2th = 1-1/R2p;%ok
P3th = 1-1/R3p;%ok

% Equilibre (theorique)
lambda1=mu*(R1p-1);lambda2=beta2*(1-1/R2p);lambda3=beta3*(1-1/R3p);
gamma3t2 = gamma3p-s3*s2*rho;
gamma1t2 = gamma1p-s1*s2*rho;
gamma2t3 = gamma2p-s3*s2*rho;
gamma2t1 = gamma2p-s1*s2*rho;
gamma3t1 = gamma3p-s1*s3*rho;
gamma1t3 = gamma1p-s1*s3*rho;
gamma23t = s2*s3*rho;
gamma12t= s1*s2*rho;
gamma23 = gamma2p+gamma3p-gamma23t;
gamma12 = gamma1p+gamma2p-gamma12t;
gamma123= gamma1+gamma2+gamma3+(1-(1-s1)*(1-s2)*(1-s3))*rho;
gamma3t= gamma3+s3*(1-s1)*(1-s2)*rho;

S12 = b/mu/R1p*(lambda1+gamma2p+mu)/(lambda1+beta2);
S13 = b/mu/R1p*(lambda1+gamma3p+mu)/(lambda1+beta3);
I1_3 = b/mu*lambda1/(beta2+gamma1t2)*(gamma2t1/beta1+mu/b*S12);%ok
I1_2 = b/mu*lambda1/(beta3+gamma1t3)*(gamma3t1/beta1+mu/b*S13);%ok
I2_3 = b/(mu*R1p)-S12;%ok
I12_3= b*lambda1/(mu*beta1)-I1_3;%ok

Sth  = (b+gamma3t2*S12+gamma2t3*S13+gamma23t*b/mu/R1p)/(lambda1+lambda2+lambda3+gamma23+mu);
I2th = S13-Sth;
I3th = S12-Sth;
I23th= b/(mu*R1p)-S12-S13+Sth;
I1th =I1_2+(-lambda2*I1_2-gamma3t*I12_3-lambda1*S13+lambda1*Sth)/(lambda2+lambda3+gamma12+mu+gamma3t);
I12th=(lambda2*I1_2+gamma3t*I12_3+lambda1*S13-lambda1*Sth)/(lambda2+lambda3+gamma12+mu+gamma3t);
I13th=0;
IR2th=0;
IR3th=0;
R1th=0;
I123th=0%(lambda2*I1_3+lambda1*I2_3-(gamma12+mu)*I12_3)/(gamma123+gamma12-gamma3t);
IR23th=0;

[Sth,I1th,I2th,I3th,I12th,I23th,I13th,I123th,R1th,IR2th,IR3th,IR23th]
T

% prevalence
S23   = b/mu*(gamma2t3/R3p+gamma3t2/R2p+mu+gamma23t)/(beta2+beta3-mu-gamma23t); %T(1)+T(2)
P123th = 1-mu/b*(S23-I1th);
P123num = (sum(T)-T(1)-T(9))/sum(T);%ok


% %% I123+IR23 doit etre egale a I123 dans le modele SIS^3 si gamma1=0,s1=0? (ça ne marche pas)
% I23_SISIS2 = T(8)+T(12);%I123+IR23
% Y0 = [S0; I10; I20; I30; I120; I230; I130; I1230];
% [ts,Ys] = ode45(@(t,Y) SIS3(t,Y,b,beta1,beta2,beta3,0,gamma2,gamma3,0,s2,s3,rho,mu,'frequency'),tspan,Y0, options);
% Ys(end,:)
% Y0 = [S0; I10; I20; I30; I120; I230; I130; I1230;R10;IR2;IR3;IR23];
% [ts,Ys] = ode45(@(t,Y) SIRSIS2(t,Y,b,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,rho,mu,'frequency'),tspan,Y0, options);
% Tt = Ys(end,:)
% [Tt(2)+Tt(9),Tt(7)+Tt(11),Tt(5)+Tt(10),Tt(8)+Tt(12)]
%
% %calcul de I123 SIS^3 si gamma1=0,s1=0
% lambda10 = mu*(beta1/mu-1);
% gamma23t=s2*s3*rho;
% gamma2t3 = gamma2p-s3*s2*rho;
% gamma3t2 = gamma3p-s3*s2*rho;
% lambda2=beta2*(1-1/R2p);lambda3=beta3*(1-1/R3p);
% S23   = b/mu*(gamma2t3/R3p+gamma3t2/R2p+mu+gamma23t)/(beta2+beta3-mu-gamma23t); %T(1)+T(2)
% I23_1 = b/mu*(1-1/R2p-1/R3p)+S23; %I23+I123, T(6)+T(8)
% R1=beta1/mu;R1p=R1; gamma2t1=gamma2p;gamma12t=0;gamma1t2=0;
% S12   = b/mu*(gamma1t2/R2p+gamma2t1/R1p+mu+gamma12t)/(beta1+beta2-mu-gamma12t); %T(1)+T(4)
% I12_3 = b/mu*(1-1/R1p-1/R2p)+S12; %I12+I123 T(5)+T(8)
% gamma1t3=0;gamma3t1=gamma3p;gamma13t=0;
% S13   = b/mu*(gamma1t3/R3p+gamma3t1/R1p+mu+gamma13t)/(beta1+beta3-mu-gamma13t); %T(1)+T(3)
% I13_2 = b/mu*(1-1/R1p-1/R3p)+S13; %I13+I123
% gamma123=gamma2+gamma3+(1-(1-s2)*(1-s3))*rho;
% I123_SIS3 = (lambda10*I23_1+lambda2*I13_2+lambda3*I12_3)/(lambda10+lambda2+lambda3+gamma123+mu)

%% Utility function
alpha1=beta1/s1*(1-1/R1);
alpha2=beta2/s2*(1-1/R2);
alpha3=beta3/s3*(1-1/R3);
maxalpha=max([alpha1,alpha2,alpha3]);
rho=0:maxalpha/1000:1.3*maxalpha;
c=-0.5;

[U,P] = U_SIRSIS2(rho,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c);
[U123,P123] = U123_SIRSIS2(rho,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c);
[U12,P12] = U_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
[U13,P13] = U_SIRSIS7(rho,beta1,beta3,gamma1,gamma3,s1,s3,b,mu,c);
[U23,P23] = U_SIS2(rho,beta2,beta3,gamma2,gamma3,s2,s3,b,mu,c);

gamma1p=gamma1+s1*rho;
gamma2p=gamma2+s2*rho;
gamma3p=gamma3+s3*rho;
R1p = beta1./(gamma1p+mu);
R2p = beta2./(gamma2p+mu);
R3p = beta3./(gamma3p+mu);

U1 = max(rho.*(mu/beta1*(R1p-1)-c),0);
U2 = max(rho.*(1-1./R2p-c),0);
U3 = max(rho.*(1-1./R3p-c),0);

c11 = -mu/beta1*(1-1/R1);
c01 = mu/beta1*(1-R1);
c22 = 1/R2-1;
c02 = 1-1/R2;
c33 = 1/R3-1;
c03 = 1-1/R3;
[U120,c012] = U12_SIRSIS7(0,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,0); %c120 = P12(0)
[U130,c013] = U12_SIRSIS7(0,beta1,beta3,gamma1,gamma3,s1,s3,b,mu,0); %c130 = P13(0)
[U230,c023] = U12_SIS2(0,beta3,beta2,gamma3,gamma2,s3,s2,b,mu,0); %need to be checked
[U120,P120,c112] = U12_SIRSIS7(alpha1,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,0);
[U120,P120,c212] = U12_SIRSIS7(alpha2,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,0);
[U130,P130,c113] = U12_SIRSIS7(alpha1,beta1,beta3,gamma1,gamma3,s1,s3,b,mu,0);
[U130,P130,c313] = U12_SIRSIS7(alpha3,beta1,beta3,gamma1,gamma3,s1,s3,b,mu,0);
[U230,P230,c223] = U12_SIS2(alpha2,beta2,beta3,gamma2,gamma3,s2,s3,b,mu,0);
[U230,P230,c323] = U12_SIS2(alpha3,beta2,beta3,gamma2,gamma3,s2,s3,b,mu,0);
[U1230,c0123,dU123] = U123_SIRSIS2(0,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,0)
[U1230,P1230,c1123] = U123_SIRSIS2(alpha1,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,0)
[U1230,P1230,c2123] = U123_SIRSIS2(alpha2,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,0)
[U1230,P1230,c3123] = U123_SIRSIS2(alpha3,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,0)

clear U120 U130 U230 P120 P130 P230 P1230 U1230
%% U's as a function of rho
figure
plot(rho,U1,':',rho,U2,':',rho,U3,':')
hold on
plot(rho,U12,'--',rho,U23,'--',rho,U13,'--',rho,U123)
plot(rho,U,'LineWidth',1.5)
legend('U_1','U_2','U_3','U_{12}','U_{23}','U_{13}','U_{123}','U')

title([{'Utility U of the SIR$\times$SIS$^2$ model'},...
    {['$\beta_1=$',num2str(round(beta1,2)), ' $\beta_2=$',num2str(round(beta2,2)),' $\beta_3=$',num2str(round(beta2,2)),' $\gamma_1=$',num2str(round(gamma1,2)),' $\gamma_2=$',...
    num2str(round(gamma2,2)),' $\gamma_3=$',num2str(round(gamma3,2)), ' $s_1=$', num2str(s1),' $s_2=$', num2str(s2),' $s_3=$', num2str(s3), ' $\mu=$', num2str(mu)]}, ...
    {[' $\mathtt R_1(0)=$' num2str(round(R1,2)), ' $\mathtt R_2(0)=$' num2str(round(R2,2)),' $\mathtt R_3(0)=$' num2str(round(R3,2)),...
    ' $\rho_1\prime=$',num2str(round(alpha1,2)),' $\rho_2\prime=$',num2str(round(alpha2,2)),' $\rho_3\prime=$',num2str(round(alpha3,2)),' $c=$',num2str(round(c,2))]},...
    {[' $c_0^1=$',num2str(round(c01,2)),' $c_1^1=$',num2str(round(c11,2)),' $c_0^2=$',num2str(round(c02,2)),' $c_2^2=$',num2str(round(c22,2)),...
    ' $c_0^3=$',num2str(round(c03,2)),' $c_3^3=$',num2str(round(c33,2)),...
    ' $c_0^{12}=$',num2str(round(c012,2)),' $c_1^{12}=$',num2str(round(c112,2)), ' $c_2^{12}=$',num2str(round(c212,2)),...
    ' $c_0^{23}=$',num2str(round(c023,2))]},{[' $c_2^{23}=$',num2str(round(c223,2)), ' $c_3^{23}=$',num2str(round(c323,2)),...
    ' $c_0^{13}=$',num2str(round(c013,2)),' $c_1^{13}=$',num2str(round(c113,2)), ' $c_3^{13}=$',num2str(round(c313,2)),...
    ' $c_0^{123}=$',num2str(round(c0123,2)),' $c_1^{123}=$',num2str(round(c1123,2)), ' $c_2^{123}=$',num2str(round(c2123,2)),' $c_3^{123}=$',num2str(round(c3123,2))]}],...
    'Interpreter','latex')
xlabel("$\rho$","fontweight","bold",'Interpreter','latex')
ylabel("$U(\rho)$","fontweight","bold",'Interpreter','latex')



%% comparaison entre les deux façons de trouver rhohat
figure
plot(vecC,vecRhohat, vecC, vecRhotothat)

%% beta in fonction of gamma^{-1} for a certain given set of parameters (SIR,SIS models)
gamma0 = 6:-0.01:1/3;
gamma1=gamma0;
mu = 1/35;
P1 = 16.1/100*9/100;
beta1 = mu./(mu./(gamma1+mu)-P1);
R1 = beta1./(gamma1+mu);
s1=1;
alpha1 = beta1/s1.*(1-1./R1).*(R1>=1);
plot(1./gamma1(R1>=1),alpha1(R1>=1))

hold on;
gamma2=gamma0;
mu = 1/35;
P2 = 5.3/100;
R2 = 1./(1-P2);
beta2 = R2.*(gamma2+mu);
s2=1;
alpha2 = beta2/s2.*(1-1./R2).*(R2>=1);
hold on
plot(1./gamma2,alpha2)

gamma3=gamma0;
mu = 1/35;
P3 = 6.6/100;
R3 = 1./(1-P3);
beta3 = R3.*(gamma3+mu);
s3=1;
alpha3 = beta3/s3.*(1-1./R3).*(R3>=1);
hold on
plot(1./gamma3,alpha3)
plot(0:3/100:3,0.0814*ones(1,101), 'k--')

title(['$\Pi_{HIV}=$',num2str(round(P1,3)),' $\Pi_{ch}=$',num2str(P2),' $\Pi_{s}=$',num2str(P3)],'Interpreter','latex')
ylim([0 1])
xlabel('Time between infection and diagnosis $\gamma(0)^{-1}$ in years','Interpreter','latex')
ylabel('$\rho\prime$','Interpreter','latex')
legend('$\rho_{HIV}\prime$', '$\rho_{ch}\prime$','$\rho_{s}\prime$', 'Interpreter', 'latex')

%plot R
figure
plot(1./gamma1(R1>1),R1(R1>1),1./gamma2,R2*ones(length(gamma2),1),1./gamma3,R3*ones(length(gamma2),1))
ylim([0.90;2])
xlabel('Time between infection and diagnosis $\gamma(0)^{-1}$ in years','Interpreter','latex')
ylabel('$\mathtt R(0)$','Interpreter','latex')
legend('$\mathtt R_{HIV}(0)$', '$\mathtt R_{ch}(0)$','$\mathtt R_{s}(0)$', 'Interpreter','latex')


%% PLOT RHOHAT
script_SIRSIS2_plot;