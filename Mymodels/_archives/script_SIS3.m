% SIS^3 model - combined testing strategy
clear all
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);
gamma3=gamma1+gamma2/2;
beta3=(beta1+beta2)/3;
s3=0.9;

R1=beta1/(gamma1+mu);
R2=beta2/(gamma2+mu);
R3=beta3/(gamma3+mu);

% Populations initiales
S0      = b/mu*1.1;
I10     = b/mu*0.1;
I20     = b/mu*0.1;
I30     = b/mu*0.01;    %I_12 (coinfection)
I120     = 0;
I230 = 0;
I130 = 0;
I1230     = 0.;

% Parametres du systeme d'ODE 
tspan = 0:0.1:100000;
Y0 = [S0; I10; I20; I30; I120; I230; I130; I1230];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) SIS3(t,Y,b,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,rho,mu,'frequency'),tspan,Y0, options);
T = Ys(end,:)
Nnum = sum(T);

P1num = (T(2)+T(5)+T(7)+T(8))/Nnum;
P2num = (T(3)+T(5)+T(6)+T(8))/Nnum;
P3num = (T(4)+T(6)+T(7)+T(8))/Nnum;
P12num = 1-(T(1)+T(4))/Nnum;
P23num = 1-(T(1)+T(2))/Nnum;
P13num = 1-(T(1)+T(3))/Nnum;
P123num = 1 - T(1)/Nnum;

%%
gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;
gamma3p = gamma3+s3*rho;
R1p = beta1/(gamma1p+mu);
R2p = beta2/(gamma2p+mu);
R3p = beta3/(gamma3p+mu);
N = b/mu;

%1 disease equ
P1 = 1-1/R1p;
P2 = 1-1/R2p;
P3 = 1-1/R3p;

% 2 diseases equ.
gamma12t = s1*s2*rho;gamma23t = s3*s2*rho;gamma13t = s1*s3*rho;

gamma1t2 = gamma1p-gamma12t;gamma2t1 = gamma2p-gamma12t;
S12= b/mu*(gamma1t2/R2p+gamma2t1/R1p+mu+gamma12t)/(beta1+beta2-mu-gamma12t);
P12 = 1 - S12/N;

gamma2t3 = gamma2p-gamma23t;gamma3t2 = gamma3p-gamma23t;
S23= b/mu*(gamma2t3/R3p+gamma3t2/R2p+mu+gamma23t)/(beta2+beta3-mu-gamma23t);
P23 = 1 - S23/N;

gamma1t3 = gamma1p-gamma13t;gamma3t1 = gamma3p-gamma13t;
S13= b/mu*(gamma1t3/R3p+gamma3t1/R1p+mu+gamma13t)/(beta1+beta3-mu-gamma13t);
P13 = 1 - S13/N;

% 3 diseases equ. %test
% systeme a resoudre %see draft of 15/02/22 pages 2-3
lambda1=beta1*(1-1/R1p);
lambda2=beta2*(1-1/R2p);
lambda3=beta3*(1-1/R3p);
S=T(1);I1=T(2);I2=T(3);I3=T(4);I12=T(5);I23=T(6);I13=T(7);I123=T(8);N=b/mu;
gamma2t3*(I2+I12)+gamma3t2*(I3+I13) - ((lambda2+lambda3+mu+gamma23t)*S23-b-N*gamma23t)
%equivalently
gamma2t3*(-S+I12)+gamma3t2*(-S+I13) - (-gamma2t3*S13-gamma3t2*S12+(lambda2+lambda3+mu+gamma23t)*S23-b-N*gamma23t)

%% Solutions theoriques
gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;
gamma3p = gamma3+s3*rho;
R1p = beta1/(gamma1p+mu);
R2p = beta2/(gamma2p+mu);
R3p = beta3/(gamma3p+mu);
N = b/mu;
lambda1=beta1*(1-1/R1p);
lambda2=beta2*(1-1/R2p);
lambda3=beta3*(1-1/R3p);
gamma12t = s1*s2*rho;
gamma13t = s1*s3*rho;
gamma23t = s2*s3*rho;
gamma1t2 = gamma1p-s1*s2*rho;
gamma1t3 = gamma1p-s1*s3*rho;
gamma2t1 = gamma2p-s1*s2*rho;
gamma2t3 = gamma2p-s3*s2*rho;
gamma3t1 = gamma3p-s1*s3*rho;
gamma3t2 = gamma3p-s3*s2*rho;

%I123
S13   = b/mu*(gamma1t3/R3p+gamma3t1/R1p+mu+gamma13t)/(beta1+beta3-mu-gamma13t); %T(1)+T(3)
I13_2 = b/mu*(1-1/R1p-1/R3p)+S13; %I13+I123 (prevalence des coinfectes 1 et 3 I13, voir SIS^2)
S23   = b/mu*(gamma2t3/R3p+gamma3t2/R2p+mu+gamma23t)/(beta2+beta3-mu-gamma23t); %T(1)+T(2)
I23_1 = b/mu*(1-1/R2p-1/R3p)+S23; %I23+I123, T(6)+T(8)
S12   = b/mu*(gamma1t2/R2p+gamma2t1/R1p+mu+gamma12t)/(beta1+beta2-mu-gamma12t); %T(1)+T(4)
I12_3 = b/mu*(1-1/R1p-1/R2p)+S12; %I12+I123 T(5)+T(8)
gamma123 = gamma1+gamma2+gamma3+(1-(1-s1)*(1-s2)*(1-s3))*rho;

I123th = (lambda1*I23_1+lambda2*I13_2+lambda3*I12_3)/(lambda1+lambda2+lambda3+gamma123+mu)

%S123
Sth = S12+S23+S13+b/mu*(1-1/R1p-1/R2p-1/R3p)-I123th
Sth-T(1)

