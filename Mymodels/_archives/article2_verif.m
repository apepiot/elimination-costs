% Verification des résultats theoriques presentes dans l'article 2
clear all; close all;
S0      = 50;
I10     = 5;
I20     = 3;
I30     = 2;    %I_12 (coinfection)
I40     = 0;
R10     = 0.;

% Definition des parametres 
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true)
rho1=rand(1,1)*rho;rho2=rand(1,1)*rho;

%% Targeting testing strategy (script_SIRSIS_V5)
tspan = 0:1:10000;
Y0 = [S0; I10; I20; I30; I40; R10];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_V2(t,Y,b,beta1,beta2,s1,s2,gamma1,gamma2,mu,rho1,rho2,'frequency'),tspan,Y0, options);
T1=Ys(end,:)

R1 = beta1/(gamma1+mu);
R2 = beta2/(gamma2+mu);

gamma1p = gamma1+s1*rho1;
gamma2p = gamma2+s2*rho2;
R1p = beta1/(gamma1p+mu);
R2p = beta2/(gamma2p+mu);
Lambda1 = mu*(R1p-1);
Lambda2 = beta2*(1-1/R2p);

if R1>1 && R2<=1
    %ES1 
    S_ES1  = b/(mu*R1p);
    I1_ES1 = b/beta1*(R1p-1);
    R1_ES1 = gamma1p/mu*(R1p-1);
    ES1 = [S_ES1,I1_ES1,0,0,0,R1_ES1];
    ES=ES1
end
if R1<=1 && R2>1
    %ES2
    S_ES2  = b/(mu*R2p);
    I2_ES2 = b/mu*(1-1/R2p);
    ES2 = [S_ES2,0,I2_ES2,0,0,0];
    ES=ES2
end
if R1>1 & R2>1
%ES12
    S_ES12  = b/(mu*R1p)*(Lambda1+gamma2p+mu)/(Lambda1+beta2);
    I1_ES12 = b/mu*Lambda1/(beta2+gamma1p)*(gamma2p/beta1+mu/b*S_ES12);
    I2_ES12 = b/(mu*R1p)-S_ES12;
    I12_ES12= b*Lambda1/(mu*beta1)-I1_ES12;
    R1_ES12 = b/(mu*R2p)-S_ES12-I1_ES12;
    IR2_ES12= b/mu*(Lambda1*gamma1p/(beta1*mu)-1/R2p)+S_ES12 + I1_ES12;
    ES12 = [S_ES12,I1_ES12,I2_ES12,I12_ES12,IR2_ES12,R1_ES12]
    ES = ES12;
end
sum(ES)-b/mu
(T1-ES)./T1

%% Combined testing strategy (script_SIRSIS_V7)
clear rho1 rho2
tspan = 0:1:10000;
Y0 = [S0; I10; I20; I30; I40; R10];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_V7(t,Y,b,beta1,beta2,s1,s2,rho,gamma1,gamma2,mu,'frequency'),tspan,Y0, options);
T2=Ys(end,:)

R1 = beta1/(gamma1+mu);
R2 = beta2/(gamma2+mu);

gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;
R1p = beta1/(gamma1p+mu);
R2p = beta2/(gamma2p+mu);
Lambda1 = mu*(R1p-1);
Lambda2 = beta2*(1-1/R2p);
gamma12t = s1*s2*rho;
gamma1t  = gamma1p-s1*s2*rho;
gamma2t  = gamma2p-s1*s2*rho;

if R1>1 && R2<=1
    %ES1 
    S_ES1  = b/(mu*R1p);
    I1_ES1 = b/beta1*(R1p-1);
    R1_ES1 = gamma1p/mu*(R1p-1);
    ES1 = [S_ES1,I1_ES1,0,0,0,R1_ES1];
    ES=ES1
end
if R1<=1 && R2>1
    %ES2
    S_ES2  = b/(mu*R2p);
    I2_ES2 = b/mu*(1-1/R2p);
    ES2 = [S_ES2,0,I2_ES2,0,0,0];
    ES=ES2
end
if R1>1 & R2>1
%ES12
    S_ES12  = b/(mu*R1p)*(Lambda1+gamma2p+mu)/(Lambda1+beta2);
    I1_ES12 = b/mu*Lambda1/(beta2+gamma1t)*(gamma2t/beta1+mu/b*S_ES12);
    I2_ES12 = b/(mu*R1p)-S_ES12;
    I12_ES12= b*Lambda1/(mu*beta1)-I1_ES12;
    R1_ES12 = b/(mu*R2p)-S_ES12-I1_ES12;
    IR2_ES12= b/mu*(Lambda1*gamma1p/(beta1*mu)-1/R2p)+S_ES12+I1_ES12;
    ES12 = [S_ES12,I1_ES12,I2_ES12,I12_ES12,IR2_ES12,R1_ES12]
    ES = ES12;
end
sum(ES)-b/mu
(T2-ES)./T2


%% Utility function of combined testing
clear all; close all;
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true)
s1=1;
s2=1;
R1 = beta1/(gamma1+mu);
R2 = beta2/(gamma2+mu);
alpha1=beta1/s1*(1-1/R1);
alpha2=beta2/s2*(1-1/R2);

rho=0:min(alpha1,alpha2)/1000:max(alpha1,alpha2);

gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;
R1p = beta1./(gamma1p+mu);
R2p = beta2./(gamma2p+mu);
Lambda1 = mu*(R1p-1);
Lambda2 = beta2*(1-1./R2p);
gamma12t = s1*s2*rho;
gamma2t  = gamma2p-s1*s2*rho;
P12 = Lambda1./(beta2+gamma1p-gamma12t).*(gamma2t/beta1+(Lambda1+gamma2p+mu)./(R1p.*(Lambda1+beta2)))+1-1./R2p;
U12 = max(rho.*P12,0);

P1 = mu/beta1*(R1p-1);
P2 = 1-1./R2p;
minalpha = min(alpha1,alpha2);
P = P12.*(rho<minalpha) + P1.*(rho<alpha1 & rho>=alpha2) + P2.*(rho>=alpha1 & rho<alpha2);
U = rho.*P;

figure(1)
plot(rho,max(rho.*P12,0))
hold on;
plot(rho,max(rho.*P1,0));
plot(rho,max(rho.*P2,0));
plot(rho,U,'Linewidth',1.5)
maxU = max(U);
plot(minalpha*ones(1,100),0:maxU*1.1/99:maxU*1.1,'--')

%argmaxU12
clear rho;
maxalpha=max(alpha1,alpha2);
fun = @(rho) U_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,0);
argmaxU12 = fminbnd(fun,0,maxalpha)
maxU = - U_SIRSIS7(argmaxU12,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,0);

%rho/U12=0 i.e. rho/P12=0
rhoU120 = fzero(fun,(alpha1+alpha2)/2)%(alpha1+alpha2)/2

c=0;
rhohat1 = beta1./(R1*s1).*(sqrt(R1*mu./(mu-beta1*c))-1);
rhohat2 = alpha2/2;

%plot
legend('$U_{12}$','$U_1$','$U_2$','$U$','$min(\rho_1\prime,\rho_2\prime)$','Interpreter','latex')
title([{'SIR $\times$ SIS : $U(\rho)=\rho\Pi$'},...
      {['$\beta_1=$',num2str(beta1),' $\beta_2=$',num2str(beta2),...
      ' $\gamma_1(0)=$',num2str(round(gamma1,2)),' $\gamma_2(0)=$',num2str(round(gamma2,2)),...
      ' $s_1$=', num2str(s1),' $s_2$=', num2str(s2),...
      ' $\mu=$', num2str(round(mu,2)),' $\mathtt{R}_1(0)$=', num2str(round(R1,2)),' $\mathtt{R}_2(0)$=', num2str(round(R2,2))]},...
      {[' $\rho_1\prime$=', num2str(round(alpha1,2)),' $\rho_2\prime$=', num2str(round(alpha2,2)),...
      ' $\hat\rho_1=$', num2str(round(rhohat1,2)),' $\hat\rho_2=$', num2str(round(rhohat2,3)),...
      ' $\hat\rho_{12}=$', num2str(round(argmaxU12,3)), ' $\rho\prime\prime=$',num2str(round(rhoU120,3))]}],'Interpreter','latex')

%comparison between argmaxU12 and rho" (=rhoU120)
rhoU120/2 - argmaxU12