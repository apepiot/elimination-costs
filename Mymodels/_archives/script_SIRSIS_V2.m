%% Script mymodel SIRxSIS (V2 - strategy 1 - U=rho1+rho2*P12 ? pas clair)
clear all
close all
%% Parametres
% Populations initiales
S0      = 50;
I10     = 5;
I20     = 3;
I30     = 2;    %I_12 (coinfection)
I40     = 0;
R10     = 0.;

% Definition des parametres 
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters('true', 'true')
rho1 = rand()*rho;
rho2 = rand()*rho;

% Parametres du systeme d'ODE 
tspan = 0:1:10000;
Y0 = [S0; I10; I20; I30; I40; R10];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_V2(t,Y,b,beta1,beta2,s1,s2,gamma1,gamma2,mu,rho1,rho2,'frequency'),tspan,Y0, options);
Ys(end,:)
%%
% Plot
figure(1)
%couleurs = {[0.9 0.9 0.9]; [0.8 0.8 0.8]; [0.6 0.6 0.6];[0.4 0.4 0.4]; [0.2 0.2 0.2];[0.4 0.4 0.4]}
plot(ts,Ys)%./(sum(Ys')'));
title('Evolution');
xlabel('Time t');
ylabel('Y');


% Y0 = [0.95;0;0;0;0;0];
% [R,F,V] = R_V2(Y0,b,beta1,beta2,e1,e2,gamma1,gamma2,mu);
% R1 = R(1); R2 = R(2); 

%% Frequency-dependent transmission case

%Equilibrium of 1
gamma1p = gamma1+s1*rho1;
gamma2p = gamma2+s2*rho2;
SE1 = (gamma1p + mu)/beta1*b/mu;
IE1 = b/beta1*(beta1/(gamma1p+mu) - 1);
RE1 = gamma1p*b/beta1/mu*(beta1/(gamma1p+mu) - 1);
E1 = [SE1,IE1, 0, 0, 0, RE1];

%Equilibrium of 2
SE2 = b/mu*(gamma2p+mu)/beta2;
IE2 = b/mu - SE2;
E2 = [SE2, 0, IE2, 0,0,0];

%Coequilibrium
gamma1p = gamma1+s1*rho1;
gamma2p = gamma2+s2*rho2;
R1p = beta1/(gamma1p+mu);
R2p = beta2/(gamma2p+mu);
SE12 = b/mu*(mu*(R1p-1)+gamma2p+mu)/(R1p*(mu*(R1p-1)+beta2))
I2E12 = b/mu/R1p - SE12;
I1E12 = b*(R1p-1)/(beta2+gamma1p)*(gamma2p/beta1+mu/b*SE12)
I12E12 = b/beta1*(R1p-1)-I1E12; 
R1E12 = b/mu/R2p -SE12- I1E12;
IR2E12 = b/mu*(1-1/R2p)-I2E12-I12E12


P_tot = mu/b*( b*(R01 -1)/(beta2+gamma1)*(gamma2/beta1 + 1/R01*(mu*(R01-1)+gamma2+mu)/(mu*(R01-1)+beta2)) + b/mu*(1-1/R02))


%% Density dependent transmission case
% Equilibria (to verify) (density dependent case)
E0 = [b/mu; 0;0;0;0;0]';
Ntot = b/mu%sum(Ys(end,:));
E2 = [0,1-(mu+gamma2)/beta2, (mu + gamma2)/beta2,0,0,0];
E3 = [S0/R1, b/mu-S0/R1, 0,0,0,0]
IE1 = (beta1*b- gamma1*mu - mu^2)/(beta1*(gamma1+mu));
SE1 = b/mu - IE1*(gamma1+mu)/mu;
SE1 = (gamma1 + mu)/beta1;
RE1 = gamma1*IE1/mu;
RE1 = -gamma1*(gamma1*mu + mu^2-beta1*b)/(mu*beta1*(gamma1+mu));
RE1 = Ntot - SE1 - IE1;
E1 = [SE1, IE1, 0,0,0,RE1];

IE2 = (beta2*b/mu - gamma2 - mu)/beta2;
SE2 = b/mu - IE2;
SE2 = (gamma2+mu)/beta2;

E2 = [SE2, 0, IE2, 0,0,0];

Ys(end,:)%/sum(Ys(end,:))


%% Utility function
%% Parametres
% Populations initiales
S0      = 50;
I10     = 5;
I20     = 3;
I30     = 2;    %I_12 (coinfection)
I40     = 0;
R10     = 0.;

% Definition des parametres 
b      = 0.3; %susceptibles birthrate
beta1   = 4;
beta2   = 2;
e1      = beta1;
e2      = beta2;
mu      = 0.01; %natural death rate
gamma1  = 0.4;
gamma2  = 0.20;

%% U2
vecDelta2 = 0:0.005:5;
vecDelta1 = vecDelta;
s1 = 0.95;
s2 = 0.9;
[U,vecPrev, vecR10,vecR20] = utilityFunctionSIRxSIS(beta1,beta2,gamma1,gamma2,vecDelta1, vecDelta2, s1, s2, b,mu,e1,e2);
U2 = vecDelta2.*vecPrev/b*mu;


figure(2)
plot(vecDelta2, vecPrev)
title([{'Prevalence of the SISxSIR model'},...
    {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
    num2str(gamma2), ' s1=', num2str(s1),' s2=', num2str(s2), '\mu=', num2str(mu), ' \pi=',  num2str(b), ' R_0^1=' num2str(round(beta1/(gamma1+mu),2)), ' R_0^2=' num2str(round(beta2/(gamma2+mu),2)) ]}])
xlabel("Voluntary-testing rate \delta_2","fontweight","bold")
ylabel("I(\delta_2)","fontweight","bold")


figure(3)
plot(vecDelta2, U2)
title([{'Utility of the SISxSIR model'},...
    {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
    num2str(gamma2), ' s1=', num2str(s1),' s2=', num2str(s2), '\mu=', num2str(mu), ' \pi=',  num2str(b), ' R_0^1=' num2str(round(beta1/(gamma1+mu),2)), ' R_0^2=' num2str(round(beta2/(gamma2+mu),2)) ]}])
xlabel("Voluntary-testing rate \delta_2","fontweight","bold")
ylabel("U(\delta_2)","fontweight","bold")
hold on;


[Umax, imax] = max(U2);
delta2Num = vecDelta2(imax)

%Verification par la théorie 
R01 = vecR10;
R02 = vecR20;
%R02 = beta2/(gamma2+mu);
gamma1p = gamma1 + s1*vecDelta1;
A = 2*(mu*(R01-1)./(beta2+gamma1p).*(1/beta1 + 1./(R01.*(mu*(R01-1)+beta2))) - 1/beta2);
B = vecPrev(1)/b*mu;
B = (b*(R01 -1)./(beta2+gamma1).*(gamma2/beta1 + 1./R01.*(mu*(R01-1)+gamma2+mu)./(mu*(R01-1)+beta2)) + b/mu.*(1-1./R02))/b*mu;
-B/A/s2
res = -(mu*((b*((gamma2 + mu)./beta2 - 1))./mu - (b*(gamma2./beta1 + ((gamma1 + mu + s1*vecDelta1).*(gamma2 + mu + mu.*(beta1./(gamma1 + mu + s1*vecDelta1) - 1)))./(beta1*(beta2 + mu.*(beta1./(gamma1 + mu + s1*vecDelta1) - 1)))).*(beta1./(gamma1 + mu + s1*vecDelta1) - 1))./(beta2 + gamma1)))./(b*s2*(2/beta2 - (2*mu*(1/beta1 + (gamma1 + mu + s1*vecDelta1)/(beta1*(beta2 + mu.*(beta1./(gamma1 + mu + s1*vecDelta1) - 1)))).*(beta1./(gamma1 + mu + s1*vecDelta1) - 1))./(beta2 + gamma1 + s1*vecDelta1)))
 

%3.7
%% Delta1 = Delta2
s1= 0.8;
s2= 0.7;
delta = 0.8;

% Parametres du systeme d'ODE 
tspan = 0:1:10000;
Y0 = [S0; I10; I20; I30; I40; R10];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_V2(t,Y,b,beta1, beta2,e1,e2,gamma1+s1*delta,gamma2+s2*delta,mu, 'frequency'),tspan,Y0, options);

%verif de la prevalence
sum(Ys(end, 2:5))*mu/b
% Par la theorie
gamma1p = gamma1 + s1*delta;
gamma2p = gamma2 + s2*delta;
R1 = beta1./(gamma1p + mu);
R2 = beta2./(gamma2p + mu);
% Prevalence
PI1 = (mu*(R1-1)./(beta2+gamma1p)).*(gamma2p/beta1 + (mu*(R1-1)+gamma2p+mu)./(R1.*(mu*(R1-1)+beta2)));
PI2 = 1-1./R2;
P12 = PI1 + PI2;


% Prvalence et utility function
vecDelta = 0.996:0.00000005:0.998;
vecDelta2 = vecDelta;
vecDelta1 = vecDelta;


[U,vecPrev, vecR10,vecR20] = utilityFunctionSIRxSIS(beta1,beta2,gamma1,gamma2,vecDelta1, vecDelta2, s1, s2, b,mu,e1,e2);
%U3 = vecDelta.*vecPrev;

figure(3)
plot(vecDelta, U)
title([{'Utility of the SISxSIR model'},...
    {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
    num2str(gamma2), ' s1=', num2str(s1),' s2=', num2str(s2), '\mu=', num2str(mu), ' \pi=',  num2str(b), ' R_0^1=' num2str(round(beta1/(gamma1+mu),2)), ' R_0^2=' num2str(round(beta2/(gamma2+mu),2)) ]}])
xlabel("Voluntary-testing rate \delta","fontweight","bold")
ylabel("U(\delta)","fontweight","bold")
hold on;


[Umax, imax] = max(U);
deltaNum = vecDelta(imax)

derivTh = (mu.*((gamma2 + vecDelta.*s2)./beta1 + ((gamma1 + mu + vecDelta.*s1).*(gamma2 + mu + vecDelta.*s2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))./(beta1.*(beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))).*(beta1./(gamma1 + mu + vecDelta.*s1) - 1))./(beta2 + gamma1 + vecDelta.*s1) - (2.*vecDelta.*s2)./beta2 - (gamma2 + mu + vecDelta.*s2)./beta2 + (vecDelta.*mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1).*(s2./beta1 + (s1.*(gamma2 + mu + vecDelta.*s2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))./(beta1.*(beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1))) + ((s2 - (beta1.*mu.*s1)./(gamma1 + mu + vecDelta.*s1).^2).*(gamma1 + mu + vecDelta.*s1))./(beta1.*(beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1))) + (mu.*s1.*(gamma2 + mu + vecDelta.*s2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))./((beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)).^2.*(gamma1 + mu + vecDelta.*s1))))./(beta2 + gamma1 + vecDelta.*s1) - (vecDelta.*mu.*s1.*((gamma2 + vecDelta.*s2)./beta1 + ((gamma1 + mu + vecDelta.*s1).*(gamma2 + mu + vecDelta.*s2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))./(beta1.*(beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))).*(beta1./(gamma1 + mu + vecDelta.*s1) - 1))./(beta2 + gamma1 + vecDelta.*s1).^2 - (beta1.*vecDelta.*mu.*s1.*((gamma2 + vecDelta.*s2)./beta1 + ((gamma1 + mu + vecDelta.*s1).*(gamma2 + mu + vecDelta.*s2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))./(beta1.*(beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))))./((beta2 + gamma1 + vecDelta.*s1).*(gamma1 + mu + vecDelta.*s1).^2) + 1
%calcule avec matlab, derive de tout dU1+dU2 ??
%le resultat ci dessus est bof, un peu approximatif
[derivThmin, imax2] = min(abs(derivTh))
deltaNum2 = vecDelta(imax2)

% Par la theorie
gamma1p = gamma1 + s1*vecDelta1;
gamma2p = gamma2 + s2*vecDelta2;
R1 = beta1./(gamma1p + mu);
R2 = beta2./(gamma2p + mu);
% Prevalence
PI1 = mu*(R1-1)./(beta2+gamma1p).*(gamma2p/beta1 + (mu*(R1-1)+gamma2p+mu)./(R1.*(mu*(R1-1)+beta2)));
PI2 = 1-1./R2;
P12 = max(PI1 + PI2,0);

plot(vecDelta, vecPrev)
hold on
plot( vecDelta, P12)


%Fonction d utilite
Uth = vecDelta.*P12;

%trouver le max theorique

syms beta1 beta2 gamma1 gamma2 delta mu b s1 s2

gamma1p = gamma1 + s1*delta;
gamma2p = gamma2 + s2*delta;

R1p = beta1./(gamma1p + mu);
R2p = beta2./(gamma2p + mu);

PI1 = (mu*(R1p-1)./(beta2+gamma1p)).*(gamma2p/beta1 + (mu*(R1p-1)+gamma2p+mu)./(R1p.*(mu*(R1p-1)+beta2)));
PI2 = 1 - 1./R2p;

UI1 = delta.*PI1;
UI2 = delta.*PI2;
Utot = UI1 + UI2;

R2 = beta2/(gamma2+mu);
dUI2 = 1-1./R2 - 2*s2/beta2*delta;
%dUI1 = exprimer en fct de sous variables


res1 =diff(UI1, delta)
%res1 = (mu*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu + delta*s1)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1))))*(beta1/(gamma1 + mu + delta*s1) - 1))/(beta2 + gamma1 + delta*s1) + (delta*mu*(beta1/(gamma1 + mu + delta*s1) - 1)*(s2/beta1 + (s1*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1))) + ((s2 - (beta1*mu*s1)/(gamma1 + mu + delta*s1)^2)*(gamma1 + mu + delta*s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1))) + (mu*s1*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))/((beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1))^2*(gamma1 + mu + delta*s1))))/(beta2 + gamma1 + delta*s1) - (delta*mu*s1*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu + delta*s1)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1))))*(beta1/(gamma1 + mu + delta*s1) - 1))/(beta2 + gamma1 + delta*s1)^2 - (beta1*delta*mu*s1*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu + delta*s1)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))))/((beta2 + gamma1 + delta*s1)*(gamma1 + mu + delta*s1)^2)

res = diff(Utot, delta);
res = dUI2 + res1
%res = (mu*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu + delta*s1)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1))))*(beta1/(gamma1 + mu + delta*s1) - 1))/(beta2 + gamma1 + delta*s1) - (2*delta*s2)/beta2 - (gamma2 + mu)/beta2 + (delta*mu*(beta1/(gamma1 + mu + delta*s1) - 1)*(s2/beta1 + (s1*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1))) + ((s2 - (beta1*mu*s1)/(gamma1 + mu + delta*s1)^2)*(gamma1 + mu + delta*s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1))) + (mu*s1*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))/((beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1))^2*(gamma1 + mu + delta*s1))))/(beta2 + gamma1 + delta*s1) - (delta*mu*s1*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu + delta*s1)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1))))*(beta1/(gamma1 + mu + delta*s1) - 1))/(beta2 + gamma1 + delta*s1)^2 - (beta1*delta*mu*s1*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu + delta*s1)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + delta*s1) - 1)))))/((beta2 + gamma1 + delta*s1)*(gamma1 + mu + delta*s1)^2) + 1

eqn = res==0;
sol = solve(eqn,delta)


%%
% Definition des parametres 
b      = 0.3; %susceptibles birthrate
beta1   = 4;
beta2   = 2;
e1      = beta1;
e2      = beta2;
mu      = 0.01; %natural death rate
gamma1  = 0.4;
gamma2  = 0.20;
s1 = 0.8;
s2 = 0.9;

delta= 0.9975;

syms beta1 beta2 gamma1 gamma2 delta mu b s1 s2

%% Essai de dU/ddelta = 0 (voir brouillon 060520 page 1) %semble OK
%delta = 1.281025;

vecRes=[];
vecDelta = -5:0.1:10;
for delta = vecDelta

gamma1p = gamma1 + s1*delta;
Gamma1 = gamma1p + mu;
gamma2p = gamma2 + s2*delta;
R1p = beta1/Gamma1;
R1 = beta1/(gamma1+mu);
%R2p = beta2/(gamma2p+mu);
R2 = beta2/(gamma2+mu);


A = (beta1-Gamma1)* Gamma1*(beta2+gamma1p) + delta*(-s1*beta1*(beta2+gamma1p) - s1*Gamma1*(beta1-Gamma1));
%     xA1 =  -(beta2*gamma1 + beta2*mu + gamma1*mu - (beta1*(beta2 + gamma1)*(gamma1 + mu)*(beta1 + beta2 - mu))^(1/2) + gamma1^2)/(beta1*s1 + beta2*s1 + gamma1*s1)
%     xA2 = -(beta2*gamma1 + beta2*mu + gamma1*mu + (beta1*(beta2 + gamma1)*(gamma1 + mu)*(beta1 + beta2 - mu))^(1/2) + gamma1^2)/(beta1*s1 + beta2*s1 + gamma1*s1)
%     -(xA2 - delta)*(xA1 -delta)*s1^2*(beta1+beta2+gamma1) 
% %[solx,parameters,conditions] = solve(eqn,delta,'ReturnConditions',true)
B1 = gamma2p/beta1*(mu*(beta1-Gamma1)+beta2*Gamma1)^2;
%     xB11 = -gamma2/s2;
%     xB12 =  (mu*(gamma1 - beta1 + mu) - beta2*(gamma1 + mu))/(beta2*s1 - mu*s1); %xB1
%     - s2/beta1*(beta2-mu)^2*s1^2*(xB11 - delta)*(xB12 -delta)^2;

B2 = Gamma1*(mu+gamma2p/R1p)*(mu*(beta1-Gamma1)+ Gamma1*beta2);
%     xB21 = (mu*(gamma1 - beta1 + mu) - beta2*(gamma1 + mu))/(beta2*s1 - mu*s1); %xB1
%     xB22 = -(gamma1 + mu)/s1; %xC04
%     xB23 =  -(gamma1*s2 + gamma2*s1 + mu*s2 + (gamma1^2*s2^2 - ... %xC31
%         2*gamma1*gamma2*s1*s2 + 2*gamma1*mu*s2^2 + gamma2^2*s1^2 - 2*gamma2*mu*s1*s2 + mu^2*s2^2 - 4*beta1*mu*s1*s2)^(1/2))/(2*s1*s2);%xC31
%     xB24 =  -(gamma1*s2 + gamma2*s1 + mu*s2 - (gamma1^2*s2^2 - ... %xC32
%         2*gamma1*gamma2*s1*s2 + 2*gamma1*mu*s2^2 + gamma2^2*s1^2 - 2*gamma2*mu*s1*s2 + mu^2*s2^2 - 4*beta1*mu*s1*s2)^(1/2))/(2*s1*s2); %xC31
%     s1^3*(beta2-mu)*s2/beta1*(xB21 -delta)*(xB22-delta)*(xB23 -delta)*(xB24-delta)
%     
%     
B = B1 + B2;
%     xB1 = -(beta2*gamma1 + beta1*mu + beta2*mu - gamma1*mu - mu^2)/(beta2*s1 - mu*s1);%xB12
%     %xB2, xB3 ?
%     Bp = (-(beta2-mu)*(xB11-delta)*(xB21-delta) + s1*(xB22-delta)*(xB23-delta)*(xB24-delta));
%     Bp = -(beta2-mu)*(xB11-delta)*(xB21-delta) + s1*(xB22-delta)*(delta^2 + (gamma1*s2+gamma2*s1+mu*s2)/(s1*s2)*delta + (beta1*mu+gamma2*(gamma1+mu))/(s1*s2) )
%     (xB23-delta)*(xB24-delta) - (delta^2 + (gamma1*s2+gamma2*s1+mu*s2)/(s1*s2)*delta + (beta1*mu+gamma2*(gamma1+mu))/(s1*s2) )
%     s1^2*s2/beta1*(beta2-mu)*(xB21-delta)*Bp %B
%     
% %     xC31f2 = (beta1*mu+gamma2*(gamma1+mu))/(s1*s2); %xC31*xC32
% %     xC31p2 = - (gamma1*s2+gamma2*s1+mu*s2)/(s1*s2);
% %     %Resolution dun systeme
% %     lambda = s1;
% %     eqn1 = a*b*c - ((beta2-mu)*gamma2/s1*xB1 + s1*xC04*xC31f2) ==0;
% %     %a = -(((gamma1 + mu)*(beta1*mu + gamma2*(gamma1 + mu)))/(s1*s2) + (gamma2*(beta2 - mu)*(beta2*gamma1 + beta1*mu + beta2*mu - gamma1*mu - mu^2))/(s1*(beta2*s1 - mu*s1)))/(b*c)
% %     eqn2 = -(b*c+a*b+a*c) - ((beta2-mu)*(-gamma2/s1+xB1) - s1*(xC04*xC31p2 + xC31f2)) ==0;
% %     eqn3 = c+b+a -(-(beta2-mu) + s1*(xC04 +  xC31f2)) ==0;   
    
C1 = s2/beta1*(mu*(beta1-Gamma1)+beta2*Gamma1)^2;
%    %xC11 = (mu*(gamma1 - beta1 + mu) - beta2*(gamma1 + mu))/(beta2*s1 - mu*s1);
%     xC11 = xB1;
%     xC12 = xC11;
%     s2/beta1*(beta2-mu)^2*s1^2*(xC11 - delta)^2 %C1

C2 = Gamma1*((gamma2*s1+2*s1*s2*delta)/beta1+s2/R1)*(mu*(beta1-Gamma1)+beta2*Gamma1);
% 	  xC21 = (mu*(gamma1 - beta1 + mu) - beta2*(gamma1 + mu))/(beta2*s1 - mu*s1);
% 	  xC22 = -(gamma1 + mu)/s1;
%     xC23 = -(beta1*((s2*(gamma1 + mu))/beta1 + (gamma2*s1)/beta1))/(2*s1*s2);% = -(gamma1*s2 + gamma2*s1 + mu*s2)/(2*s1*s2)
%     -s1^3*s2*2/beta1*(beta2-mu)*(xC21 -delta)*(xC22-delta)*(xC23-delta)

C3 = (mu+gamma2p/R1p)*mu*s1*beta1;
%     xC31 = -(gamma1*s2 + gamma2*s1 + mu*s2 + (gamma1^2*s2^2 - 2*gamma1*gamma2*s1*s2 + 2*gamma1*mu*s2^2 + gamma2^2*s1^2 - 2*gamma2*mu*s1*s2 + mu^2*s2^2 - 4*beta1*mu*s1*s2)^(1/2))/(2*s1*s2)
%     xC32 =  -(gamma1*s2 + gamma2*s1 + mu*s2 - (gamma1^2*s2^2 - 2*gamma1*gamma2*s1*s2 + 2*gamma1*mu*s2^2 + gamma2^2*s1^2 - 2*gamma2*mu*s1*s2 + mu^2*s2^2 - 4*beta1*mu*s1*s2)^(1/2))/(2*s1*s2)
%     mu*(s1*s2)*s1*(xC31 -delta)*(xC32-delta) %C3
    
C0 = delta*(beta2+gamma1p)*(beta1-Gamma1)*Gamma1;
%     xC01 = 0; 
%     xC02 = -(gamma1 - beta1 + mu)/s1
%     xC03 = -(beta2 + gamma1)/s1
%     xC04 = -(gamma1 + mu)/s1
%     res = s1^3*delta*(xC02 -delta)*(xC03 -delta)*(xC04 -delta) %C0

C12 = C1 + C2;
    %xC121 = -(beta2*gamma1 + beta1*mu + beta2*mu - gamma1*mu - mu^2)/(beta2*s1 - mu*s1); %xB1
    %xC122 = -(beta2*s2 + 3*gamma1*s2 + gamma2*s1 + 2*mu*s2 + (beta2^2*s2^2 - 2*beta2*gamma1*s2^2 + 2*beta2*gamma2*s1*s2 - 4*beta2*mu*s2^2 + gamma1^2*s2^2 - 2*gamma1*gamma2*s1*s2 + 4*gamma1*mu*s2^2 + gamma2^2*s1^2 - 4*gamma2*mu*s1*s2 + 4*mu^2*s2^2 - 8*beta1*mu*s2^2)^(1/2))/(4*s1*s2);
    %xC123 = -(beta2*s2 + 3*gamma1*s2 + gamma2*s1 + 2*mu*s2 - (beta2^2*s2^2 - 2*beta2*gamma1*s2^2 + 2*beta2*gamma2*s1*s2 - 4*beta2*mu*s2^2 + gamma1^2*s2^2 - 2*gamma1*gamma2*s1*s2 + 4*gamma1*mu*s2^2 + gamma2^2*s1^2 - 4*gamma2*mu*s1*s2 + 4*mu^2*s2^2 - 8*beta1*mu*s2^2)^(1/2))/(4*s1*s2);
C23 = C2+C3; %ne trouve pas

C13 = C1 + C3;
%     xC131 = -(2*mu^3*s2 - beta1*mu^2*s2 - 4*beta2*mu^2*s2 + 2*beta2^2*mu*s2 + 2*gamma1*mu^2*s2 + beta1*mu*(- 4*beta2^2*s1*s2 - 4*beta2*gamma1*s2^2 + 4*beta2*gamma2*s1*s2 + 8*beta2*mu*s1*s2 - 4*beta2*mu*s2^2 + gamma1^2*s2^2 - 2*gamma1*gamma2*s1*s2 + 6*gamma1*mu*s2^2 + gamma2^2*s1^2 - 6*gamma2*mu*s1*s2 - 4*mu^2*s1*s2 + 5*mu^2*s2^2 - 4*beta1*mu*s1*s2 - 4*beta1*mu*s2^2)^(1/2) + 2*beta2^2*gamma1*s2 + 2*beta1*beta2*mu*s2 + beta1*gamma1*mu*s2 + beta1*gamma2*mu*s1 - 4*beta2*gamma1*mu*s2)/(2*(s1*s2*beta2^2 - 2*s1*s2*beta2*mu + s1*s2*mu^2 + beta1*s1*s2*mu))
%     xC132 = -(2*mu^3*s2 - beta1*mu^2*s2 - 4*beta2*mu^2*s2 + 2*beta2^2*mu*s2 + 2*gamma1*mu^2*s2 - beta1*mu*(- 4*beta2^2*s1*s2 - 4*beta2*gamma1*s2^2 + 4*beta2*gamma2*s1*s2 + 8*beta2*mu*s1*s2 - 4*beta2*mu*s2^2 + gamma1^2*s2^2 - 2*gamma1*gamma2*s1*s2 + 6*gamma1*mu*s2^2 + gamma2^2*s1^2 - 6*gamma2*mu*s1*s2 - 4*mu^2*s1*s2 + 5*mu^2*s2^2 - 4*beta1*mu*s1*s2 - 4*beta1*mu*s2^2)^(1/2) + 2*beta2^2*gamma1*s2 + 2*beta1*beta2*mu*s2 + beta1*gamma1*mu*s2 + beta1*gamma2*mu*s1 - 4*beta2*gamma1*mu*s2)/(2*(s1*s2*beta2^2 - 2*s1*s2*beta2*mu + s1*s2*mu^2 + beta1*s1*s2*mu))
 
    
C = C0*(C1+C2+C3);
D = Gamma1^2/mu*(1-1/R2-2*delta*s2/beta2)*(beta2+gamma1p)^2*(mu*(beta1-Gamma1)+beta2*Gamma1)^2;
%     xD1 = (mu*(gamma1 - beta1 + mu) - beta2*(gamma1 + mu))/(beta2*s1 - mu*s1);%idem xB1
%     xD2 = xD1;
%     xD3 = - (gamma2 - beta2 + mu)/(2*s2);
%     xD4 = -(beta2 + gamma1)/s1;
%     xD5 = xD4;
%     xD6 = -(gamma1 + mu)/s1;
%     xD7 = xD6;
%     
%     res = (s1^6/mu * 2*s2/beta2 *(beta2-mu)^2)*(xD1 -delta)^2*(xD3-delta)*(xD4 -delta)^2*(xD6-delta)^2 %D

res = A*B + C +D;

vecRes = [vecRes, res];
end

plot(vecDelta, vecRes)
dU*Gamma1^2/mu*(beta2+gamma1p)^2*(mu*(beta1-Gamma1)+beta2*Gamma1)^2 %je ne sais plus ce que j'ai voulu faire ici, certainement comparer les resultats


%% On repart de la derivee dP1/ddelta (voir 050520 page 3)

vecdU = [];
vecDelta = 0:0.0005:2;
for delta=vecDelta
    dU = dU_SIRSIS(beta1,beta2,gamma1,gamma2,delta,s1, s2, mu);
    vecdU = [vecdU;dU];
end

[mindU, imin] = min(abs(vecdU))
vecDelta(imin)

plot(vecDelta, vecdU)

%resultat matlab
vecDelta = 1.2810
derivTh = (mu.*((gamma2 + vecDelta.*s2)./beta1 + ((gamma1 + mu + vecDelta.*s1).*(gamma2 + mu + vecDelta.*s2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))./(beta1.*(beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))).*(beta1./(gamma1 + mu + vecDelta.*s1) - 1))./(beta2 + gamma1 + vecDelta.*s1) - (2.*vecDelta.*s2)./beta2 - (gamma2 + mu + vecDelta.*s2)./beta2 + (vecDelta.*mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1).*(s2./beta1 + (s1.*(gamma2 + mu + vecDelta.*s2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))./(beta1.*(beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1))) + ((s2 - (beta1.*mu.*s1)./(gamma1 + mu + vecDelta.*s1).^2).*(gamma1 + mu + vecDelta.*s1))./(beta1.*(beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1))) + (mu.*s1.*(gamma2 + mu + vecDelta.*s2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))./((beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)).^2.*(gamma1 + mu + vecDelta.*s1))))./(beta2 + gamma1 + vecDelta.*s1) - (vecDelta.*mu.*s1.*((gamma2 + vecDelta.*s2)./beta1 + ((gamma1 + mu + vecDelta.*s1).*(gamma2 + mu + vecDelta.*s2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))./(beta1.*(beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))).*(beta1./(gamma1 + mu + vecDelta.*s1) - 1))./(beta2 + gamma1 + vecDelta.*s1).^2 - (beta1.*vecDelta.*mu.*s1.*((gamma2 + vecDelta.*s2)./beta1 + ((gamma1 + mu + vecDelta.*s1).*(gamma2 + mu + vecDelta.*s2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))./(beta1.*(beta2 + mu.*(beta1./(gamma1 + mu + vecDelta.*s1) - 1)))))./((beta2 + gamma1 + vecDelta.*s1).*(gamma1 + mu + vecDelta.*s1).^2) + 1;
figure(2)
plot(vecDelta, derivTh)


%% Etude des signes ?
vecDelta =-3:0.05:3;
[U,vecPrev, vecR10,vecR20] = utilityFunctionSIRxSIS(beta1,beta2,gamma1,gamma2,vecDelta, vecDelta, s1, s2, b,mu,e1,e2);
plot(vecDelta, U)

%bof 
res = A*B+C+D; %dU/dd %pas tout a fait
res1 = diff(res, delta)  %dU²/d²d
res2 = diff(res1, delta) %dU³/d³d
res3 = diff(res2, delta) %dU⁴/d⁴d
res4 = diff(res3, delta) %dU⁵/d⁵d
eqn4 = res4==0
solve(eqn4,delta)

vecDelta = -4.001:0.002:2;
plot(vecDelta, eval( (subs(res,delta,vecDelta)) ))

res5 = diff(res4, delta) %dU⁶/d⁶d
eqn5 = res5==0;
solve(eqn5, delta)
res6 = diff(res5, delta) %dU⁷/d⁷d
eqn6 = res6 ==0;
solve(eqn6, delta)
% x6 = -(360*s1^3*(2*mu*s1^2*s2 + (2*s2*(beta2*s1 - mu*s1)^2)/beta1 + 2*s1*((s2*(gamma1 + mu))/beta1 + (gamma2*s1)/beta1)*(beta2*s1 - mu*s1) - (4*s1^2*s2*(mu*(gamma1 - beta1 + mu) - beta2*(gamma1 + mu)))/beta1 + (4*s1*s2*(gamma1 + mu)*(beta2*s1 - mu*s1))/beta1) + (720*s1^4*((gamma2 + mu)/beta2 - 1)*(beta2*s1 - mu*s1)^2)/mu + (1440*s1^4*s2*(beta2*s1 - mu*s1)*(gamma1 - beta1 + mu))/beta1 + (1440*s1^4*s2*(beta2 + gamma1)*(beta2*s1 - mu*s1))/beta1 + (1440*s1^4*s2*(gamma1 + mu)*(beta2*s1 - mu*s1))/beta1 + (360*s1^2*s2*(2*s1^2*(beta2 + gamma1) + 2*beta1*s1^2)*(beta2*s1 - mu*s1))/beta1 - (2880*s1^4*s2*(mu*(gamma1 - beta1 + mu) - beta2*(gamma1 + mu))*(beta2*s1 - mu*s1))/(beta2*mu) + (2880*s1^3*s2*(beta2 + gamma1)*(beta2*s1 - mu*s1)^2)/(beta2*mu) + (2880*s1^3*s2*(gamma1 + mu)*(beta2*s1 - mu*s1)^2)/(beta2*mu))/((10080*s1^5*s2*(beta2*s1 - mu*s1))/beta1 + (10080*s1^4*s2*(beta2*s1 - mu*s1)^2)/(beta2*mu))

res7 = diff(res6, delta) %dU⁸/d⁸d
%- (10080*s1^5*s2*(beta2*s1 - mu*s1))/beta1 - (10080*s1^4*s2*(beta2*s1 - mu*s1)^2)/(beta2*mu)
res8 = diff(res7, delta) %dU⁹/d⁹d = 0

%% Essayons de minorer U %120520
vecDelta = 0:0.05:3;
[U,vecPrev, vecR10,vecR20] = utilityFunctionSIRxSIS(beta1,beta2,gamma1,gamma2,vecDelta, vecDelta, s1, s2, b,mu,e1,e2);
plot(vecDelta, U)
U2 = vecDelta.*(1-(gamma2+s2*vecDelta +mu)./beta2);
hold on;
plot(vecDelta, U2)

%
syms beta1 beta2 gamma1 gamma2 delta mu b s1 s2

gamma1p = gamma1 + s1*delta;
gamma2p = gamma2 + s2*delta;
R1p = beta1/(gamma1p+mu);
Gamma1 = gamma1p + mu;

U11 = delta*mu*(R1p -1)/(beta2+gamma1p)*gamma2p/beta1
dU11 = diff(U11, delta)

solve(dU11==0, delta)

dF = (beta2+gamma1p)*Gamma1*(beta1-Gamma1)*gamma2p/beta1 - (beta2+gamma1p)*delta*s1*gamma2p + delta*Gamma1*(beta1-Gamma1)*s1*gamma2p/beta1 + (beta2+gamma1p)*delta*Gamma1^2*s2
solve(dF==0, delta)

% non abouti...


%% Au pifometre 13/05/20
vecDelta =-0:0.0005:3;
[U,vecPrev, vecR10,vecR20] = utilityFunctionSIRxSIS(beta1,beta2,gamma1,gamma2,vecDelta, vecDelta, s1, s2, b,mu);
[UmaxNum  ,imax] = max(U)
dNum = vecDelta(imax)

%delta_chap du SIS:
R2 = beta2/(gamma2+mu);
dTh2 = beta2/(2*s2)*(1-1/R2)
UmaxTh2 = beta2/(4*s2)*(1-1/R2)^2

R1 = beta1/(gamma1+mu);
dTh1 = beta1/s1*(1/sqrt(R1)-1/R1)
UmaxTh1 = mu/s1*(1-1/sqrt(R1))^2


res = A*B+C*D;
%conditions (x7)
%s1 ~= 0 & s2 ~= 0 & beta2^2*mu + beta1*beta2^2 + beta1*mu^2 ~= beta2*mu^2 + 2*beta1*beta2*mu


vecExps1 = -3:0.5:7;
vecExps2 = -3:0.5:7;
vecExp3  = -3:0.5:7;
vecExp4  = -3:0.5:7;
vecExp5  = -3:0.5:7;
s3 = beta2^2*mu + beta1*beta2^2 + beta1*mu^2 -(beta2*mu^2 + 2*beta1*beta2*mu);
s4 = R1-1
s5 = R2-1

for exp1=vecExps1 
    for exp2=vecExps2
        for exp3 = vecExp3
            for exp4 = vecExp4
                for exp5 = vecExp5
                    res = s1^exp1*s2^exp2*s3^exp3*s4^exp4*s5^exp5;
                     if abs(res - 0.9975)<0.0001
                         res
                         exp1
                         exp2
                         exp3
                         exp4
                         exp5
                     end
                 end
    
             end
        end
    end
end


%% On trace la fonction d'utilite 14/05/20
clear all
close all

b      = 0.3; %susceptibles birthrate
beta1   = 0.2;
beta2   = 3;
mu      = 0.09; %natural death rate
gamma1  = 0.1;
gamma2  = 3;
s1 = 0.9;
s2 = 0.9;

v12inf2 = true;
t=0

%while v12inf2 && t<10000 
    %t= t+1 
    %[beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters(true, true)

    R1 = beta1/(mu+gamma1)
    R2 = beta2/(gamma2+mu)


    vecDelta = 0:0.0005:0.2;
    [U,vecPrev, vecR10,vecR20,P1,P2,P12] = utilityFunctionSIRxSIS(beta1,beta2,gamma1,gamma2,vecDelta,vecDelta,s1,s2,b,mu);
    hold off
    
    figure(1)
    plot(vecDelta, U, 'LineWidth',2)
    title([{'Utility of the SISxSIR model'},...
        {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
        num2str(gamma2), ' s1=', num2str(s1),' s2=', num2str(s2), '\mu=', num2str(mu), ' \pi=',  num2str(b), ...
        ' R_0^1=' num2str(round(beta1/(gamma1+mu),2)), ' R_0^2=' num2str(round(beta2/(gamma2+mu),2)),...
        ' \alpha_1=', num2str(round(beta1/s1*(1-1/R1),2)),...
        ' \alpha_2=', num2str(round(beta2/s2*(1-1/R2),2))]}])
    xlabel("Voluntary-testing rate \delta","fontweight","bold")
    ylabel("U(\delta)","fontweight","bold")
    hold on;
    
    plot(vecDelta, vecDelta.*max(P1,0),  '-.','LineWidth',3)
    plot(vecDelta, vecDelta.*max(P2,0),  '--','LineWidth',3)
    plot(vecDelta, vecDelta.*max(P12,0),  '.','MarkerSize',10)
    legend('U', 'U_1', 'U_2', 'U_{12}')

    alpha1 = beta1/s1*(1-1/R1)
    alpha2 = beta2/s2*(1-1/R2)

    deltaMaxth1 = beta1/s1*(1/sqrt(R1)-1/R1);
    deltaMaxth2 = beta2/(2*s2)*(1-1/R2);

    max1 = max(vecDelta.*P1);
    [max2, i2max] = max(vecDelta.*P2);
    [max12, i12max] = max(vecDelta.*P12);

    delta2Num = vecDelta(i2max)
    delta12Num = vecDelta(i12max)

    Umax1 = mu/(s1)*(1-1/sqrt(R1))^2;
    Umax2 = beta2/(4*s2)*(1-1/R2)^2;
    
%     if(delta2Num > delta12Num)
%         v12inf2 = false;
%     end
%end

% 
% plot(vecDelta, (vecR10>1).*U,  '-.','LineWidth',4)
% plot(vecDelta, (vecR20>1).*U,  '--','LineWidth',4)
% legend('U', 'U_{R_1>1}', 'U_{R_2>1}')
% 


% figure(2)
% 
% beta2 = 0:0.02:5;
% R2 = beta2./(gamma2+mu);
% Umax2 = beta2./(4*s2).*(1-1./R2).^2
% alpha2 = beta2./s2.*(1-1./R2)
% 
% plot(beta2, Umax2)



%% U1 < U2 180520

syms beta1 beta2 gamma1 gamma2 delta mu b s1 s2 R1 R2

R1p = beta1/(gamma1+mu+s1*delta);
R2p = beta2/(gamma2+mu+s2*delta);

II1 = mu/beta1*(R1p-1);
II2 = 1 - 1/R2p;

eqn = II1 - II2 ==0

solve(eqn , delta)


eqn = mu/beta1*(1-1/R1 - s1/beta1*delta) - (1/R1+s1/beta1*delta)*(1-1/R2 - s2/beta2*delta) ==0
 %-((R1^2*R2^2*beta1^2*beta2^2*s1^2 - 4*R1^2*R2^2*beta1^2*beta2*mu*s1*s2 + 2*R1^2*R2^2*beta1*beta2^2*mu*s1^2 + R1^2*R2^2*beta2^2*mu^2*s1^2 - 2*R1^2*R2*beta1^2*beta2^2*s1^2 - 2*R1^2*R2*beta1*beta2^2*mu*s1^2 + R1^2*beta1^2*beta2^2*s1^2 + 2*R1*R2^2*beta1^3*beta2*s1*s2 + 2*R1*R2^2*beta1^2*beta2*mu*s1*s2 - 2*R1*R2*beta1^3*beta2*s1*s2 + R2^2*beta1^4*s2^2)^(1/2) + R2*beta1^2*s2 + R1*beta1*beta2*s1 - R1*R2*beta1*beta2*s1 - R1*R2*beta2*mu*s1)/(2*R1*R2*beta1*s1*s2)
 % ((R1^2*R2^2*beta1^2*beta2^2*s1^2 - 4*R1^2*R2^2*beta1^2*beta2*mu*s1*s2 + 2*R1^2*R2^2*beta1*beta2^2*mu*s1^2 + R1^2*R2^2*beta2^2*mu^2*s1^2 - 2*R1^2*R2*beta1^2*beta2^2*s1^2 - 2*R1^2*R2*beta1*beta2^2*mu*s1^2 + R1^2*beta1^2*beta2^2*s1^2 + 2*R1*R2^2*beta1^3*beta2*s1*s2 + 2*R1*R2^2*beta1^2*beta2*mu*s1*s2 - 2*R1*R2*beta1^3*beta2*s1*s2 + R2^2*beta1^4*s2^2)^(1/2) - R2*beta1^2*s2 - R1*beta1*beta2*s1 + R1*R2*beta1*beta2*s1 + R1*R2*beta2*mu*s1)/(2*R1*R2*beta1*s1*s2)

 
 %% Conjecturer le Umax 180520
 clear all
 close all
 t=1;
 vecDelta=0:0.05:10;
 tableRes =[];
while t<200
    [beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters(true, true);
    [R1, R2, alpha1, alpha2, deltaMaxth1,deltaMaxth2,maxU,deltaUmax,max1,delta1Num,max2,delta2Num,max12,delta12Num] = findMaxU(beta1,beta2,gamma1,gamma2,vecDelta,s1,s2,b,mu);
    
    tableRes = [tableRes;  [R1, R2, alpha1, alpha2, deltaMaxth1,deltaMaxth2,maxU,deltaUmax,max1,delta1Num,max2,delta2Num,max12,delta12Num]];
    t = t+1;

end

T = array2table(tableRes,...
    'VariableNames',{'R1', 'R2', 'alpha1', 'alpha2', 'deltaMaxth1','deltaMaxth2','maxU','deltaUmax','max1','delta1Num','max2','delta2Num','max12','delta12Num'});
 T = T(T.R1>=1 & T.R2>=1 ,:);

T.A1 = ((T.alpha1 - T.alpha2)>0); %true if alpha1>alpha2, bcp vrai
T.A2 = (T.delta12Num <min(T.alpha1,T.alpha2)); %bcp vrai, semble etre faux quand R1 proche de 1, 


%% Approximation de U par developpement de Taylor
clear all

syms beta1 beta2 gamma1 gamma2 delta mu b s1 s2 R1 R2
gamma1p = gamma1 + s1*delta;
Gamma1 = gamma1p + mu;
gamma2p = gamma2 + s2*delta;
R1p = beta1/Gamma1;
R1 = beta1/(gamma1+mu);
%R2p = beta2/(gamma2p+mu);
R2 = beta2/(gamma2+mu);


PI1 = (mu*(R1-1)./(beta2+gamma1p)).*(gamma2p/beta1 + (mu*(R1-1)+gamma2p+mu)./(R1.*(mu*(R1-1)+beta2)));
PI2 = 1-1./R2;
P12 = PI1 + PI2;

U = delta*P12;
%dU =diff(U,delta)
dU = delta*((mu*(s2/beta1 + (s2*(gamma1 + mu))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1 + delta*s1) - (mu*s1*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1 + delta*s1)^2) - (gamma2 + mu)/beta2 + (mu*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1 + delta*s1) + 1

d2U = diff(dU,delta);
d2U = (2*mu*(s2/beta1 + (s2*(gamma1 + mu))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1 + delta*s1) - delta*((2*mu*s1*(s2/beta1 + (s2*(gamma1 + mu))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1 + delta*s1)^2 - (2*mu*s1^2*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1 + delta*s1)^3) - (2*mu*s1*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1 + delta*s1)^2

d3U = diff(d2U,delta)
d3U = delta*((6*mu*s1^2*(s2/beta1 + (s2*(gamma1 + mu))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1 + delta*s1)^3 - (6*mu*s1^3*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1 + delta*s1)^4) - (6*mu*s1*(s2/beta1 + (s2*(gamma1 + mu))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1 + delta*s1)^2 + (6*mu*s1^2*((gamma2 + delta*s2)/beta1 + ((gamma1 + mu)*(gamma2 + mu + delta*s2 + mu*(beta1/(gamma1 + mu) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1 + delta*s1)^3

%Developpement de Taylor en 0
U0  = 0;
dU0 = (mu*(gamma2/beta1 + ((gamma1 + mu)*(gamma2 + mu + mu*(beta1/(gamma1 + mu) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1) - (gamma2 + mu)/beta2 + 1;
d2U0 = (2*mu*(s2/beta1 + (s2*(gamma1 + mu))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1) - (2*mu*s1*(gamma2/beta1 + ((gamma1 + mu)*(gamma2 + mu + mu*(beta1/(gamma1 + mu) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1)^2;
d3U0 = (6*mu*s1^2*(gamma2/beta1 + ((gamma1 + mu)*(gamma2 + mu + mu*(beta1/(gamma1 + mu) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1)^3 - (6*mu*s1*(s2/beta1 + (s2*(gamma1 + mu))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1))))*(beta1/(gamma1 + mu) - 1))/(beta2 + gamma1)^2

%syms delta

taylor = U0 + dU0.*delta + d2U0.*delta.^2/2 + d3U0/6.*delta.^3;
%%
clear all
close all

b      = 0.3; %susceptibles birthrate
beta1   = 1.35;
beta2   = 2.11;
mu      = 0.05; %natural death rate
gamma1  = 0.39;
gamma2  = 1.44;
s1 = 0.83;
s2 = 0.99;

R1 = beta1/(gamma1+mu);
R2 = beta2/(gamma2+mu);

delta = -10:0.01:100;

taylor = (delta.^2.*((2.*mu*(s2./beta1 + (s2.*(gamma1 + mu))./(beta1.*(beta2 + mu*(beta1./(gamma1 + mu) - 1)))).*(beta1./(gamma1 + mu) - 1))./(beta2 + gamma1) - (2*mu*s1*(gamma2./beta1 + ((gamma1 + mu).*(gamma2 + mu + mu*(beta1./(gamma1 + mu) - 1)))./(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1)))).*(beta1/(gamma1 + mu) - 1))./(beta2 + gamma1)^2))/2 + delta.*((mu*(gamma2/beta1 + ((gamma1 + mu)*(gamma2 + mu + mu*(beta1/(gamma1 + mu) - 1)))./(beta1*(beta2 + mu*(beta1/(gamma1 + mu) - 1)))).*(beta1/(gamma1 + mu) - 1))./(beta2 + gamma1) - (gamma2 + mu)./beta2 + 1)

figure(1)
plot(delta, taylor)
hold on
[U,vecPrev, vecR10,vecR20,P1,P2,P12] = utilityFunctionSIRxSIS(beta1,beta2,gamma1,gamma2,delta,delta,s1,s2,b,mu);
plot(delta, delta.*P12)

% Approximation par degré 2 en fct des racines
alpha1 = beta1/s1*(1-1/R1)
alpha2 = beta2/s2*(1-1/R2)

aupif = delta.*(-delta + min(alpha1, alpha2)) %ca va donner U2

plot(delta, aupif)

%% Solving max U1 + U2 26/05
clear all

syms beta1 beta2 gamma1 gamma2 delta mu b s1 s2 R1 R2
gamma1p = gamma1 + s1*delta;
gamma2p = gamma2 + s2*delta;
R1p = beta1/(gamma1p+mu);
R1 = beta1/(gamma1+mu);
R2p = beta2/(gamma2p+mu);
R2 = beta2/(gamma2+mu);

U1p2 = delta*mu/beta1*(R1p-1) + delta*(1-1/R2p)

dU1p2 = diff(U1p2,delta)

res = solve(dU1p2 ==0 , delta)

 d = mu/beta1/R1-mu/beta1/R1^2 +1/R1^2 -1/R2/R1^2
 c = -mu/beta1*2*s1/R1 + 2*s1/R1 - 2*s1/R2/R1 -2*s2/beta2/R1^2;
 b = -mu/beta1*s1^2/beta1^2 + s1^2/beta1^2 - s1^2/beta1^2/R2 - 4*s2*s1/beta2/R1;
 a = -2*s2/beta2*s1^2/beta1^2;
 
dU1p2bis = a*delta^3 + b*delta^2 + c*delta +d
res = solve(dU1p2bis ==0 , delta)

%% 230620 solving rho1 > rho2
clear all
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');
gamma1 = 0:0.01:30;

R1 = beta1./(gamma1+mu);
R2 = beta2/(gamma2+mu);
rho1 = beta1/s1.*(1./sqrt(R1)- 1./R1);
alpha1 = beta1./s1.*(1-1./R1);
alpha2 = beta2./s2.*(1-1/R2);
rho2 = alpha2/2;

figure(1)
plot(gamma1,rho1-rho2,gamma1,alpha1-alpha2)


alpha2 > beta1/2/s1
%rho1>rho2

lim1 = beta1/4*(1-sqrt(1-2*s1/beta1*alpha2))^2-mu
lim2 = beta1/4*(1+sqrt(1-2*s1/beta1*alpha2))^2-mu
