%% Script mymodel SIRxSIS (V3)
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


gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;

R1 = beta1/(gamma1p +mu); %il faudrait recalculer les R
R2 = beta2/(gamma2p +mu)  %il faudrait recalculer les R

s=1;
% Parametres du systeme d'ODE 
tspan = 0:1:10000;
Y0 = [S0; I10; I20; I30; I40; R10];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_V3(t,Y,b,beta1,beta2,s,s1,s2,rho,gamma1,gamma2,mu,'frequency'),tspan,Y0, options);

Ys(end,:)

beta2/(gamma2p - s*rho+mu)%pas ça
beta2/(gamma2p + mu) %ça a l'air d'etre ca
beta2/(gamma2p - s*rho/2+mu)%pas ça
beta2/(gamma2p + s*rho/2+mu)%pas ça

%Equilibria of the version 2
S_th = b/mu*(mu*(R1-1)+gamma2p+mu+s*rho*R2*mu/b)/(R1*(mu*(R1-1)+beta2+s*rho*R2*mu/b)) %ok si rho=0, 3.694
I1_th = (R1-1)/(beta2 + gamma1)*b*(gamma2/beta1 + (mu+gamma2/R1)/(mu*(R1-1) + beta2))
A = (R2 - 1)/(gamma2 + R1*mu)/(R2 - 1 + mu*R1*R2/beta2)
I2_th = A*b/(1- gamma2*A)

%Finding the equilibria of version 2
Ys(end,:)

gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;
R1 = beta1/(gamma1p + mu); %il faudrait recalculer les R
R2 = beta2/(gamma2p + mu); %il faudrait recalculer les R

S_th = b/mu*((mu)*(R1-1)+gamma2p+mu + s*rho*R1*mu/b*Ys(end,4))/(R1*((mu)*(R1-1)+beta2+ s*rho*R1*mu/b*Ys(end,4))) 
S1 = b/mu*(gamma1p+mu)/beta1
I1_th = (R1-1)/(beta2 + gamma1p)*b*(gamma2p/beta1 + (mu+gamma2p/R1)/(mu*(R1-1) + beta2)) %chercher encore
A = (R2 - 1)/(gamma2p + R1*(mu+s*rho))/(R2 - 1 + (mu)*R1*R2/beta2); I2_th = A*b/(1-gamma2p*A) %chercher encore

%I1 + I12
b/beta1*(R1-1)+s*rho/beta1*R1*Ys(end,4) 
Ys(end,2) + Ys(end,4)

%voir page 4 200520
(b/beta1*(R1-1)-Ys(end,2))*(gamma1p+mu)/(gamma1p+mu-s*rho)

%I12, voir 2205 page 2
num = b*(beta2-mu) - (beta2-gamma2p + gamma1p)*b/beta1*(R1-1)*(beta2-mu)-(beta2-gamma2p+gamma1p)*(b-mu*R1 + gamma2p*b/mu/R1);
denom = R1*((beta2-gamma2p+gamma1p)*s*rho/beta1*(beta2-mu)-beta2*mu); %pas bon

num/denom

%verif S en fct de I12 22/05 page 1
y=Ys(end,4) %I12
num = b + gamma2p*b/(mu*R1) +s*rho*y;
denom = mu*(R1-1) + s*rho*R1*mu/b*y + beta2*(1-1/R2)+gamma2p+mu;
num/denom %presque bon ,erreur d'arrondi ?
Ys(end,:)


Lambda1 = beta1*(Ys(end,2)+Ys(end,4))*mu/b
mu*(R1-1) + s*rho*R1*mu/b*Ys(end,4) %ok
Lambda2 = beta2*(Ys(end,3)+Ys(end,4)+Ys(end,5))*mu/b
beta2*(1-1/R2) %ok

%verif de I12 en fct de S
y = Ys(end,1)
num = b -(mu*R1+beta2*(1-1/R2))*y+ gamma2p*(b/mu/R1-y)
denom = mu/b*s*rho*R1*y
num/denom %pas bon

%verif de R1+Ir2
gamma1p/mu*b/beta1*(R1-1) + s*rho*gamma1p/mu/beta1*(R1 - beta1/gamma1p)*Ys(end,4)
Ys(end, 5) + Ys(end,6) %ok

%I12 en fct de S
y = Ys(end,1);
num = b-(mu*(R1-1) + beta2)*y+gamma2p/R1*b/mu;
denom = s*rho*(R1*mu/b*y-1)
res = num/denom %ok


%% Solving I12 (see 220520 page 2 verso)
clear all

syms beta1 beta2 gamma1 gamma2 mu b s1 s2 rho R1 R2 gamma1p gamma2p s
gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;
R1 = beta1/(mu+gamma1p);
R2 = beta2/(mu+gamma2p);

S = b/mu*(mu*(R1-1)+gamma2p+mu + s*rho*R1*mu/b*y)/(R1*(mu*(R1-1)+beta2+ s*rho*R1*mu/b*y)) 
Lambda2 = beta2*(1-1/R2); 
Lambda1 = mu*(R1-1) + s*rho*R1*mu/b*y;
I1 = b/beta1*(R1-1)+s*rho/beta1*R1*y - y;%ok

dSpdI1 = b - (Lambda2+mu+gamma2p)*S - (Lambda2+mu+gamma1p)*I1 + gamma2p*b/mu/R1 + gamma2p*y;

resI12= solve(dSpdI1==0,y)
Ys(end,:)
syms S
Lambda2 = beta2*(1-1/R2); 
I12 = (b-(mu*(R1-1)+beta2)*S + gamma2p/R1*b/mu)/(s*rho*(R1*mu/b*S-1))
I12 = Ys(end,4)
S2= (b-(Lambda2+mu+gamma1p)*b/beta1*(R1-1) + gamma2p*b/mu/R1)/beta2+...
    (gamma2p-(Lambda2+mu+gamma1p)*(s*rho/(beta1)*R1-1))/beta2*I12

eqn = S-S2==0
res = solve(eqn,S)


