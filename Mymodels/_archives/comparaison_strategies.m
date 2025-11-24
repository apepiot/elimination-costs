
%% Comparaison des deux strategies
clear all; 
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');

% mu  = 0;
% s2  = 1;
% s1  = s2;
% gamma1 = beta1/3;
% gamma2 = beta2/3;

% beta1 = 6.44
% beta2 = 3.79
% gamma1=5.28
% gamma2=1.86
% s1=0.81
% mu=0.41
% s2=0.98


R1 = beta1/(gamma1+mu);
R2 = beta2/(gamma2+mu);

alpha1 = beta1/s1*(1-1/R1);
alpha2 = beta2/s2*(1-1/R2);
maxalpha = max(alpha1,alpha2);

rho = 0:maxalpha/1000:(maxalpha);

gamma1p = gamma1+s1.*rho;
gamma2p = gamma2+s2.*rho;
R1p = beta1./(gamma1p+mu);
R2p = beta2./(gamma2p+mu);

%strategy 1 (target testing)
P12_1 = max(1 - (gamma1p./R2p + gamma2p./R1p + mu)./(beta1+beta2-mu),0); %COINFECTION EQU

%strategy 2 (kits)
P12_2 = max(1 - (gamma1p./R2p + gamma2p./R1p + mu + s1*s2*rho.*(1-1./R1p-1./R2p))./(beta1+beta2-mu -s1*s2.*rho),0); %COINFECTION EQU

U12_1 = 2*rho.*P12_1;
U12_2 = rho.*P12_2;


figure()
subplot(2,2,1)
plot(rho, P12_1, rho, P12_2)
title('\Pi_{12}')
legend('\Pi_{12}^T target', '\Pi_{12}^K kit')

subplot(2,2,2)
plot(rho, P12_1-P12_2)
title('\Pi_{12}^T - \Pi_{12}^K')

subplot(2,2,3)
plot(rho, U12_1,rho, U12_2)
legend('U_{12}^T target', 'U_{12}^K kit')
title('U_{12}')

subplot(2,2,4)
plot(rho, U12_1-U12_2)
title('U_{12}^T - U_{12}^K')
sgt = sgtitle([{'Comparaison of \Pi_{12} and U_{12} of the two strategies'},...
    {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
        num2str(gamma2), ' s1=', num2str(s1),' s2=', num2str(s2), '\mu=', num2str(mu)]}, ...
        {[' R_0^1=' num2str(round(R1,2)), ' R_0^2=' num2str(round(R2,2)),...
        ' \rho_1''=',num2str(round(alpha1,2)),' \rho_2''=',num2str(round(alpha2,2))]}])

sgt.FontSize = 15;

%% symbolic
clear all 
syms beta1 beta2 gamma1 gamma2 s1 s2 mu rho

mu=0;
s2=1;
s1=s2;
%gamma1=beta1/4;
%gamma2=beta2/4;

R1 = beta1/(gamma1+mu);
R2 = beta2/(gamma2+mu);
gamma1p = gamma1+s1.*rho;
gamma2p = gamma2+s2.*rho;
R1p = beta1./(gamma1p+mu);
R2p = beta2./(gamma2p+mu);

P12_1 = 1 - (gamma1p./R2p + gamma2p./R1p + mu)./(beta1+beta2-mu); %COINFECTION EQU (target)
P12_2 = 1 - (gamma1p./R2p + gamma2p./R1p + mu + s1*s2*rho.*(1-1./R1p-1./R2p))./(beta1+beta2-mu-s1*s2*rho); %COINFECTION EQU (kit)

U12_1 = 2*rho.*P12_1;
U12_2 = rho.*P12_2;


%% Difference of U
f = U12_1 - U12_2
solve(f==0,rho)

 %% Difference of Pi
g = P12_1 - P12_2
sol = solve(g==0,rho)

% Case mu=0, s1=s2=1;
%              0
%  beta1 - gamma1 = alpha1
%  beta2 - gamma2 = alpha2
alpha1 = beta1/s1*(1-1/R1);
alpha2 = beta2/s2*(1-1/R2);
g -  rho.*(rho-alpha1).*(rho-alpha2)./(beta1.*beta2.*(beta1+beta2-rho)); %OK

%General case
%                          0
%  -(gamma1 - beta1 + mu)/s1
%  -(gamma2 - beta2 + mu)/s2

Pi120 = 1 - (gamma1/R2+gamma2/R1+mu)/(beta1+beta2-mu);
lambda = s1*s2/(beta1*beta2)*(1 + ((1-Pi120)*R1*R2 -1)/((R1-1)*(R2-1)));
g - s1*s2*lambda*rho.*(rho-alpha1).*(rho-alpha2)./(beta1+beta2-mu-rho*s1*s2)



%% 19/08
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true');


R1    = beta1/(gamma1+mu);
R2    = beta2/(gamma2+mu);
Pi120 = 1 - (gamma1/R2+gamma2/R1+mu)/(beta1+beta2-mu);

f = Pi120*R1*R2 + 2 - R1 - R2
