%% Script mymodel SIRxSIS (V7) - strategy 2 february 2021
% U = rho Pi ou Pi=I1+I2+I12+IR2 - I12-> R1
clear all
close all

%% Parametres
% Definition des parametres
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);

R1p=beta1/(gamma1+mu+s1*rho)%
R2p=beta2/(gamma2+mu+s2*rho)%

% Populations initiales
S0      = b/mu*1.1;
I10     = b/mu*0.1;
I20     = b/mu*0.1;
I30     = b/mu*0.01;    %I_12 (coinfection)
I40     = 0;
R10     = 0.;

% Parametres du systeme d'ODE
tspan = 0:0.1:10000;
Y0 = [S0; I10; I20; I30; I40; R10];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_V7(t,Y,b,beta1,beta2,s1,s2,rho,gamma1,gamma2,mu,'frequency'),tspan,Y0, options);
T = Ys(end,:)

sum(T(2:5))/b*mu

%comparaison avec la theorie
gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;
gamma2t = gamma2p-s1*s2*rho;
R1p = beta1/(gamma1p+mu);
R2p = beta2/(gamma2p+mu);
Lambda1 = mu*(R1p-1);
Lambda2 = beta2*(1-1/R2p);
P12 = 1-1/R2p + Lambda1/R1p*((gamma2t/(gamma1p+mu)+(Lambda1+gamma2p+mu)/(Lambda1+beta2))/(Lambda2+gamma1p+mu+gamma2t))

%% equilibria/refaire les verifs, erreurs de l'ordre 10^-6
gamma1p = gamma1+s1*rho;gamma2p=gamma2+s2*rho;
R1p=beta1/(gamma1p+mu);R2p=beta2/(gamma2p+mu);
lambda1 = mu*(R1p-1);lambda2=beta2*(1-1/R2p);
gamma2t = gamma2p-s1*s2*rho;

%Coinfection equilibrium : ES12 (a changer, ces resultats proviennent de la
%strategie 1)
SES12   = b/mu*(mu*(R1p-1) + gamma2p+mu)/(mu*(R1p-1)+beta2)/R1p;%ok
I2ES12  = b/(mu*R1p)-SES12;
%I12ES12 = b/beta1*(R1p-1)-I1ES12; strategy 1
%I12ES12 = mu*(R1p-1)*((beta2*(1-1/R2p)+gamma1p+mu)*b/mu/beta1-SES12)/(beta2-s1*s2*rho+gamma1p);%seems ok
I12ES12 = b/mu*lambda1/R1p*((lambda2+gamma1p+mu)/(gamma1p+mu) - (lambda1+mu+gamma2p)/(lambda1+beta2))/(beta2+gamma1p-s1*s2*rho)
%I1ES12  = b/mu*lambda1/(beta2+gamma1p)*(gamma2p/beta1+mu/b*SES12); %strategy 1
%I1ES12 = b/beta1*(R1p-1)-I12ES12; %seems ok
I1ES12 = b/mu*lambda1/R1p*(gamma2t/(gamma1p+mu)+(lambda1+gamma2p+mu)/(lambda1+beta2))/(lambda2+gamma1p+mu+gamma2t);
IR2ES12 = b/mu*(1-1/R2p)-I2ES12-I12ES12;
R1ES12  = b/(mu*R2p)-SES12 -I1ES12;
N = b/mu;

[SES12, I1ES12,I2ES12,I12ES12,IR2ES12,R1ES12]-T
%% solution par matlab
clear all
syms beta1 beta2 gamma1 gamma2 s1 s2 rho mu b S I1 I2 I12 IR2 R1
s1=1;s2=1;
gamma1p = gamma1+s1*rho;gamma2p=gamma2+s2*rho;
R1p = beta1/(gamma1p+mu); R2p = beta2/(gamma2p+mu);
lambda1 = mu*(R1p-1); lambda2=beta2*(1-1/R2p);
gamma1t = gamma1p-s1*s2*rho;
gamma2t = gamma2p-s1*s2*rho;
gamma12t = s1*s2*rho;
gamma12 = gamma1p+gamma2p-s1*s2*rho;

eqn1 = b - (lambda1 + lambda2)*S + gamma2p*I2 - mu*S==0; %S
eqn2 = lambda1*S - lambda2*I1 + gamma2t*I12 - (gamma1p + mu)*I1==0; %I1
eqn3 = lambda2*S - lambda1*I2 - (gamma2p + mu)*I2==0; %I2
eqn4 = lambda2*I1 + lambda1*I2 - (gamma12 + mu)*I12==0;%I12
eqn5 = gamma1t*I12 + lambda2*R1  - (gamma2p + mu)*IR2==0; %IR2
eqn6 = gamma1p*I1 - lambda2*R1 + gamma2p*IR2 + gamma12t*I12 - mu*R1==0; % R1

sol = solve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6],[S,I1,I2,I12,IR2,R1])

S = (b*(gamma1 + mu + rho*s1)*(gamma1*gamma2 + beta1*mu + gamma2*mu + rho^2*s1*s2 + gamma1*rho*s2 + gamma2*rho*s1 + mu*rho*s2))/(beta1^2*mu^2 - beta1*mu^3 + beta1*beta2*mu^2 - beta1*gamma1*mu^2 - beta1*mu^2*rho*s1 + beta1*beta2*gamma1*mu + beta1*beta2*mu*rho*s1);
I1 = -(b*(gamma1 - beta1 + mu + rho*s1)*(gamma1^2*gamma2 + beta1*mu^2 + gamma1^2*rho*s2 + gamma2*rho^2*s1^2 + rho^3*s1^2*s2 + beta2*gamma1*gamma2 + beta1*gamma1*mu + beta1*gamma2*mu + beta2*gamma2*mu + gamma1*gamma2*mu + beta2*gamma1*rho*s2 + beta2*gamma2*rho*s1 + 2*gamma1*gamma2*rho*s1 + beta1*mu*rho*s1 + beta1*mu*rho*s2 + beta2*mu*rho*s2 + gamma1*mu*rho*s2 + gamma2*mu*rho*s1 + beta2*rho^2*s1*s2 + 2*gamma1*rho^2*s1*s2 + mu*rho^2*s1*s2 + mu^2*rho*s1*s2 - beta2*rho^2*s1^2*s2 + mu*rho^2*s1^2*s2 - beta2*gamma1*rho*s1*s2 - beta1*mu*rho*s1*s2 - beta2*mu*rho*s1*s2 + gamma1*mu*rho*s1*s2))/((beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu - beta1*gamma1*mu + beta1*beta2*rho*s1 - beta1*mu*rho*s1)*(rho^2*s1^2 + beta2*gamma1 + beta2*mu + gamma1*mu + gamma1^2 - rho^2*s1^2*s2 + beta2*rho*s1 + 2*gamma1*rho*s1 + mu*rho*s1 - gamma1*rho*s1*s2 - mu*rho*s1*s2));
I2 = -(b*(gamma1 + mu + rho*s1)^2*(gamma2 - beta2 + mu + rho*s2))/(beta1^2*mu^2 - beta1*mu^3 + beta1*beta2*mu^2 - beta1*gamma1*mu^2 - beta1*mu^2*rho*s1 + beta1*beta2*gamma1*mu + beta1*beta2*mu*rho*s1);
I12 = (b*(gamma1 - beta1 + mu + rho*s1)*(gamma2 - beta2 + mu + rho*s2)*(rho^2*s1^2 + beta2*gamma1 + beta1*mu + beta2*mu + gamma1*mu + gamma1^2 + beta2*rho*s1 + 2*gamma1*rho*s1 + mu*rho*s1))/((beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu - beta1*gamma1*mu + beta1*beta2*rho*s1 - beta1*mu*rho*s1)*(rho^2*s1^2 + beta2*gamma1 + beta2*mu + gamma1*mu + gamma1^2 - rho^2*s1^2*s2 + beta2*rho*s1 + 2*gamma1*rho*s1 + mu*rho*s1 - gamma1*rho*s1*s2 - mu*rho*s1*s2));
IR2 = (b*(gamma1 - beta1 + mu + rho*s1)*(gamma2 - beta2 + mu + rho*s2)*(beta2*gamma1^3 + beta2^2*gamma1^2 + beta2*rho^3*s1^3 + beta2^2*rho^2*s1^2 + beta1*gamma1*mu^2 + beta1*gamma1^2*mu + beta2*gamma1^2*mu + beta2^2*gamma1*mu + 3*beta2*gamma1^2*rho*s1 + 2*beta2^2*gamma1*rho*s1 + beta1*mu^2*rho*s1 + beta2^2*mu*rho*s1 + 3*beta2*gamma1*rho^2*s1^2 + beta1*mu*rho^2*s1^2 + beta2*mu*rho^2*s1^2 - beta2*rho^3*s1^3*s2 + beta1*beta2*gamma1*mu - 2*beta2*gamma1*rho^2*s1^2*s2 - beta1*mu*rho^2*s1^2*s2 - 2*beta2*mu*rho^2*s1^2*s2 + beta1*beta2*mu*rho*s1 + 2*beta1*gamma1*mu*rho*s1 + 2*beta2*gamma1*mu*rho*s1 - beta2*gamma1^2*rho*s1*s2 - beta1*mu^2*rho*s1*s2 - beta2*mu^2*rho*s1*s2 - beta1*gamma1*mu*rho*s1*s2 - 2*beta2*gamma1*mu*rho*s1*s2))/(beta2*mu*(beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu - beta1*gamma1*mu + beta1*beta2*rho*s1 - beta1*mu*rho*s1)*(rho^2*s1^2 + beta2*gamma1 + beta2*mu + gamma1*mu + gamma1^2 - rho^2*s1^2*s2 + beta2*rho*s1 + 2*gamma1*rho*s1 + mu*rho*s1 - gamma1*rho*s1*s2 - mu*rho*s1*s2));
R1 = -(b*(gamma1 - beta1 + mu + rho*s1)*(beta2^2*gamma1^2*gamma2 + beta1*gamma1^2*mu^2 + beta2*gamma1^3*gamma2 + beta1*gamma1*mu^3 + beta2^2*gamma2*rho^2*s1^2 + beta1*mu^2*rho^2*s1^2 + beta2^2*rho^3*s1^2*s2 - beta2*rho^4*s1^3*s2^2 + beta1*gamma1*gamma2*mu^2 + beta1*gamma1^2*gamma2*mu + beta2*gamma1^2*gamma2*mu + beta2^2*gamma1*gamma2*mu + beta2*gamma1^3*rho*s2 + beta1*mu^3*rho*s1 + beta2^2*gamma1^2*rho*s2 + beta2*gamma2*rho^3*s1^3 + beta2*rho^4*s1^3*s2 + 3*beta2*gamma1*gamma2*rho^2*s1^2 + beta1*gamma2*mu*rho^2*s1^2 + beta2*gamma2*mu*rho^2*s1^2 + 3*beta2*gamma1^2*rho^2*s1*s2 + 2*beta2^2*gamma1*rho^2*s1*s2 + 3*beta2*gamma1*rho^3*s1^2*s2 - beta2*gamma2*rho^3*s1^3*s2 + beta1*mu^2*rho^2*s1*s2 + beta1*mu*rho^3*s1^2*s2 + beta2^2*mu*rho^2*s1*s2 + beta2^2*mu^2*rho*s1*s2 + beta2*mu*rho^3*s1^2*s2 + beta1*beta2*gamma1*gamma2*mu - beta2*gamma1^2*rho^2*s1*s2^2 - 2*beta2*gamma1*rho^3*s1^2*s2^2 - beta1*mu^2*rho^2*s1*s2^2 - beta1*mu^2*rho^2*s1^2*s2 - beta1*mu*rho^3*s1^2*s2^2 - beta2*mu^2*rho^2*s1*s2^2 - beta2*mu^2*rho^2*s1^2*s2 + beta2^2*mu*rho^2*s1^2*s2 - 2*beta2*mu*rho^3*s1^2*s2^2 + 3*beta2*gamma1^2*gamma2*rho*s1 + 2*beta2^2*gamma1*gamma2*rho*s1 + 2*beta1*gamma1*mu^2*rho*s1 + beta1*gamma1*mu^2*rho*s2 + beta1*gamma2*mu^2*rho*s1 + beta1*gamma1^2*mu*rho*s2 + beta2*gamma1^2*mu*rho*s2 + beta2^2*gamma1*mu*rho*s2 + beta2^2*gamma2*mu*rho*s1 - beta1*mu^3*rho*s1*s2 - beta2*mu^3*rho*s1*s2 - beta2*gamma1^2*gamma2*rho*s1*s2 + beta1*beta2*mu*rho^2*s1*s2 + beta1*beta2*mu^2*rho*s1*s2 + 2*beta1*gamma1*mu*rho^2*s1*s2 - beta1*gamma1*mu^2*rho*s1*s2 - beta1*gamma2*mu^2*rho*s1*s2 + 2*beta2*gamma1*mu*rho^2*s1*s2 - beta2*gamma1*mu^2*rho*s1*s2 + beta2^2*gamma1*mu*rho*s1*s2 - beta2*gamma2*mu^2*rho*s1*s2 - 2*beta2*gamma1*gamma2*rho^2*s1^2*s2 - beta1*gamma1*mu*rho^2*s1*s2^2 - beta1*gamma2*mu*rho^2*s1^2*s2 - 2*beta2*gamma1*mu*rho^2*s1*s2^2 - 2*beta2*gamma2*mu*rho^2*s1^2*s2 + beta1*beta2*gamma1*mu*rho*s2 + beta1*beta2*gamma2*mu*rho*s1 + 2*beta1*gamma1*gamma2*mu*rho*s1 + 2*beta2*gamma1*gamma2*mu*rho*s1 - beta1*gamma1*gamma2*mu*rho*s1*s2 - 2*beta2*gamma1*gamma2*mu*rho*s1*s2))/(beta2*mu*(beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu - beta1*gamma1*mu + beta1*beta2*rho*s1 - beta1*mu*rho*s1)*(rho^2*s1^2 + beta2*gamma1 + beta2*mu + gamma1*mu + gamma1^2 - rho^2*s1^2*s2 + beta2*rho*s1 + 2*gamma1*rho*s1 + mu*rho*s1 - gamma1*rho*s1*s2 - mu*rho*s1*s2));

[S,I1,I2,I12,IR2,R1]-T

%% utility on reprend les resultats matlab precedents pour l'instant
clear all; close all;

isbetween = 1; %rhohat between rhohat1 and rhohat2
maxIter = 1; iter=1;

while(isbetween & iter<=maxIter)
    [beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);
    s1=1;s2=1;
    R2=beta2/(gamma2+mu);R1=beta1/(gamma1+mu);
    alpha1=beta1/s1*(1-1/R1);
    alpha2=beta2/s2*(1-1/R2);
    
    c=0;
    k=0; U=[];prev=[];
    vecRho=0:(max(alpha1,alpha2)/1000):max(alpha1,alpha2);
    for rho=vecRho
        k=k+1;
        I1 = -(b*(gamma1 - beta1 + mu + rho*s1)*(gamma1^2*gamma2 + beta1*mu^2 + gamma1^2*rho*s2 + gamma2*rho^2*s1^2 + rho^3*s1^2*s2 + beta2*gamma1*gamma2 + beta1*gamma1*mu + beta1*gamma2*mu + beta2*gamma2*mu + gamma1*gamma2*mu + beta2*gamma1*rho*s2 + beta2*gamma2*rho*s1 + 2*gamma1*gamma2*rho*s1 + beta1*mu*rho*s1 + beta1*mu*rho*s2 + beta2*mu*rho*s2 + gamma1*mu*rho*s2 + gamma2*mu*rho*s1 + beta2*rho^2*s1*s2 + 2*gamma1*rho^2*s1*s2 + mu*rho^2*s1*s2 + mu^2*rho*s1*s2 - beta2*rho^2*s1^2*s2 + mu*rho^2*s1^2*s2 - beta2*gamma1*rho*s1*s2 - beta1*mu*rho*s1*s2 - beta2*mu*rho*s1*s2 + gamma1*mu*rho*s1*s2))/((beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu - beta1*gamma1*mu + beta1*beta2*rho*s1 - beta1*mu*rho*s1)*(rho^2*s1^2 + beta2*gamma1 + beta2*mu + gamma1*mu + gamma1^2 - rho^2*s1^2*s2 + beta2*rho*s1 + 2*gamma1*rho*s1 + mu*rho*s1 - gamma1*rho*s1*s2 - mu*rho*s1*s2));
        I2 = -(b*(gamma1 + mu + rho*s1)^2*(gamma2 - beta2 + mu + rho*s2))/(beta1^2*mu^2 - beta1*mu^3 + beta1*beta2*mu^2 - beta1*gamma1*mu^2 - beta1*mu^2*rho*s1 + beta1*beta2*gamma1*mu + beta1*beta2*mu*rho*s1);
        I12 = (b*(gamma1 - beta1 + mu + rho*s1)*(gamma2 - beta2 + mu + rho*s2)*(rho^2*s1^2 + beta2*gamma1 + beta1*mu + beta2*mu + gamma1*mu + gamma1^2 + beta2*rho*s1 + 2*gamma1*rho*s1 + mu*rho*s1))/((beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu - beta1*gamma1*mu + beta1*beta2*rho*s1 - beta1*mu*rho*s1)*(rho^2*s1^2 + beta2*gamma1 + beta2*mu + gamma1*mu + gamma1^2 - rho^2*s1^2*s2 + beta2*rho*s1 + 2*gamma1*rho*s1 + mu*rho*s1 - gamma1*rho*s1*s2 - mu*rho*s1*s2));
        IR2 = (b*(gamma1 - beta1 + mu + rho*s1)*(gamma2 - beta2 + mu + rho*s2)*(beta2*gamma1^3 + beta2^2*gamma1^2 + beta2*rho^3*s1^3 + beta2^2*rho^2*s1^2 + beta1*gamma1*mu^2 + beta1*gamma1^2*mu + beta2*gamma1^2*mu + beta2^2*gamma1*mu + 3*beta2*gamma1^2*rho*s1 + 2*beta2^2*gamma1*rho*s1 + beta1*mu^2*rho*s1 + beta2^2*mu*rho*s1 + 3*beta2*gamma1*rho^2*s1^2 + beta1*mu*rho^2*s1^2 + beta2*mu*rho^2*s1^2 - beta2*rho^3*s1^3*s2 + beta1*beta2*gamma1*mu - 2*beta2*gamma1*rho^2*s1^2*s2 - beta1*mu*rho^2*s1^2*s2 - 2*beta2*mu*rho^2*s1^2*s2 + beta1*beta2*mu*rho*s1 + 2*beta1*gamma1*mu*rho*s1 + 2*beta2*gamma1*mu*rho*s1 - beta2*gamma1^2*rho*s1*s2 - beta1*mu^2*rho*s1*s2 - beta2*mu^2*rho*s1*s2 - beta1*gamma1*mu*rho*s1*s2 - 2*beta2*gamma1*mu*rho*s1*s2))/(beta2*mu*(beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu - beta1*gamma1*mu + beta1*beta2*rho*s1 - beta1*mu*rho*s1)*(rho^2*s1^2 + beta2*gamma1 + beta2*mu + gamma1*mu + gamma1^2 - rho^2*s1^2*s2 + beta2*rho*s1 + 2*gamma1*rho*s1 + mu*rho*s1 - gamma1*rho*s1*s2 - mu*rho*s1*s2));
        prev(k) = (I1+I2+I12+IR2)/b*mu;
        U(k) = rho*(prev(k)+c);
    end
    
    
    [Umax,imax] = max(U);
    rhomaxU = vecRho(imax);
    
    rhohat2 = beta2/2/s2*(1-1/R2);
    rhohat1 = beta1/R1/s1*(sqrt(R1)-1);
    
    isbetween = (rhomaxU>=min(rhohat1,rhohat2)) & (rhomaxU<=max(rhohat1,rhohat2));
    iter=iter+1
end
P1 = mu/beta1.*(beta1./(s1*vecRho+gamma1+mu)-1);
U1 = vecRho.*(P1+c);
P2 = 1-(s2*vecRho+gamma2+mu)/beta2;
U2 = vecRho.*(P2+c);
title([{'Utility of the SIRxSIS model (strategy 2)'},...
    {['$\beta_1$=',num2str(beta1), ' $\beta_2$=',num2str(beta2),' $\gamma_1(0)$=',num2str(round(gamma1,2)),' $\gamma_2(0)$=',...
    num2str(round(gamma2,2)), ' $s_1$=', num2str(s1),' $s_2$=', num2str(s2), ' $\mu$=', num2str(mu), ' $\pi$=',  num2str(b),...
    ' $\mathtt R_1(0)$=' num2str(round(R1,2)), ' $\mathtt R_2(0)=$' num2str(round(R2,2))]},...
    {[' $\rho_1\prime$=' num2str(round(alpha1,2)),' $\rho_2\prime$=' num2str(round(alpha2,2)),...
    ' $\hat \rho$=' num2str(round(rhomaxU,2)),...
    ' $\hat \rho_{1}$=' num2str(round(rhohat1,2)),' $\hat \rho_{2}$=' num2str(round(rhohat2,2)),...
    ' $c$=',num2str(round(c,2)),...
    ]}],...
    'Interpreter','latex')
xlabel('$\rho$',"fontweight","bold",'Interpreter','latex')
ylabel('$U(\rho)$',"fontweight","bold",'Interpreter','latex')
hold on
plot(vecRho, max(U1,0))
plot(vecRho, max(U2,0))
plot(vecRho,max(max(max(U,U1),U2),0))

legend('U_1', 'U_2','U')

figure(2)
plot(vecRho,prev)

%% Finding hatRho
clear all;
syms beta1 beta2 gamma1 gamma2 s1 s2 rho mu b
s1=1;s2=1;
gamma2=1;
gamma1=gamma2;
beta2=1;
beta1=beta2;
mu=1;


gamma1p = gamma1+s1*rho;gamma2p=gamma2+s2*rho;
R1p=beta1./(gamma1p+mu); R2p=beta2./(gamma2p+mu);
lambda1 = mu*(R1p-1);lambda2=beta2.*(1-1./R2p);
gamma2t = gamma2p-s1*s2*rho;

N = b/mu;
P2 = 1-1./R2p;
I1ES12 = b/mu.*lambda1./R1p.*(gamma2t./(gamma1p+mu)+(lambda1+gamma2p+mu)./(lambda1+beta2))./(lambda2+gamma1p+mu+gamma2t);%strategy2
%%cas ou s1=s2=1
%I1ES12 = -(b*(gamma1 - beta1 + mu + rho)*(gamma1^2*gamma2 + beta1*mu^2 + 2*gamma1*rho^2 + gamma1^2*rho + gamma2*rho^2 + 2*mu*rho^2 + mu^2*rho + rho^3 + beta2*gamma1*gamma2 + beta1*gamma1*mu + beta1*gamma2*mu + beta2*gamma2*mu + beta2*gamma2*rho + gamma1*gamma2*mu + 2*gamma1*gamma2*rho + beta1*mu*rho + 2*gamma1*mu*rho + gamma2*mu*rho))/((beta2*gamma1 + beta2*mu + beta2*rho + gamma1*mu + gamma1*rho + gamma1^2)*(beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu + beta1*beta2*rho - beta1*gamma1*mu - beta1*mu*rho));

P12 = P2 + I1ES12./N;

U12 = rho.*P12;

dU12 = diff(U12,rho);
solve(dU12==0,rho)

%% Finding rhoHat with directly defining prevalence12
clear all;close all;
syms beta1 beta2 gamma1 gamma2 s1 s2 rho mu b
s1=1;
s2=1;
beta1=1;
gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;
gamma2t = gamma2p-s1*s2*rho;
R1p = beta1/(gamma1p+mu);
R2p = beta2/(gamma2p+mu);
Lambda1 = mu*(R1p-1);
Lambda2 = beta2*(1-1/R2p);

P12 = 1-1/R2p + Lambda1/R1p*((gamma2t/(gamma1p+mu)+(Lambda1+gamma2p+mu)/(Lambda1+beta2))/(Lambda2+gamma1p+mu+gamma2t));

U12 = rho.*P12;

dU12 = diff(U12,rho);
solve(dU12==0,rho)

%% argmax of rho*I1/N
clear all;close all;
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);
vecRho = -5:0.01:10;
vecu1 = [];
vecdu1=[];
vecdu11=[];


for rho=vecRho
    I1 = -(b*(gamma1 - beta1 + mu + rho*s1).*(gamma1^2*gamma2 + beta1*mu^2 + gamma1^2*rho*s2 + gamma2*rho^2*s1^2 + rho^3*s1^2*s2 + beta2*gamma1*gamma2 + beta1*gamma1*mu + beta1*gamma2*mu + beta2*gamma2*mu + gamma1*gamma2*mu + beta2*gamma1*rho*s2 + beta2*gamma2*rho*s1 + 2*gamma1*gamma2*rho*s1 + beta1*mu*rho*s1 + beta1*mu*rho*s2 + beta2*mu*rho*s2 + gamma1*mu*rho*s2 + gamma2*mu*rho*s1 + beta2*rho^2*s1*s2 + 2*gamma1*rho^2*s1*s2 + mu*rho^2*s1*s2 + mu^2*rho*s1*s2 - beta2*rho^2*s1^2*s2 + mu*rho^2*s1^2*s2 - beta2*gamma1*rho*s1*s2 - beta1*mu*rho*s1*s2 - beta2*mu*rho*s1*s2 + gamma1*mu*rho*s1*s2))/((beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu - beta1*gamma1*mu + beta1*beta2*rho*s1 - beta1*mu*rho*s1)*(rho^2*s1^2 + beta2*gamma1 + beta2*mu + gamma1*mu + gamma1^2 - rho^2*s1^2*s2 + beta2*rho*s1 + 2*gamma1*rho*s1 + mu*rho*s1 - gamma1*rho*s1*s2 - mu*rho*s1*s2));
    i1 = I1*mu/b;
    u1 = rho*i1;
    vecu1=[vecu1,u1];
    
    du1 = (mu*rho*(beta1*beta2 - beta1*mu)*(gamma1 - beta1 + mu + rho)*(gamma1^2*gamma2 + beta1*mu^2 + 2*gamma1*rho^2 + gamma1^2*rho + gamma2*rho^2 + 2*mu*rho^2 + mu^2*rho + rho^3 + beta2*gamma1*gamma2 + beta1*gamma1*mu + beta1*gamma2*mu + beta2*gamma2*mu + beta2*gamma2*rho + gamma1*gamma2*mu + 2*gamma1*gamma2*rho + beta1*mu*rho + 2*gamma1*mu*rho + gamma2*mu*rho))/((beta2*gamma1 + beta2*mu + beta2*rho + gamma1*mu + gamma1*rho + gamma1^2)*(beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu + beta1*beta2*rho - beta1*gamma1*mu - beta1*mu*rho)^2) - (mu*rho*(gamma1^2*gamma2 + beta1*mu^2 + 2*gamma1*rho^2 + gamma1^2*rho + gamma2*rho^2 + 2*mu*rho^2 + mu^2*rho + rho^3 + beta2*gamma1*gamma2 + beta1*gamma1*mu + beta1*gamma2*mu + beta2*gamma2*mu + beta2*gamma2*rho + gamma1*gamma2*mu + 2*gamma1*gamma2*rho + beta1*mu*rho + 2*gamma1*mu*rho + gamma2*mu*rho))/((beta2*gamma1 + beta2*mu + beta2*rho + gamma1*mu + gamma1*rho + gamma1^2)*(beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu + beta1*beta2*rho - beta1*gamma1*mu - beta1*mu*rho)) - (mu*rho*(gamma1 - beta1 + mu + rho)*(beta2*gamma2 + 2*gamma1*gamma2 + beta1*mu + 2*gamma1*mu + gamma2*mu + 4*gamma1*rho + 2*gamma2*rho + 4*mu*rho + gamma1^2 + mu^2 + 3*rho^2))/((beta2*gamma1 + beta2*mu + beta2*rho + gamma1*mu + gamma1*rho + gamma1^2)*(beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu + beta1*beta2*rho - beta1*gamma1*mu - beta1*mu*rho)) - (mu*(gamma1 - beta1 + mu + rho)*(gamma1^2*gamma2 + beta1*mu^2 + 2*gamma1*rho^2 + gamma1^2*rho + gamma2*rho^2 + 2*mu*rho^2 + mu^2*rho + rho^3 + beta2*gamma1*gamma2 + beta1*gamma1*mu + beta1*gamma2*mu + beta2*gamma2*mu + beta2*gamma2*rho + gamma1*gamma2*mu + 2*gamma1*gamma2*rho + beta1*mu*rho + 2*gamma1*mu*rho + gamma2*mu*rho))/((beta2*gamma1 + beta2*mu + beta2*rho + gamma1*mu + gamma1*rho + gamma1^2)*(beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu + beta1*beta2*rho - beta1*gamma1*mu - beta1*mu*rho)) + (mu*rho*(beta2 + gamma1)*(gamma1 - beta1 + mu + rho)*(gamma1^2*gamma2 + beta1*mu^2 + 2*gamma1*rho^2 + gamma1^2*rho + gamma2*rho^2 + 2*mu*rho^2 + mu^2*rho + rho^3 + beta2*gamma1*gamma2 + beta1*gamma1*mu + beta1*gamma2*mu + beta2*gamma2*mu + beta2*gamma2*rho + gamma1*gamma2*mu + 2*gamma1*gamma2*rho + beta1*mu*rho + 2*gamma1*mu*rho + gamma2*mu*rho))/((beta2*gamma1 + beta2*mu + beta2*rho + gamma1*mu + gamma1*rho + gamma1^2)^2*(beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu + beta1*beta2*rho - beta1*gamma1*mu - beta1*mu*rho));
    vecdu1=[vecdu1,du1];
    
    R1p = beta1/(gamma1+s1*rho+mu);
    Lambda1 = mu*(R1p-1);
    du11 = R1p^2*(Lambda1+beta2)^2*(gamma1+mu+s1*rho)^5*du1;
    vecdu11=[vecdu11,du11];
end
plot(vecRho,vecdu1)
hold on
plot(vecRho,vecdu11)
%%
clear all
syms beta1 beta2 gamma1 gamma2 s1 s2 rho mu b
s1=1;s2=1;
I1 = -(b*(gamma1 - beta1 + mu + rho*s1).*(gamma1^2*gamma2 + beta1*mu^2 + gamma1^2*rho*s2 + gamma2*rho^2*s1^2 + rho^3*s1^2*s2 + beta2*gamma1*gamma2 + beta1*gamma1*mu + beta1*gamma2*mu + beta2*gamma2*mu + gamma1*gamma2*mu + beta2*gamma1*rho*s2 + beta2*gamma2*rho*s1 + 2*gamma1*gamma2*rho*s1 + beta1*mu*rho*s1 + beta1*mu*rho*s2 + beta2*mu*rho*s2 + gamma1*mu*rho*s2 + gamma2*mu*rho*s1 + beta2*rho^2*s1*s2 + 2*gamma1*rho^2*s1*s2 + mu*rho^2*s1*s2 + mu^2*rho*s1*s2 - beta2*rho^2*s1^2*s2 + mu*rho^2*s1^2*s2 - beta2*gamma1*rho*s1*s2 - beta1*mu*rho*s1*s2 - beta2*mu*rho*s1*s2 + gamma1*mu*rho*s1*s2))/((beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu - beta1*gamma1*mu + beta1*beta2*rho*s1 - beta1*mu*rho*s1)*(rho^2*s1^2 + beta2*gamma1 + beta2*mu + gamma1*mu + gamma1^2 - rho^2*s1^2*s2 + beta2*rho*s1 + 2*gamma1*rho*s1 + mu*rho*s1 - gamma1*rho*s1*s2 - mu*rho*s1*s2));
i1 = I1*mu/b;
u1 = rho*i1;
du1=diff(u1,rho);

R1p = beta1/(gamma1+s1*rho+mu);
Lambda1 = mu*(R1p-1);
du11 = R1p^2*(Lambda1+beta2)^2*(gamma1+mu+s1*rho)^5*du1;
solve(du11==0,rho)


%% difference between I1 from strategy 1 and strategy 2 (13/10)
clear all;
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);
s1=1;s2=1;

gamma1p = gamma1+s1*rho;gamma2p=gamma2+s2*rho;
R1p = beta1./(gamma1p+mu); R2p=beta2./(gamma2p+mu);
lambda1 = mu*(R1p-1);lambda2=beta2.*(1-1./R2p);
gamma2t = gamma2p-s1*s2*rho;
I1ES12 = b/mu.*lambda1./R1p.*(gamma2t./(gamma1p+mu)+(lambda1+gamma2p+mu)./(lambda1+beta2))./(lambda2+gamma1p+mu+gamma2t);%strategy2

SES12_s1   = b/mu*(mu*(R1p-1) + gamma2p+mu)/(mu*(R1p-1)+beta2)/R1p;
I1ES12_s1  = b*(R1p-1)/(beta2+gamma1p)*(gamma2p/beta1+mu/b*SES12_s1);

I1ES12 - I1ES12_s1

%% argmax U - solution of the maximization problem (numerically, for a set of parameters)
clear all; close all;

rhobetweenrho1andrho2=true;
kMax=100; k=0; cond=1;

while rhobetweenrho1andrho2 && cond && k<kMax %cond:rho''/2>rhohat
    k=k+1
    [beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);
    s1=1;s2=1;
    clear rho;
    c=0;
    
    R1=beta1/(gamma1+mu); R2=beta2/(gamma2+mu);
    alpha1=beta1/s1*(1-1/R1); alpha2=beta2/s2*(1-1/R2);
    minalpha=min(alpha1,alpha2);
    maxalpha=max(alpha1,alpha2);
    
    fun = @(rho) U_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    %contraintes lineaires
    %A=[-1 1]'; b =[0 alpha]';
    %vecRhomax1(i) = fmincon(fun,alpha/2, A,b);
    %vecRhomax1(i) = min(max(fmincon(fun,0),0),alpha)
    rhohat = fminbnd(fun,0, maxalpha)
    
    rhohat2 = alpha2/2+c*beta2/(2*s2) %remplacer - par + ?
    rhohat1 = beta1*(sqrt(R1*mu/(mu-beta1*c))-1)/(R1*s1) %remplacer - par + ?
    
    (rhohat1+rhohat2)/2;
    
    vecRho = 0:minalpha/1000:maxalpha;
    U = max(-fun(vecRho),0);
    plot(vecRho, U)
    Umax=max(U);
    
    % conditions
    rhobetweenrho1andrho2 = (rhohat <= max(rhohat1,rhohat2) & rhohat>=min(rhohat1,rhohat2));
    
    [u0,i] = min(U(2:end));
    rhoprimeprime = vecRho(i) %rho such that U=0
    cond = rhoprimeprime/2>=rhohat;
end
%%
%%GRAPHIQUE
hold on
limit_y = 1.1*max(Umax);
xlim([0 maxalpha])
ylim([0 limit_y])
plot(alpha1*ones(1,100) , 0:limit_y/99:limit_y ,'--','Color',[0.6,0.6,0.6],'LineWidth',2)
plot(alpha2*ones(1,100) , 0:limit_y/99:limit_y ,'--','Color',[0.6,0.6,0.6],'LineWidth',2)

em = 0.05*Umax;
text(rhohat1,0,'$|$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
text(rhohat1,-em,'$\hat\rho_1$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
text(rhohat2,0,'$|$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
text(rhohat2,-em,'$\hat\rho_2$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
text(alpha1,0,'$|$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
text(alpha1,-em,'$\rho_1\prime$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
text(alpha2,0,'$|$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
text(alpha2,-em,'$\rho_2\prime$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')

text(rhohat,0,'$|$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
text(rhohat,-em,'$\hat\rho$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
text((rhohat1+rhohat2)/2,0,'$|$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')
text((rhohat1+rhohat2)/2,-em,'$\frac{\hat\rho_1 + \hat\rho_2}{2}$','Interpreter','latex','FontSize',20, 'FontWeight','bold','HorizontalAlignment', 'center')


%% trouver rho tel que U12=0 i.e. Pi12=0
clear all;
syms beta1 beta2 gamma1 gamma2 s1 s2 rho mu b
s1=1;s2=1;
%gamma1=0;
mu=0;

gamma1p = gamma1+s1*rho;gamma2p=gamma2+s2*rho;
R1p=beta1./(gamma1p+mu); R2p=beta2./(gamma2p+mu);
lambda1 = mu*(R1p-1);lambda2=beta2.*(1-1./R2p);
gamma2t = gamma2p-s1*s2*rho;

N = b/mu;
P2 = 1-1./R2p;
I1ES12 = b/mu.*lambda1./R1p.*(gamma2t./(gamma1p+mu)+(lambda1+gamma2p+mu)./(lambda1+beta2))./(lambda2+gamma1p+mu+gamma2t);%strategy2
%%cas ou s1=s2=1
%I1ES12 = -(b*(gamma1 - beta1 + mu + rho)*(gamma1^2*gamma2 + beta1*mu^2 + 2*gamma1*rho^2 + gamma1^2*rho + gamma2*rho^2 + 2*mu*rho^2 + mu^2*rho + rho^3 + beta2*gamma1*gamma2 + beta1*gamma1*mu + beta1*gamma2*mu + beta2*gamma2*mu + beta2*gamma2*rho + gamma1*gamma2*mu + 2*gamma1*gamma2*rho + beta1*mu*rho + 2*gamma1*mu*rho + gamma2*mu*rho))/((beta2*gamma1 + beta2*mu + beta2*rho + gamma1*mu + gamma1*rho + gamma1^2)*(beta1^2*mu - beta1*mu^2 + beta1*beta2*gamma1 + beta1*beta2*mu + beta1*beta2*rho - beta1*gamma1*mu - beta1*mu*rho));
P12 = P2 + I1ES12./N;

solve(P12==0,rho)
%-> s1=1=s2 : 4 roots (not sucessfully computed)
%-> s1=1=s2, gamma2=0 : 3 roots (not sucessfully computed)

%P12=0 if la quantite suivante vaut 0 (voir les calculs du 7 janvier)
%calculs a verifier
%a = gamma2*(mu*beta1-mu*(gamma1p+mu)+beta2)*(mu*(beta1-gamma1p-mu)) + (gamma1p+mu)*(mu*beta1-mu*(gamma1p+mu)+gamma2p+mu)*mu*(beta1-gamma1p-mu) + ...
%    beta1*(gamma1p+mu)*(beta2+gamma1)*(mu*beta1-(gamma1p+mu)*mu+beta2)*(1-1/R2p);

%% CONJECTURE rho"/2 = argmaxU12
clear all; close all;
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);
s1=1;s2=1;
R1=beta1/(gamma1+mu);
R2=beta2/(gamma2+mu);
alpha1=beta1/s1*(1-1/R1);
alpha2=beta2/s2*(1-1/R2);

% argmaxU12
maxalpha = max(alpha1,alpha2);
clear rho;
fun = @(rho) U_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,0);
argmaxU = fminbnd(fun,0,maxalpha)

%rho/U12=0 i.e. rho/P12=0
rhoU120 = fzero(fun,(alpha1+alpha2)/2)

%comparison between argmaxU12 and rho" (=rhoU120)
argmaxU*2


%% Finding c1 and c2
%c1 such that argmaxU12 = rho'(min) (i.e. alpha(min))
%c2 such that argmaxU12 = rho'(max) (i.e. alpha(max)
%c3 such that argmaxU12 = 0 (i.e. U<0 for all rho>0)













%% plot: rhohat en fonction de c
clear all;
close all;
condition =1
tmax = 1000
t=0;
while (condition && t<tmax)
    t=t+1
    alpha1=1;alpha2=0;
    %while(alpha2<alpha1)
        %[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);
        %beta1=3.44;beta2=1.71;gamma1=1.32;gamma2=0.13;s1=1;s2=1;mu=1.224;b=0.16;
        %beta1=0.42;beta2=2.18;gamma1=0.37;gamma2=1.71;s1=1;s2=1;mu=0.03; b=0.03; %s1=0.67;s2=0.42
        %Paris HSH 
        mu = 1/35; PHIV = 16.1/100*9/100;gamma1 = 1/2.7;beta1 = mu./(mu./(gamma1+mu)-PHIV);R1 = beta1/(gamma1+mu);s1=1;
        Pch = 5.3/100;R2 = 1/(1-Pch);gamma2 = 1/1.5;beta2=R2*(gamma2+mu);s2=1;
        b=10;
    
        R1=beta1/(gamma1+mu);
        R2=beta2/(gamma2+mu);
            
        %      R1 = 1 + abs(random(makedist('normal','mu',0,'sigma',1),1,1)*0.2/1.96); % 95% des valeurs de R1 sont comprises entre 1 et 1.10;
        %      R2 = R1*(1+abs(random(makedist('normal','mu',0,'sigma',1),1,1)*0.4/1.96));
        %      gamma1 = 1/2.7; %in years
        %      gamma2 = 1/(7/12); %environ 7 mois
        %      s1=1;s2=1;
        %      mu = 1/35;
        %      b = mu;
        %      beta1 = (gamma1+mu)*R1;
        %      beta2 = (gamma2+mu)*R2;
        
        alpha1 = beta1/s1*(1-1/R1);
        alpha2 = beta2/s2*(1-1/R2);
    %end
    maxalpha = max(alpha1,alpha2);
    minalpha = min(alpha1,alpha2);
    
    %c1, c=c12_1/argmaxU12=alpha1
    c_alpha1=(beta1*((gamma1 + mu)/beta1 - 1)*(s2/beta2 - (mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)*((s2 - s1*s2)/beta1 + ((s2 - (beta1*mu*s1)/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))^2)*(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1))) + (s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1))) + (mu*s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/((beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1))^2*(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)))))/(beta2 + gamma1 - beta1*((gamma1 + mu)/beta1 - 1) + beta1*s2*((gamma1 + mu)/beta1 - 1)) + (mu*(s1 - s1*s2)*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)*((gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta1 + ((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)))))/(beta2 + gamma1 - beta1*((gamma1 + mu)/beta1 - 1) + beta1*s2*((gamma1 + mu)/beta1 - 1))^2 + (beta1*mu*s1*((gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta1 + ((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)))))/((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))^2*(beta2 + gamma1 - beta1*((gamma1 + mu)/beta1 - 1) + beta1*s2*((gamma1 + mu)/beta1 - 1)))))/s1 - (gamma2 + mu - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta2 + (mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)*((gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta1 + ((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)))))/(beta2 + gamma1 - beta1*((gamma1 + mu)/beta1 - 1) + beta1*s2*((gamma1 + mu)/beta1 - 1)) + 1;
    %c2, c=c12_2/argmaxU12=alpha2
    c_alpha2=(mu*((gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1))/beta1 + ((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))))*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))/(beta2 + gamma1 + beta2*s1*((gamma2 + mu)/beta2 - 1) - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - (gamma2 + mu - beta2*((gamma2 + mu)/beta2 - 1))/beta2 + (beta2*((gamma2 + mu)/beta2 - 1)*(s2/beta2 - (mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1)*((s2 - s1*s2)/beta1 + ((s2 - (beta1*mu*s1)/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)^2)*(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))) + (s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))) + (mu*s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/((beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))^2*(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2))))/(beta2 + gamma1 + beta2*s1*((gamma2 + mu)/beta2 - 1) - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) + (mu*((gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1))/beta1 + ((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))))*(s1 - s1*s2)*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))/(beta2 + gamma1 + beta2*s1*((gamma2 + mu)/beta2 - 1) - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)^2 + (beta1*mu*s1*((gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1))/beta1 + ((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1)))))/((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)^2*(beta2 + gamma1 + beta2*s1*((gamma2 + mu)/beta2 - 1) - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2))))/s2 + 1;
    
    c12min = min(c_alpha1,c_alpha2);
    c12max = max(c_alpha1,c_alpha2);
    P20 = 1-1./R2;
    lambda10 = mu*(R1-1);
    i1ES120 = lambda10./(beta2+gamma1)*(gamma2/beta1+(lambda10+gamma2+mu)/R1/(lambda10+beta2));
    P120 = P20 + i1ES120;
    c120 = P120;
    
    i=1;
    vecC =(c12min-(c120-c12min)/2):(c120-c12min)/500:(c120+(c120-c12min)/2);
    vecRhomax1 = []; vecRhomax2 = [];vecRhomax3 = [];
    vecRhomaxU12=[];vecRhomaxU1=[];vecRhomaxU2=[];
    for c=vecC
        %fun1 = @(rho) U1_SIRSIS7(rho,beta1,gamma1,s1,mu,c);
        %rhomax1 = min(max(fmincon(fun1,minalpha),0),maxalpha);
        c11 = -mu/beta1*(1-1/R1); c10 = mu/beta1*(R1-1);
        c22 = 1/R2-1; %c such that argmax U2 = alpha2
        c20 = -c22; %c such that argmax U2 = 0
        
        rhomax1 = (beta1/(R1*s1)*(sqrt(R1*mu/(mu+beta1*c))-1))*(c>c11 & c<c10) + alpha1*(c<=c11);
        %fun2 = @(rho) U2_SIRSIS7(rho,beta2,gamma2,s2,mu,c);
        %rhomax2 = min(max(fmincon(fun2,minalpha),0),maxalpha);
        rhomax2 = beta2/(2*s2)*(1-1/R2-c)*(c>c22 & c<c20) + alpha2*(c<=c22);
        
        fun12 = @(rho) U12_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c); %-U12
        options = optimoptions('fmincon','Display','none');
        rhomax12 = fmincon(fun12,minalpha,[],[],[],[],[],[],[], options);
        
        %selon la deuxieme definition, voir les notes du 17 janvier 2022
        vecRhomax2(i) = rhomax12*(c<c120 & c>c12max) +...
            rhomax1*(c>c12min & c<=c12max & alpha1>alpha2) +...
            rhomax2*(c>c12min & c<=c12max & alpha1<alpha2) +...
            alpha1*(c<=c12min & alpha1>alpha2) +...
            alpha2*(c<=c12min & alpha1<alpha2); %manque le cas d'egalite entre  alpha1 et alpha2
        
        vecRhomaxU12(i) = rhomax12;
        vecRhomaxU1(i)  = rhomax1;
        vecRhomaxU2(i)  = rhomax2;
        
        % selon la premiere definition, voir les notes du 17 janvier 2022
        % ici on cherche argmaxU, donc ou bien le max est en rho12 ou en rho1/rho2
        Umaxrho1  = rhomax1*(mu/beta1*(beta1/(gamma1+s1*rhomax1+mu)-1)-c); %maxU1, i.e. U1 en hatrho1
        Umaxrho2  = rhomax2*(1-(gamma2+s2*rhomax2+mu)/beta2-c); %maxU2, i.e. U2 en hatrho2
        Umaxrho12 = -U12_SIRSIS7(rhomax12,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
        
        if(alpha1>alpha2)
            if(Umaxrho1<Umaxrho12 && rhomax1~=alpha1)
                vecRhomax1(i) = rhomax12;
            else
                vecRhomax1(i) = rhomax1;
            end
        else
            if(Umaxrho2<Umaxrho12 && rhomax2~=alpha2)
                vecRhomax1(i) = rhomax12;
            else
                vecRhomax1(i) = rhomax2;
            end
        end
        
        i=i+1;
    end
    
    
    %plot(vecC,vecRhomax1,vecC,vecRhomax2,vecC,vecRhomax3)
    legend('fmincon','fminbnd','max')
    %ylim([0,1.1*minalpha])
    xlim([(c12min-(c120-c12min)/2),(c120+(c120-c12min)/2)])
    
    % %thresholds cost from disease 1 and disease 2
    % c11=-mu/beta1*(1-1/R1);c10=mu/beta1*(R1-1);
    % c22=1/R2-1;c20=-c21;
    
    %recherche de c pour lequel Rhomax devient maxalpha
    vecTemp = (vecRhomax1-maxalpha==0);
    imaxalpha = sum(vecTemp); %index where rhomax=maxalpha
    cmaxalpha = vecC(imaxalpha);
    clear vecTemp ;
    
    %recherche de c pour lequel Rhomax devient rhohat1 ou rhohat2
    if(alpha1<alpha2)
        vecTemp = (vecRhomax1-vecRhomaxU2==0);
        iswitch = sum(vecTemp); %index where rhomax=rhohat2 (switch between rho2 and rho12)
        cswitch = vecC(iswitch);
        clear vecTemp ;
    end
%     if(alpha1>=alpha2)
%         vecTemp = (vecRhomax1-vecRhomaxU1==0);
%         iswitch = sum(vecTemp); %index where rhomax=rhohat2 (switch between rho2 and rho12)
%         cswitch = vecC(iswitch);
%         clear vecTemp ;
%     end
    %condition = (cswitch>c120/2);
    condition = 0;
end
%%
%plot
close all;
figure(1)
plt.rho = plot_paper2_procedure1(max(vecRhomax1,0), c12min, c12max, c120, alpha1,alpha2, vecC,c11,c10,c22,c20,cmaxalpha,cswitch)
hold on;
%legend('off')
plt.rho12 = plot(vecC,vecRhomaxU12,'b--');
plt.rho1 = plot(vecC,vecRhomaxU1,'r--');
%plt.argmaxU = plot(vecC,vecRhomax2,'g--');
plt.rho2 = plot(vecC,vecRhomaxU2,'y--')
%legend([plt.rho12,plt.rho1,plt.argmaxU],'$\hat\rho_{12}$','$\hat\rho_1$','first peak','Interpreter','latex')
legend([plt.rho,plt.rho12,plt.rho1,plt.rho2],'$\hat\rho(c)$','$\hat\rho_{12}(c)$','$\hat\rho_1(c)$','$\hat\rho_2(c)$','Interpreter','latex')

title([{'$\hat\rho(c)$'},...
    {['$\beta_1$=',num2str(round(beta1,2)), ' $\beta_2$=',num2str(round(beta2,2)),' $\gamma_1(0)$=',num2str(round(gamma1,2)),' $\gamma_2(0)$=',...
    num2str(round(gamma2,2)), ' $s_1$=', num2str(s1),' $s_2$=', num2str(s2), ' $\mu$=', num2str(round(mu,2)), ' $\pi$=',  num2str(round(b,2)),...
    ' $\mathtt R_1(0)$=' num2str(round(R1,2)), ' $\mathtt R_2(0)=$' num2str(round(R2,2))]},...
    {[' $\rho_1\prime$=' num2str(round(alpha1,2)),' $\rho_2\prime$=' num2str(round(alpha2,2)),...
    ' $c_0^{12}$=', num2str(round(c120,2)),' $c_1^{12}$=', num2str(round(c12min,2)),' $c_2^{12}$=', num2str(round(c12max,2)),...
    ' $c_0^1$=', num2str(round(c10,2)),' $c_1^1$=', num2str(round(c11,2)),' $c_0^2$=', num2str(round(c20,2)),' $c_2^2$=', num2str(round(c22,2)),' $c_s$=', num2str(round(cswitch,2))...
    ]}],...
    'Interpreter','latex')
if(0)
    figure(2)
    plot_paper2_procedure1(vecRhomax2, c1, c12max, c120, alpha1,alpha2, vecC,c11,c10,c21,c22)
    hold on;
    legend('off')
    plt.rho12 = plot(vecC,vecRhomaxU12,'b--');
    plt.rho1 = plot(vecC,vecRhomaxU1,'r--');
    plt.argmaxU = plot(vecC,vecRhomax1,'y--');
    %plot(vecC,vecRhomaxU2,'y')
    legend([plt.rho12,plt.rho1,plt.argmaxU],'$\hat\rho_{12}(c)$','$\hat\rho_1(c)$','argmax U','Interpreter','latex')
    
    title([{'$\hat\rho(c)$'},...
        {['$\beta_1$=',num2str(beta1), ' $\beta_2$=',num2str(beta2),' $\gamma_1(0)$=',num2str(round(gamma1,2)),' $\gamma_2(0)$=',...
        num2str(round(gamma2,2)), ' $s_1$=', num2str(s1),' $s_2$=', num2str(s2), ' $\mu$=', num2str(mu), ' $\pi$=',  num2str(b)]},...
        {[' $\mathtt R_1(0)$=' num2str(round(R1,2)), ' $\mathtt R_2(0)=$' num2str(round(R2,2)),...
        ' $\rho_1\prime$=' num2str(round(alpha1,2)),' $\rho_2\prime$=' num2str(round(alpha2,2)),...
        ' $c_0^{12}$=', num2str(round(c120,2)),' $c_1^{12}$=', num2str(round(c1,2)),' $c_2^{12}$=', num2str(round(c12max,2)),...
        ' $c_0^1$=', num2str(round(c10,2)),' $c_1^1$=', num2str(round(c11,2))...
        ]}],...
        'Interpreter','latex')
    
end

%%
figure(3)
cminalpha=c10;
plt.rho = plot_paper2_procedure1inv(max(vecRhomax1,0), c12min, c12max, c120, alpha1,alpha2, vecC,c11,c10,c22,c20,cmaxalpha,cminalpha,cswitch)
hold on;
%legend('off')
plt.rho12 = plot(vecRhomaxU12,vecC,'b--');
plt.rho1 = plot(vecRhomaxU1,vecC,'r--');
%plt.argmaxU = plot(vecC,vecRhomax2,'g--');
plt.rho2 = plot(vecRhomaxU2,vecC,'y--')
%legend([plt.rho12,plt.rho1,plt.argmaxU],'$\hat\rho_{12}$','$\hat\rho_1$','first peak','Interpreter','latex')
legend([plt.rho,plt.rho12,plt.rho1,plt.rho2],'$\hat\rho(c)$','$\hat\rho_{12}(c)$','$\hat\rho_1(c)$','$\hat\rho_2(c)$','Interpreter','latex')

title([{'$\hat\rho(c)$'},...
    {['$\beta_1$=',num2str(round(beta1,2)), ' $\beta_2$=',num2str(round(beta2,2)),' $\gamma_1(0)$=',num2str(round(gamma1,2)),' $\gamma_2(0)$=',...
    num2str(round(gamma2,2)), ' $s_1$=', num2str(s1),' $s_2$=', num2str(s2), ' $\mu$=', num2str(round(mu,2)), ' $\pi$=',  num2str(round(b,2)),...
    ' $\mathtt R_1(0)$=' num2str(round(R1,2)), ' $\mathtt R_2(0)=$' num2str(round(R2,2))]},...
    {[' $\rho_1\prime$=' num2str(round(alpha1,2)),' $\rho_2\prime$=' num2str(round(alpha2,2)),...
    ' $c_0^{12}$=', num2str(round(c120,2)),' $c_1^{12}$=', num2str(round(c12min,2)),' $c_2^{12}$=', num2str(round(c12max,2)),...
    ' $c_0^1$=', num2str(round(c10,2)),' $c_1^1$=', num2str(round(c11,2)),' $c_0^2$=', num2str(round(c20,2)),' $c_2^2$=', num2str(round(c22,2)),' $c_s$=', num2str(round(cswitch,2))...
    ]}],...
    'Interpreter','latex')
%% We add \rhohat1 and \rhohat2 of strategy 1
%rhohats
epsilon=0.;
vecC1 = epsilon.*vecC;
vecC2 = vecC-vecC1;
rho1max = 0.*(vecC1>=c12)+...
    (beta1.*(sqrt(R1*mu./(mu+beta1*vecC1))-1)./R1/s1).*(vecC1>c11 & vecC1<c12)+...
    alpha1.*(vecC1<=c11);
rho2max = 0.*(vecC2>=c22)+...
    (beta2./(2*s2).*(1-1/R2-vecC2)).*(vecC2>c21 & vecC2<c22)+...
    alpha2.*(vecC2<=c21);

plot(vecC,rho1max,'r','Linewidth', 2,'DisplayName','$\hat\rho_1(c_1)$')
plot(vecC,rho2max,'y','Linewidth', 2,'DisplayName','$\hat\rho_2(c_2)$')
hl = legend('show');
set(hl, 'Interpreter','latex')
title([{['$\hat\rho(c)$ and $\hat\rho_1(c_1),\hat\rho_2(c_2)$ where $c=c_1+c_2$ such that $c_1=\epsilon c$ and $c_2=(1-\epsilon) c$ with $\epsilon=$',num2str(epsilon)]},...
    {['$\beta_1$=',num2str(beta1), ' $\beta_2$=',num2str(beta2),' $\gamma_1(0)$=',num2str(round(gamma1,2)),' $\gamma_2(0)$=',...
    num2str(round(gamma2,2)), ' $s_1$=', num2str(s1),' $s_2$=', num2str(s2), ' $\mu$=', num2str(mu), ' $\pi$=',  num2str(b)]},...
    {[' $\mathtt R_1(0)$=' num2str(round(R1,2)), ' $\mathtt R_2(0)=$' num2str(round(R2,2)),...
    ' $\rho_1\prime$=' num2str(round(alpha1,2)),' $\rho_2\prime$=' num2str(round(alpha2,2)),...
    ' $c_0^{12}$=', num2str(round(c0,2)),' $c_1^{12}$=', num2str(round(c1,2)),' $c_2^{12}$=', num2str(round(c2,2)),...
    ' $c_0^1$=', num2str(round(c12,2)),' $c_1^1$=', num2str(round(c11,2))...
    ]}],...
    'Interpreter','latex')

%% plot U
vecC0 = [1.1*c1,c1,(c1+c11)/2,c2,(c11+c2)/2,c11,c2/2,0,c12/2,c12,(c12+c0)/2,c0,1.1*c0];
[vecC,isorted] = sort(vecC0);
vecCnames0 = {'1.1*c_1^{12}','c_1^{12}','(c_1^{12}+c_1^1)/2','c_2^{12}',...
    '(c_1^1+c_2^{12})/2','c_1^1','c_2^{12}/2','0','c_0^1/2','c_0^1','(c_0^1+c_0^{12})/2','c_0^{12}','1.1*c_0^{12}'};
vecCnames = vecCnames0(isorted);

figure()
for c=vecC
    vecRho=0:(maxalpha/100):maxalpha*1.1;
    %R1p = beta1./(gamma1+s1*vecRho+mu);R2p=beta2./(gamma2+s2*vecRho+mu);
    P1 = mu/beta1.*(beta1./(s1*vecRho+gamma1+mu)-1);
    U1 = vecRho.*(P1+c);
    P2 = 1-(s2*vecRho+gamma2+mu)/beta2;
    U2 = vecRho.*(P2+c);
    U = -U_SIRSIS7(vecRho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    plot(vecRho,max(U,0))
    hold on
end
legend(vecCnames)
text(alpha1,-0.01,'$\rho_1\prime$','Interpreter','latex','FontSize',22, 'FontWeight','bold','HorizontalAlignment', 'center')

% title([{'Utility of the SIRxSIS model (strategy 2)'},...
%             {['$\beta_1$=',num2str(beta1), ' $\beta_2$=',num2str(beta2),' $\gamma_1(0)$=',num2str(round(gamma1,2)),' $\gamma_2(0)$=',...
%             num2str(round(gamma2,2)), ' $s_1$=', num2str(s1),' $s_2$=', num2str(s2), ' $\mu$=', num2str(mu), ' $\pi$=',  num2str(b),...
%             ' $\mathtt R_1(0)$=' num2str(round(R1,2)), ' $\mathtt R_2(0)=$' num2str(round(R2,2))]},...
%             {[' $\rho_1\prime$=' num2str(round(alpha1,2)),' $\rho_2\prime$=' num2str(round(alpha2,2)),...
%             ' $\hat \rho$=' num2str(round(rhomaxU,2)),...
%             ' $\hat \rho_{1}$=' num2str(round(rhohat1,2)),' $\hat \rho_{2}$=' num2str(round(rhohat2,2)),...
%             ' $c$=',num2str(round(c,2)),...
%             ]}],...
%             'Interpreter','latex')
xlabel('$\rho$',"fontweight","bold",'Interpreter','latex')
ylabel('$U(\rho)$',"fontweight","bold",'Interpreter','latex')
%plot(vecRho, max(U1,0))
%plot(vecRho, max(U2,0))
%plot(vecRho,U,'Linewidth',2)

%legend('U_{12}','U_1', 'U_2','U')

%% Computing c2 and c3
clear all;close all;
syms beta1 beta2 gamma1 gamma2 s1 s2 rho mu b c

gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;
gamma2t = gamma2p-s1*s2*rho;
gamma1t = gamma1p-s1*s2*rho;

R1p = beta1/(gamma1p+mu);
R2p = beta2/(gamma2p+mu);
Lambda1 = mu*(R1p-1);
Lambda2 = beta2*(1-1/R2p);

R1=beta1/(gamma1+mu);
R2=beta2/(gamma2+mu);
alpha1 = beta1/s1*(1-1/R1);
alpha2 = beta2/s2*(1-1/R2);

P12 = 1-1/R2p + Lambda1/(beta2+gamma1t)*(gamma2t/beta1+(Lambda1+gamma2p+mu)/R1p/(Lambda1+beta2));

%
U120 = rho.*P12;
dU120 = diff(U120,rho);
rho=alpha2;
(mu*((gamma2 + rho*s2 - rho*s1*s2)/beta1 + ((gamma1 + mu + rho*s1)*(gamma2 + mu + rho*s2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1))))*(beta1/(gamma1 + mu + rho*s1) - 1))/(beta2 + gamma1 + rho*s1 - rho*s1*s2) - rho*(s2/beta2 - (mu*(beta1/(gamma1 + mu + rho*s1) - 1)*((s2 - s1*s2)/beta1 + (s1*(gamma2 + mu + rho*s2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1))) + ((s2 - (beta1*mu*s1)/(gamma1 + mu + rho*s1)^2)*(gamma1 + mu + rho*s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1))) + (mu*s1*(gamma2 + mu + rho*s2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))/((beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1))^2*(gamma1 + mu + rho*s1))))/(beta2 + gamma1 + rho*s1 - rho*s1*s2) + (mu*((gamma2 + rho*s2 - rho*s1*s2)/beta1 + ((gamma1 + mu + rho*s1)*(gamma2 + mu + rho*s2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1))))*(s1 - s1*s2)*(beta1/(gamma1 + mu + rho*s1) - 1))/(beta2 + gamma1 + rho*s1 - rho*s1*s2)^2 + (beta1*mu*s1*((gamma2 + rho*s2 - rho*s1*s2)/beta1 + ((gamma1 + mu + rho*s1)*(gamma2 + mu + rho*s2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu + rho*s1) - 1)))))/((gamma1 + mu + rho*s1)^2*(beta2 + gamma1 + rho*s1 - rho*s1*s2))) - (gamma2 + mu + rho*s2)/beta2 + 1


%% plot dynamically U in function of C
clear all;close all;
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);
R1=beta1/(gamma1+mu);
R2=beta2/(gamma2+mu);
alpha1=beta1/s1*(1-1/R1);
alpha2=beta2/s2*(1-1/R2);

maxalpha=max(alpha1,alpha2);
minalpha=min(alpha1,alpha2);

c_alpha1=(beta1*((gamma1 + mu)/beta1 - 1)*(s2/beta2 - (mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)*((s2 - s1*s2)/beta1 + ((s2 - (beta1*mu*s1)/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))^2)*(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1))) + (s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1))) + (mu*s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/((beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1))^2*(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)))))/(beta2 + gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1) + (mu*(s2 - s1*s2)*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)*((gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta1 + ((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)))))/(beta2 + gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)^2 + (beta1*mu*s1*((gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta1 + ((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)))))/((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))^2*(beta2 + gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))))/s1 - (gamma2 + mu - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta2 + (mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)*((gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta1 + ((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)))))/(beta2 + gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1) + 1;
c_alpha2=(beta2*((gamma2 + mu)/beta2 - 1)*(s2/beta2 - (mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1)*((s2 - s1*s2)/beta1 + ((s2 - (beta1*mu*s1)/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)^2)*(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))) + (s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))) + (mu*s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/((beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))^2*(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2))))/(beta2 + gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1)) + (mu*((gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1))/beta1 + ((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))))*(s2 - s1*s2)*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))/(beta2 + gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1))^2 + (beta1*mu*s1*((gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1))/beta1 + ((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1)))))/((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)^2*(beta2 + gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1)))))/s2 - (gamma2 + mu - beta2*((gamma2 + mu)/beta2 - 1))/beta2 + (mu*((gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1))/beta1 + ((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))))*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))/(beta2 + gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1)) + 1;
c1 = min(c_alpha1,c_alpha2);c3 = max(c_alpha1,c_alpha2);
P20 = 1-1./R2; lambda10 = mu*(R1-1); i1ES120 = lambda10./(beta2+gamma1)*(gamma2/beta1+(lambda10+gamma2+mu)/R1/(lambda10+beta2));P120 = P20 + i1ES120;
c2 = P120;

i=1;
vecC =(c1-(c2-c1)/2):(c2-c1)/10:(c2+(c2-c1)/2);
vecRho=0:(maxalpha/1000):maxalpha*1.1;
R1p = beta1./(gamma1+s1*vecRho+mu);R2p=beta2./(gamma2+s2*vecRho+mu);
for c=vecC
    figure(1)
    plot(vecRho,vecRho.*(max(mu/beta1.*(R1p-1),0)-c))
    hold on
    
    %ylim([0 1.5*max(1-1/R2,mu/beta1*(R1-1))])
    xlim([vecRho(1),vecRho(end)])
    plot(vecRho,vecRho.*(max(1-1./R2p,0)-c))
    [moinsU]=U12_SIRSIS7(vecRho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    plot(vecRho,moinsU)
    legend('U_1','U_2','U')
    title(['c=',num2str(round(c,2)), ' in the range [',num2str(round(vecC(1),3)),',',num2str(round(vecC(end),2)),']'])
    pause(0.1)
    %hold off
end


%% 07/04 : on cherche les possibilités de comaptibilité entre c11,c22 et alpha1, alpha2
tabcond=[];
for i=1:1000
    [beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);
    s1=1; s2=1;
    R1 = beta1/(gamma1+mu);
    R2 = beta2/(gamma2+mu);
    c11 = mu/beta1*(1/R1-1);
    c22 = 1/R2-1;
    alpha1 = beta1/s1*(1-1/R1);
    alpha2 = beta2/s2*(1-1/R2);
    cond1 = alpha1<alpha2;
    cond2=c11<c22;
    R2star = 1/(mu/beta1*(1/R1-1)+1);
    cond3 = R2>R2star;
    gamma2star = (R1-1)/(R2-1)*(gamma1+mu)-mu;
    cond4 = gamma2<gamma2star;
    tabcond = [tabcond;[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho,cond1,cond2,cond3,cond4]];
    if(~cond1&cond2)
        i
        %break;
       
    end

end

soustab = tabcond(tabcond(:,12)==0,:) 
t = num2cell(soustab(1,1:8));
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = deal(t{:});

%% zones combined testing 15/04/22

clear all;close all;
len=100;
gamma2=365/30; %IST
gamma1=365/14; %VIH
s1=1;s2=1;
mu=1/35;

R10 = linspace(1.1,5,len);
R20 = linspace(1.1,5,len);

[R1,R2] = meshgrid(R10,R20);

BETA1 = R1.*(gamma1+mu);
BETA2 = R2.*(gamma2+mu);
ALPHA1 = BETA1./s1.*(1-1./R1);
ALPHA2 = BETA2./s2.*(1-1./R2);
c11   = -mu./BETA1.*(1-1./R1); c01 = mu./BETA1.*(R1-1);
c22   = 1./R2-1; c02 = -c22;

%Z = ALPHA1>ALPHA2;
%plot(R10,yR2)
%surf(R1,R2,double(Z))
N = 2;
CEL1=[];
        
%calculer c1_12 et c2_12;
for i=1:(len)
    i
    for j=1:(len)
        c11_n = c11(i,j); c01_n = c01(i,j);
        c22_n = c22(i,j); c02_n = c02(i,j);
        mincnn = min(c11_n,c22_n); maxc0n = max(c01_n,c02_n);
        maxcnn = max(c11_n,c22_n);
        %interval of c
        vecC = linspace(1.5*mincnn,maxc0n,50);
        diff = 0;step=1;
        while (diff<step)
            
            step = vecC(2) - vecC(1);
            [tab,tabco,tabcn] = findRhohat(1,1,0,[BETA2(i,j),gamma2,1],[BETA1(i,j),gamma1,1],[],mu,5,vecC);        
            vecAlpha = [ALPHA2(i,j),ALPHA1(i,j)]; %attention ordre des alpha
            pres = 0.1*abs(vecC(1)-vecC(end));
            cs1dis = findThresholds(1,1,0, tab, vecAlpha, [BETA2(i,j),BETA1(i,j)], [gamma2,gamma1], [s2,s1], 5,mu,vecC);

            if(size(cs1dis,2)==0)%did not find cs1dis in the interval vecC
                %then cs1dis<min(cnn);
                cs1dis=vecC(1);
            end
            
            if(ALPHA1(i,j)<ALPHA2(i,j))
                %alpha1=alphaj %j : first disease eliminated
                cjj = c11_n;
            else
                cjj = c22_n;
            end
            diff = abs(cs1dis-cjj); %affiner vecC
            vecC = linspace(1.1*min(cs1dis,cjj),1.1*max(cs1dis,cjj),10);
            if(diff<step)
                break;
            end
       end
       CEL1(i,j) = cs1dis; %store costs for which 1 disease has been elminated

    end
end

yR2 = (gamma1+mu)/(gamma2+mu)*(R10-1)+1;
%%
%which disease elminates the other ?
% e.g. if c11<c112 and alpha1<alpha2 then combined testing eliminates disease 1
% more easily thanks to disease 2
close all;
n=15; % nb de subdvisions colormap
titre = ['IST : ', num2str(round(1/gamma2,2)),' years and VIH 2 : ', num2str(round(1/gamma1,2)),' years'];

dis2drives1 = (ALPHA1<=ALPHA2).*(c11-CEL1)./c11;%a voir
dis1drives2 = (ALPHA1>=ALPHA2).*(c22-CEL1)./c22;
tot = dis2drives1+dis1drives2; minv = min(min(tot));maxv=max(max(tot));
dis2drives1((ALPHA1>ALPHA2) & dis2drives1==0)=NaN;
dis1drives2((ALPHA1<ALPHA2) & dis1drives2==0)=NaN;

%coordinates to plot text
yR10end = (gamma1+mu)/(gamma2+mu)*(R10(end)-1)+1;
xR20end = (gamma2+mu)/(gamma1+mu)*(R20(end)-1)+1;
if(R10(end)>xR20end) %voir les notes du 14/04 (triangle superieur)
   xG1 = (min(R10)+min(R10)+xR20end)/3;%pas tout a fait ça avec R10(end)
   yG1 = (min(R20)+yR10end+max(R20))/3;
   xG2 = (min(R10)+max(R10)+max(R10)+xR20end)/4;
   yG2 = (min(R20)+min(R20)+max(R20)+max(R20))/4;
else
   xG1 = (min(R10)+min(R10)+R10(end)+R10(end))/4;%quadrilatere superieur
   yG1 = (min(R20)+max(R20)+max(R20)+yR10end)/4;
   xG2 = (min(R10)+max(R10)+max(R10))/3;
   yG2 = (min(R20)+min(R20)+yR10end)/3;
end

%plot : diff between cjj and cjij
figure(1) 
surf1=surf(R1,R2,tot);
surf1.EdgeColor = 'none';
maxv=max(max(tot));
%text(xG1,yG1,maxv,'\rho\prime_1<\rho\prime_2','HorizontalAlignment','center','FontSize',20)
%text(xG2,yG2,maxv,'\rho\prime_1>\rho\prime_2','HorizontalAlignment','center','FontSize',20)
new_cb = subdivisedColormap([[0.8,0.8,0.8];[1,1,1];[0.85,0.325,0.098]],n, 'quad'); %2^n+1
new_cb2 = new_cb(ceil((2^n+1)*0.5*(1-abs(minv)/maxv)):end,:);
colormap(new_cb2)
colorbar;
hold on;
plot3(R10,yR2,maxv*ones(length(R10),1),'-k');
ylim([1,max(R20)]);
title([{'Is combined testing better than specific testing ?'},...
    {'$(c_{j2} - c_{j1\times2})/c_{j2}$ such that $\rho_j\prime<\rho_i\prime$'},{titre}], 'Interpreter','latex')
xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')
%%
%plot : areas of alpha1 vs alpha2
figure(2)
surf2=surf(R1,R2,double(ALPHA1>ALPHA2));
surf2.EdgeColor = 'none';
text(xG1,yG1,4,'\rho\prime_1<\rho\prime_2','HorizontalAlignment','center','FontSize',20)
text(xG2,yG2,4,'\rho\prime_1>\rho\prime_2','HorizontalAlignment','center','FontSize',20);
colormap([[1,1,1];[0.8,0.8,0.8]]);
hold on;
plot3(R10,yR2,ones(length(R10),1),'-k');
ylim([1,max(R20)]);
title([{'$\rho_1$ vs $\rho_2$'},{titre}], 'Interpreter','latex')
xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')
%%
%plot : beta1 vs beta2
figure(3)
surf3 = surf(R1,R2,double(BETA1>BETA2));
surf3.EdgeColor = 'none';
hold on;
plot3(R10,yR2,maxv*ones(length(R10),1),'-k');
text(xG1,yG1,4,'$\beta_1<\beta_2$','HorizontalAlignment','center','FontSize',20,'Interpreter','latex','Color','w') %order beta a verif
text(xG2,yG2,4,'$\beta_1>\beta_2$','HorizontalAlignment','center','FontSize',20,'Interpreter','latex')
ylim([1,max(R20)]);
title([{'$\beta_1$ vs $\beta_2$'},{titre}], 'Interpreter','latex')
xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')
%%
%plot :  gradient of change
figure(4)
%colorbars
greyred = subdivisedColormap([[0.8,0.8,0.8];[1,1,1];[0.85,0.325,0.098]],n, 'quad'); %2^n+1
%put 0 as a median of the colorbar
amp = maxv-minv; %amplitude des valeurs
%maxv/amplitude %pourcentage de valeurs +
%minv/amplitude %pourcentage de valeurs -
greyred2 = greyred(ceil((2^n+1)*0.5*(1-abs(minv)/maxv)):end,:);

%surf1=surf(R1,R2,tot)
surf1=surf(R1,R2,dis2drives1)
surf1.EdgeColor = 'none';
colormap(greyred2)
cb = colorbar;
cb.YTick = [ceil(minv/100)*100 0 floor(maxv/10)*10];
cb.YTickLabel = {num2str(ceil(minv/100)*100), '0', num2str(floor(maxv/10)*10)};
freezeColors %adding to the second area where aplha1>alpha2
title([{'Is combined testing better than specific testing ?'},...
    {'$(c_{j2} - c_{j1\times2})/c_{j2}$ such that $\rho_j\prime<\rho_i\prime$'},{titre}], 'Interpreter','latex')
xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')
hold on
%%
surf2=surf(R1,R2,dis1drives2);
surf2.EdgeColor = 'none';
cb2 = colorbar('southoutside');
maxv_2 = max(max(dis1drives2));
minv_2 = min(min(dis1drives2));
greyyellow = subdivisedColormap([[0.8,0.8,0.8];[1,1,1];[0.929,0.8,0.125]],n, 'quad'); %2^n+1
%greyyellow2 = greyyellow(1:(end-ceil((2^n+1)*0.5*(1-maxv_2/abs(minv_2)))),:);
%colormap(greyyellow2)
%cb2.YTick = [round(minv_2,3), 0,  round(maxv_2*100)/100];
%cb2.YTickLabel = {[num2str(round(minv_2,2)*100),'%'], '0', [num2str(floor(maxv_2*100)),'%']};
hold on;
plot3(R10,yR2,maxv*ones(length(R10),1),'-k','LineWidth',1.5);
ylim([1,max(R20)]);

%%
figure(5)
tot2 = (ALPHA1<ALPHA2 & (C11<CEL1))*1 +...
    (ALPHA1<ALPHA2 & (C11>=CEL1))*2 +...
    (ALPHA1>=ALPHA2 & (C22<CEL1))*4 +...
    (ALPHA1>=ALPHA2 & (C22>CEL1))*3
surf4=surf(R1,R2,tot2);
surf4.EdgeColor = 'none';
new_cb2 = [[0.85,0.325,0.098];[0.8,0.8,0.8];[0.8,0.8,0.8];[0.9290,0.6940,0.125]]; %2^n+1
colormap(new_cb2)
%text(xG1,yG1,4,'\rho\prime_1<\rho\prime_2','HorizontalAlignment','center','FontSize',20)
%text(xG2,yG2,4,'\rho\prime_1>\rho\prime_2','HorizontalAlignment','center','FontSize',20);
hold on;
plot3(R10,yR2,4*ones(length(R10),1),'-k','LineWidth',1.5);
ylim([1,max(R20)]);
title([{'Is combined testing better than specific testing ?'},{titre}], 'Interpreter','latex')
xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')



%% Reprise, septembre 2024
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);

beta1=3*beta1;
beta1=beta1/5;
R1p=beta1/(gamma1+mu+s1*rho)%
R2p=beta2/(gamma2+mu+s2*rho)%

% Populations initiales
S0      = b/mu*1.1;
I10     = b/mu*0.1;
I20     = b/mu*0.1;
I30     = b/mu*0.01;    %I_12 (coinfection)
I40     = 0;
R10     = 0.;

% Parametres du systeme d'ODE
tspan = 0:0.1:10000;
Y0 = [S0; I10; I20; I30; I40; R10];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_V7(t,Y,b,beta1,beta2,s1,s2,rho,gamma1,gamma2,mu,'frequency'),tspan,Y0, options);
T = Ys(end,:)

R1     = beta1./(gamma1+mu);
R2     = beta2./(gamma2+mu);
alpha1 = beta1/s1*(1-1/R1);
alpha2 = beta2/s2*(1-1/R2);
gamma1p = gamma1+s1*rho;
gamma2p = gamma2+s2*rho;
gamma12t = s1*s2*rho;
gamma1t = gamma1p-s1*s2*rho;
gamma2t = gamma2p-s1*s2*rho;
R1p     = beta1./(gamma1p+mu);
R2p     = beta2./(gamma2p+mu);
Lambda1 = mu*(R1p-1);
SES12   = b/mu*(Lambda1 + gamma2p+mu)./(Lambda1+beta2)./R1p;
%I1ES12  = b/mu*Lambda1./(beta2+gamma1p-gamma12t).*(gamma2t./beta1 + mu/b.*SES12);
I1ES12  = b/mu*Lambda1./(beta2+gamma1t).*(gamma2t./beta1 + mu/b.*SES12);
%Lambda2 = beta2*(1-1/R2p)


%U_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c)
