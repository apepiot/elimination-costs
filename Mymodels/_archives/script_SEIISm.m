clear all;

%there is only one disease SICR and all the others are SIIS
m = 2; %nbdiseases SEIIS (max m=4)
syms mu Lambda1 Lambda2 Lambda3 Lambda4 nu1 nu2 nu3 nu4 eps1 eps2 eps3 eps4 sigma1 sigma2 sigma3 sigma4 gamma1 gamma2 gamma3 gamma4 b rho;

nDis = m;
dis = 1:nDis;
nbBoxesSEIIS = 4;

%parameters of the SEIIS
sigmas = [sigma1,sigma2,sigma3,sigma4];
gammas = [gamma1,gamma2,gamma3,gamma4];
Lambdas = [Lambda1,Lambda2,Lambda3,Lambda4];
nus = [nu1,nu2,nu3,nu4];
epss = [eps1,eps2,eps3,eps4];

nbCompartments = nbBoxesSEIIS^m;
% creating all the compartments
% a compartment is defined bcy a n-tuple (x1,x2,x3,...xn) 
% xi takes value in 1 for S,2 for I,3 for J and concerns infection i          
x = 1:nbBoxesSEIIS; %1:S,2:E,3:IA,4:IS

%disease 1 is a sicr, disease2 is a siiis, the m others are seiis
if m==1
   compartments = combPerso(x);
elseif m==2
    compartments = combPerso(x,x);
elseif m==3
    compartments = combPerso(x,x,x);
elseif m==4
    compartments = combPerso(x,x,x,x);
end

%matrix containing all the flow rates initialization
M = sym(zeros(nbCompartments) - mu*eye(nbCompartments));

for j=dis %selon la maladie j, on recupere la matrice associee   
    %if we look for progression of disease 3 and more (SIJS)    
    if  m==1
        otherDiseaseStatesConstant = NA;
    elseif m==2
        otherDiseaseStatesConstant = x';
    elseif m==3
        otherDiseaseStatesConstant = combPerso(x,x);
    elseif m==4
        otherDiseaseStatesConstant = combPerso(x,x,x);
    end

    mX = ODESEIISv2(0,0,0,0,Lambdas(j),epss(j),nus(j),gammas(j),sigmas(j),0,0,0);
    
    %il faut assigner cette matrice dans M, où 1 varie de 1 à 3 et oùxs
    %toutes les autres maladies sont dans un état constant (ex : I3,I13,I3J1)
    for k=otherDiseaseStatesConstant'
        index = find(sum(compartments(:,dis(dis~=j))==k',2)==(nDis-1));
        M(index, index) = M(index, index) + mX;
    end
end

% adding the rho rate to all the compartments except those where ppl are only infected and 
% symptomatics for the STI :chlamydia, gono,etc. (e.g. Jc, Jg, RJc, RJg)
indexS  = find(sum(compartments,2)==nDis);
%below : not infected by HIV (1) & not infected by syphilis (1) & (sympotmatic with
%the infection/s for which the indiviudal is infected (4))
%version 2
%indexInfSTI = find(sum(compartments(:,1:end),2)==4*m); 

%version 3
if m==0
    indexInfSTI=[];
elseif m==1
    indexInfSTI = find(compartments(:,1)==4);
elseif m==2
    indexInfSTI = find(compartments(:,1)==4 | compartments(:,2)==4);
end

for k=1:nbCompartments
    if (k~=indexS && ~ismember(k,indexInfSTI)) %if notS and if (infected by HIV or by syphilis or asymptomatic for another sti)
        M(k,k) = M(k,k) - rho;
        M(indexS,k) =  M(indexS,k) + rho;
    end
end

% matrix to ODE system
[C,dC,eqn] = matToODE(nbCompartments,M);
dC.'
eqn.'
[C;compartments'].'


%% verifications
%clear all;
[beta1,beta2,gamma1,gamma2,nu1,nu2,sigma1,sigma2,rho] = random_parameters(true, true);
eps1=0.5;eps2=0.2;
mu=1/35;b=5;%rho=0;
R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rho))./((mu + sigma1+rho).*(gamma1 + mu + nu1).*(mu + nu1 + rho)); %ok
R2 = (beta2*sigma2*(gamma2*(1-eps2) + mu + nu2 + eps2*rho))./((mu + sigma2+rho).*(gamma2 + mu + nu2).*(mu + nu2 + rho)); %ok

tspan = 0:.1:1000;
Y0 = (b/mu)/80*ones(320,1)';
options = odeset('AbsTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys1] = ode45(@(t,Y) ODE_SEIIS2_3(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                 beta2,nu2,eps2,sigma2,gamma2,mu,b,rho),tspan,ones(16,1), options);
T = Ys1(end,:);
%IST 1
indexS1 = find(compartments(:,1)==1); indexE1 = find(compartments(:,1)==2); indexIA1 = find(compartments(:,1)==3); indexIS1 = find(compartments(:,1)==4);
[sum(T(indexS1)),sum(T(indexE1)),sum(T(indexIA1)),sum(T(indexIS1))]
%comparaison avec les resultats theoriques
S1  = b/(mu*max(R1,1));
lambda1 = max(beta1*(R1-1)*(sigma1+rho+mu)/(beta1+(sigma1+rho+mu)*R1),0);
E1  = lambda1*S1/(sigma1+mu+rho);
IA1  = (1-eps1)*sigma1*E1/(nu1+rho+mu);
IS1  = eps1*sigma1*E1/(gamma1+nu1+mu);
[S1,E1,IA1,IS1] %OK

%IST 2
indexS2 = find(compartments(:,2)==1); indexE2 = find(compartments(:,2)==2); indexIA2 = find(compartments(:,2)==3); indexIS2 = find(compartments(:,2)==4);
[sum(T(indexS2)),sum(T(indexE2)),sum(T(indexIA2)),sum(T(indexIS2))]
%comparaison avec les resultats theoriques
S2  = b/(mu*max(R2,1));
lambda2 = max(beta2*(R2-1)*(sigma2+rho+mu)/(beta2+(sigma2+rho+mu)*R2),0);
E2  = lambda2*S2/(sigma2+mu+rho);
IA2  = (1-eps2)*sigma2*E2/(nu2+rho+mu);
IS2  = eps2*sigma2*E2/(gamma2+nu2+mu);
[S2,E2,IA2,IS2] %OK


%% comparaison des prevalences
vecRho=0:0.05:5; c=0;
souslaracine = @(gamma0,beta,nu,p,sigma) (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2);
alpha1 = max((beta1*eps1*sigma1 + sqrt(souslaracine(gamma1,beta1,nu1,eps1,sigma1)))/(2*(gamma1+mu+nu1)) - (2*mu+nu1+sigma1)/2,0);
alpha2 = (beta2*eps2*sigma2 + sqrt(souslaracine(gamma2,beta2,nu2,eps2,sigma2)))/(2*(gamma2+mu+nu2)) - (2*mu+nu2+sigma2)/2;

indexS2 = find(compartments(:,2)==1); indexE2 = find(compartments(:,2)==2); indexIA2 = find(compartments(:,2)==3); indexIS2 = find(compartments(:,2)==4);

P12=zeros(size(vecRho,2),1);
P1=zeros(size(vecRho,2),1);
P2=zeros(size(vecRho,2),1);
P12_2=zeros(size(vecRho,2),1);

i=1;
for rho=vecRho  
    [~,~,P12(i)] = U12_SEIIS2_v2(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,b,mu,rho,c,alpha1,alpha2);
    tspan = 0:1:500;
    Y0 = ones(16,1);
    options = odeset('RelTol',1e-3,'Stats','off');%,'OutputFcn',@odeplot);
    [~,ES] = ode45(@(t,Y) ODE_SEIIS2_3(t,Y,beta1,nu1,eps1,sigma1,gamma1,...
                                 beta2,nu2,eps2,sigma2,gamma2,mu,b,rho),tspan,Y0, options);
    T = ES(end,:);
    P12_2(i) = sum([sum(T(indexE2)),sum(T(indexIA2)),sum(T(indexIS2))])/(b/mu);
    %[~,~,P1(i)] = U_SEIISv2(eps1, beta1, sigma1, gamma1, mu, nu1, alpha2, c,alpha1);
    [~,~,P2(i)] = U_SEIISv2(eps2, beta2, sigma2, gamma2, mu, nu2, rho, c,alpha2);
    i=i+1;
end

plot(vecRho,P12,'LineWidth',2)
hold on
plot(vecRho,P2,'.','LineWidth',1.5)
plot(vecRho,P12_2,'LineWidth',1.,'Color','k')
legend('$\Pi_{12}$','$\Pi_2$','$\Pi_2$ from the coinf. mod.','Interpreter','latex')
xlabel('$\rho$','Interpreter','latex')
ylabel('prevalence')

%% continuity of U ?
clear all;
[beta1,beta2,gamma1,gamma2,nu1,nu2,sigma1,sigma2,~] = random_parameters(true, true);
eps1=0.5;eps2=0.2;
mu=1/35;b=5;rho=0;
R1 = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rho))./((mu + sigma1+rho).*(gamma1 + mu + nu1).*(mu + nu1 + rho)) %ok
R2 = (beta2*sigma2*(gamma2*(1-eps2) + mu + nu2 + eps2*rho))./((mu + sigma2+rho).*(gamma2 + mu + nu2).*(mu + nu2 + rho)) %ok
souslaracine = @(gamma0,beta,nu,p,sigma) (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2);
alpha1 = max((beta1*eps1*sigma1 + sqrt(souslaracine(gamma1,beta1,nu1,eps1,sigma1)))/(2*(gamma1+mu+nu1)) - (2*mu+nu1+sigma1)/2,0);
alpha2 = (beta2*eps2*sigma2 + sqrt(souslaracine(gamma2,beta2,nu2,eps2,sigma2)))/(2*(gamma2+mu+nu2)) - (2*mu+nu2+sigma2)/2;

vecC=0.4:0.1:0.8;    vecRho = 0:0.005:max(alpha1,alpha2)*1.1;
for c=vecC
    U = U_SEIIS2_v2(beta1,gamma1,nu1,sigma1,eps1,beta2,gamma2,nu2,sigma2,eps2,b,mu,vecRho,0.6,alpha1,alpha2);
    plot(vecRho,U)
    hold on
end