%% script HCV model (SEIR extend) = additionnal testing in I
close all; 
clear all;

[beta,zeta,gamma,sigma,s,omega,b,mu,rho] = random_parameters('false', 'false');
gammap=gamma+s*rho;

R0 = beta*(omega*sigma+mu+gamma)/((sigma+mu)*(mu+gamma))%pas rho contenu dedans % a verifier theoriquement

pop0=[0.8,0.01,0.01,0,0];
MaxTime = 1000;
[t, pop] = ode45(@(t,y) HCV2(t, y, beta, gammap,sigma,zeta,omega,b,mu,'frequency'),[0 MaxTime],pop0);
%b/mu
%sum(pop(end,:))
pop(end,:)

%% theoretically
R0 = beta*(omega*sigma+mu+gammap)/((sigma+gammap+mu)*(mu+gammap))
Sth = b/(mu*R0)
Ith = b/(sigma+mu+gammap)*(1-1/R0)
Cth = b/beta*(R0-1)*omega*sigma/(omega*sigma+mu+gammap) %idem HCV version 1
Ath = b/mu*((1-omega)*sigma/(sigma+mu+gammap)*(1-1/R0) + zeta*gammap*(R0-1)/(mu+zeta)/beta)
Tth = b*gammap/beta/(mu+zeta)*(R0-1)

%% utility
clear all; close all;
syms beta gamma s rho sigma omega mu b
gammap = gamma+s*rho;
Rp = beta.*(omega*sigma+mu+gammap)./((sigma+gammap+mu).*(mu+gammap))
R0 = beta.*(omega*sigma+mu+gamma)./((sigma+gamma+mu).*(mu+gamma))
%prevalence 
P = mu/beta.*(Rp-1);
c = 0;
U = rho.*(P+c);

dU = diff(U,rho);
solve(dU==0,rho)

%(mu*((beta*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)) - 1))/beta - (mu*rho*((beta*s*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)^2) - (beta*s)/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)) + (beta*s*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)^2*(gamma + mu + sigma + rho*s))))/beta
dU = mu/beta*(Rp-1 + s*rho/(gamma+mu + s*rho)*(beta*sigma*(1-omega)/(sigma+gamma+mu+s*rho)^2 -Rp))

a = gamma+mu;
dU*(a+s*rho)^2*(sigma+a+s*rho)^2
mu/beta*((a+s*rho)*beta*(omega*sigma+a+s*rho)*(sigma+a+s*rho) - (a+s*rho)^2*(sigma+a+s*rho)^2 + s*rho*beta*sigma*(1-omega)*(a+s*rho) - s*rho*beta*(omega*sigma+a+s*rho)*(sigma+a+s*rho))

solve(dU*(gamma+mu+s*rho)^2*(sigma+gamma+mu+s*rho)^2==0,rho)

rho = ((beta-2*(gamma+mu)-sigma)/2 + sqrt((beta-sigma)^2+4*beta*omega*sigma)/2)/s

dR = diff(Rp,rho)
%(beta*s)/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)) - (beta*s*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)*(gamma + mu + sigma + rho*s)^2) - (beta*s*(gamma + mu + omega*sigma + rho*s))/((gamma + mu + rho*s)^2*(gamma + mu + sigma + rho*s))
%s/(gamma+mu+s*rho)*(beta*sigma*(1-omega)/(sigma+gamma+mu+s*rho)^2 - Rp)

%% solution of the game
%if omega=1
rho = beta/R0*(sqrt(R0)-1)/s;

%% plot utility
alpha = ((beta-sigma)/2 - (gamma+mu) + sqrt(((beta-sigma)/2)^2 + beta*omega*sigma))/s;

rho=0:alpha/200:1000*alpha;
gammap = gamma+s*rho;
Rp = beta.*(omega*sigma+mu+gammap)./((sigma+gammap+mu).*(mu+gammap));
R0 = beta.*(omega*sigma+mu+gamma)./((sigma+gamma+mu).*(mu+gamma));


%prevalence 
P = mu/beta.*(Rp-1);
%c = - (mu*((beta*(beta/2 - sigma/2 + omega*sigma + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))) - 1))/beta - (mu*((beta*s*(beta/2 - sigma/2 + omega*sigma + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))^2) - (beta*s)/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))) + (beta*s*(beta/2 - sigma/2 + omega*sigma + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))^2*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))))*(gamma - beta/2 + mu + sigma/2 - ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/(beta*s)
c=mu/beta;
U = rho.*(P+c);
figure(1)
plot(rho,U)

%% on cherche a encadrer rhohat
close all; 
clear all;

t=0;Tmax = 100000; cond=1;%rhohat<alpha/2
while(t<Tmax & cond)
    [beta,zeta,gamma,sigma,s,omega,b,mu,rho] = random_parameters('false', 'false');
    R0 = beta*(omega*sigma+mu+gamma)/((sigma+mu+gamma)*(mu+gamma));%pas rho contenu dedans % a verifier theoriquement

    if (R0>1)
        clear rho;
        alpha = ((beta-sigma)/2 - (gamma+mu) + sqrt(((beta-sigma)/2)^2 + beta*omega*sigma))/s;

        rho=0:alpha/500:alpha;
        gammap = gamma+s*rho;
        Rp = beta.*(omega*sigma+mu+gammap)./((sigma+gammap+mu).*(mu+gammap));

        %utility
        P = mu/beta.*(Rp-1);
        U = rho.*P;

        %maxU
        [maxU,imax] = max(U); rhomax = rho(imax);
        
        rho=alpha;
        gammap = gamma+s*rho;
        Rp = beta.*(omega*sigma+mu+gammap)./((sigma+gammap+mu).*(mu+gammap));
        dU = mu/beta*(Rp-1 + s*rho/(gamma+mu + s*rho)*(beta*sigma*(1-omega)/(sigma+gamma+mu+s*rho)^2-Rp));

        %cond = rhomax<=max(alpha/2,0);
        cond = dU<0;
        t=t+1
    end
end


%% recherche de rhohat pour un set de parametres en fonction de c
clear all;close all;
%beta=0.0029*2.1+0.059*0.23; 
s=1;
gamma = 0.1*1/15; sigma = 52/8; omega=0.33; mu=1/75;
%[beta,zeta,gamma,sigma,s,omega,b,mu,rho] = random_parameters('true', 'true');

%R0 = beta.*(omega*sigma+mu+gamma)./((sigma+gamma+mu).*(mu+gamma));
R0 = 3.03;
beta = R0*((sigma+gamma+mu).*(mu+gamma))/(omega*sigma+mu+gamma);

%
alpha = ((beta-sigma)/2 - (gamma+mu) + sqrt(((beta-sigma)/2)^2 + beta*omega*sigma))/s;
% rho=0:0.1:100;
% moinsU = U_HCV2(rho, beta, gamma, s, sigma, omega, mu, c);
% plot(rho,moinsU)

% c/dU(alpha,c)==0
c2 = - (mu*((beta*(beta/2 - sigma/2 + omega*sigma + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))) - 1))/beta - (mu*((beta*s*(beta/2 - sigma/2 + omega*sigma + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))^2) - (beta*s)/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))) + (beta*s*(beta/2 - sigma/2 + omega*sigma + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))^2*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))))*(gamma - beta/2 + mu + sigma/2 - ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/(beta*s)
%attention ci dessous : alpha=alpha*s
% (mu*(1/(gamma+mu + sigma + alpha) - s/(alpha+gamma+mu+sigma*omega) + s/(alpha+gamma+mu))*(-alpha))/(beta*s)
c1 = mu/beta*(1-R0)

i=1; vecC =(c1-(c2-c1)/2):(c2-c1)/1000:(c2+(c2-c1)/2);
vecRhomax1 = []; vecRhomax2 = [];vecRhomax3 = [];
for c=vecC    
    fun = @(rho) U_HCV2(rho, beta, gamma, s, sigma, omega, mu, c);
    %contraintes lineaires
    %A=[-1 1]'; b =[0 alpha]';
    %vecRhomax1(i) = fmincon(fun,alpha/2, A,b);
    vecRhomax1(i) = min(max(fmincon(fun,0),0),alpha);
    vecRhomax2(i) = fminbnd(fun,0, alpha);
    
    %vieille methode
    vecRho = 0:0.0001:alpha;
    moinsU = U_HCV2(vecRho, beta, gamma, s, sigma, omega, mu, c);
    [Umax,imax] = max(-moinsU);
    vecRhomax3(i) = vecRho(imax);
    
    i=i+1; %c
end

plot(vecC,vecRhomax1,vecC,vecRhomax2,vecC,vecRhomax3)
legend('fmincon','fminbnd','max')
ylim([0,1.1*alpha])
xlim([(c1-(c2-c1)/2),(c2+(c2-c1)/2)])

%on ajoute rhohat du SIR
R0SIR = beta/(gamma+mu);c1SIR = (1-R0SIR)*mu/beta;c2SIR = mu/beta*(1-1/R0SIR);
alphaSIR = beta/s*(1-1./R0SIR);
rhohatSIR = (beta/s.*(sqrt(mu./((mu-beta*vecC).*R0SIR))-1./R0SIR)).*(vecC<c2SIR & vecC>c1SIR)+...
    alphaSIR.*(vecC>=c2SIR);
hold on;
plot(vecC,rhohatSIR)
legend('fmincon','fminbnd','max','SIR')

% plot of the results from the previous section %plot for the paper1 with nice blue area etc.
%plot_paper1_procedure(vecRhomax1, c1, c2, alpha, vecC)


%26/11 new plot U=rho(Pi-c)
figure()
plot_paper1_procedure2(vecRhomax1, -c2, -c1, alpha, -vecC)


%% on cherche rhohat par extrapolation de SIR
%syms rho
vecOmega = 0:0.001:1; vecDU=[];c=0;vecRhomax=[];i=1;
for omega=vecOmega
    R0 = beta.*(omega*sigma+mu+gamma)./((sigma+gamma+mu).*(mu+gamma));
    %rho = beta*(sqrt(R0)-1)/(R0*s);
    %gammap = gamma+s*rho;
    %Rp = beta.*(omega*sigma+mu+gammap)./((sigma+gammap+mu).*(mu+gammap));
    %vecDU = [vecDU, mu/beta*(Rp-1 + s*rho/(gamma+mu + s*rho)*(beta*sigma*(1-omega)/(sigma+gamma+mu+s*rho)^2 -Rp))];
    
    alpha = ((beta-sigma)/2 - (gamma+mu) + sqrt(((beta-sigma)/2)^2 + beta*omega*sigma))/s;
    
    vecRho = 0:alpha/1000:alpha;
    moinsU = U_HCV2(vecRho, beta, gamma, s, sigma, omega, mu, c);
    [Umax,imax] = max(-moinsU);
    vecRhomax(i) = vecRho(imax);
    i=i+1;
end

%plot(vecOmega,vecDU)
%%
close all;
plot(vecOmega,vecRhomax)
hold on;
R0 = beta.*(vecOmega*sigma+mu+gamma)./((sigma+gamma+mu).*(mu+gamma));
extrapRho = beta*(sqrt(R0 - vecOmega.*(1-vecOmega).*(vecOmega-w0))-1)./(R0*s);
plot(vecOmega, extrapRho)

figure(2) 
plot(vecOmega,vecRhomax./extrapRho)

w0 = (gamma+mu)/sigma*((sigma+gamma+mu)/beta -1);

%% comparaison entre SIR et HCV2
close all; clear all;
[beta,zeta,gamma,sigma,s,omega,b,mu,rho] = random_parameters('true', 'true');
R0 = beta.*(omega*sigma+mu+gamma)./((sigma+gamma+mu).*(mu+gamma))
%comparaison des c2
C2SIR = mu/beta*(1-(gamma+mu)/beta)
%C2HCV2 = - (mu*((beta*(beta/2 - sigma/2 + omega*sigma + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))) - 1))/beta - (mu*((beta*s*(beta/2 - sigma/2 + omega*sigma + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))^2) - (beta*s)/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))) + (beta*s*(beta/2 - sigma/2 + omega*sigma + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/((beta/2 - sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))^2*(beta/2 + sigma/2 + ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2))))*(gamma - beta/2 + mu + sigma/2 - ((beta/2 - sigma/2)^2 + beta*omega*sigma)^(1/2)))/(beta*s)

%solve(C2SIR-C2HCV2==0,omega)

omega=0:0.001:1;
R0 = beta.*(omega*sigma+mu+gamma)./((sigma+gamma+mu).*(mu+gamma));
alpha = ((beta-sigma)/2 - gamma-mu +sqrt(((beta-sigma)/2)^2 +beta*omega*sigma))/s;
gammap = gamma+s*alpha;
vecC2HCV2 = (mu/beta.*alpha.*(1./(gammap+mu+sigma) + 1./(gammap+mu) - 1./(gammap+mu+omega.*sigma)))*s;
plot(omega,vecC2HCV2)
hold on
plot(omega,mu/beta.*(1-1./R0))

figure()
plot(omega,vecC2HCV2-mu/beta.*(1-1./R0))

figure()
plot(omega, (gamma+mu+omega.*sigma)./(gammap+mu+omega.*sigma))

c2ter = mu/beta^2.*2.*s.*alpha./(gammap+mu+omega.*sigma).*(gammap+mu+omega.*sigma +(sigma-beta)/2-omega.*sigma)
figure()
c2quater = mu/beta.^2*2.*s.*alpha./(gammap+mu+omega.*sigma).*(((beta-sigma)/2)^2+beta.*omega.*sigma).^(1/2)
c2quint = mu/beta.^2*2.*s.*alpha.*(1+ ((sigma-beta)/2 - omega*sigma)/(gammap+mu+omega*sigma))
%% c2 < mu/beta
[beta,zeta,gamma,sigma,s,omega,b,mu,rho] = random_parameters('true', 'true');
alpha = ((beta-sigma)/2 - gamma-mu +sqrt(((beta-sigma)/2)^2 +beta*omega*sigma))/s;
gammap = gamma+s*alpha;
R0 = beta.*(omega*sigma+mu+gamma)./((sigma+gamma+mu).*(mu+gamma))
c2 = (mu/beta.*alpha.*(1./(gammap+mu+sigma) + 1./(gammap+mu) - 1./(gammap+mu+omega.*sigma)))*s
%c2bis = mu/beta*(gamma+mu+omega.*sigma)./(gammap+mu+omega.*sigma).*((s*alpha).^2./(beta.*(gamma+mu+omega.*sigma)) +1-1./R0)
mu/beta



%% STABILITE DU MODELE
clear all; close all;
syms beta omega sigma mu b gamma ksi gamma0 s

beta=5; omega=0.5; sigma=0.2; mu=0.01; b=0.6; ksi=0.2; gamma0=0.3; s=0.8; %pif

alpha = (beta-sigma)/2/s - (gamma0+mu)/s + sqrt((beta-sigma)^2/4+beta*omega*sigma)/s
R0 = beta*(omega*sigma+gamma0+mu)/((gamma0+mu)*(sigma+gamma0+mu))


gamma=gamma0+s*0.98*alpha;

%DFS
s0 =1; %S0/N0
i0 = 0; %I0/N0
c0=0;

%ES
Rp = beta*(omega*sigma+gamma+mu)/((gamma+mu)*(sigma+gamma+mu))
s0 = 1/Rp; %S0/N0
i0 = mu/(sigma+gamma+mu)*(1-1/Rp); %I0/N0
c0 = mu/beta*(Rp-1)*omega*sigma/(omega*sigma+gamma+mu); %C0/N0


M = [-(beta*(i0+c0)+mu), -beta*s0, -beta*s0,0,0;...
    beta*(i0+c0), beta*s0 - (sigma+gamma+mu), beta*s0,0,0;...
    0, omega*sigma, -(mu+gamma),0,0;...
    0, (1-omega)*sigma, 0, -mu, ksi;...
    0,gamma, gamma, 0, -(mu+ksi)]

eig(M)
det(M)


%% c1 en fonction de beta
clear all; close all;
[beta,zeta,gamma,sigma,s,omega,b,mu,rho] = random_parameters('true', 'true');
clear beta
beta = 0:0.1:100;
R = beta.*(omega*sigma+gamma+mu)/((sigma+gamma+mu)*(gamma+mu));
salpha = (beta-sigma)./2 - (gamma+mu) + sqrt(((beta-sigma)./2).^2 + beta.*omega*sigma);
gammap = gamma+salpha;
vecC1 = (-mu./beta*2.*salpha.*sqrt((beta-sigma).^2./4+beta.*omega.*sigma)./(beta.*(gammap+mu+omega*sigma))).*(R>1);

plot(beta,vecC1)
title([{'$c_1$ in function of $\beta$ when $\mathtt R(0)>1$'},...
            {[' $\gamma(0)$=',num2str(round(gamma,2)),' $\omega=$', num2str(round(omega,2)),' $\sigma=$', num2str(round(sigma,2)),...
            ' $s$=', num2str(s),' $\mu$=', num2str(mu), ' $\pi$=',  num2str(b)]}],...
            'Interpreter','latex')



