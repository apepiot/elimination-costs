% SEIISv2 (STI)
%par rapport à SEIIS : rho added from E to S
clear all
syms b lambda beta p sigma gamma0 rho mu nu S E IS IA

%beta=0.5;p=0.3;sigma=1/12;gamma0=6/7;rho=0.5;mu=1/35;nu=1/5; b=2; 

N = b/mu;
lambda = beta*(E+IA)/N;

dS  = b - lambda*S + rho*E + (gamma0+nu)*IS + (rho+nu)*IA - mu*S;
dE  = lambda*S - (sigma+rho)*E - mu*E;
dIA = (1-p)*sigma*E - (rho+nu)*IA - mu*IA; %3:IA
dIS = p*sigma*E - (gamma0+nu)*IS - mu*IS; %4:IS

eqn1 = dS ==0;
eqn2 = dE==0;
eqn3 = dIS==0;
eqn4 = dIA==0;
sol = solve([eqn1,eqn2,eqn3,eqn4], [S,E,IS,IA]);

%% Case lambda = beta*(IS+IA)/N;
%S 
Smat = (b*mu^3 + b*mu*rho^2 + 2*b*mu^2*rho + b*nu*rho^2 + b*nu^2*rho + b*mu^2*sigma + b*nu^2*sigma + b*gamma0*mu^2 + b*gamma0*rho^2 + b*mu*nu^2 + 2*b*mu^2*nu + b*gamma0*mu*nu + 2*b*gamma0*mu*rho + b*gamma0*nu*rho + b*gamma0*mu*sigma + b*gamma0*nu*sigma + b*gamma0*rho*sigma + 3*b*mu*nu*rho + 2*b*mu*nu*sigma + b*mu*rho*sigma + b*nu*rho*sigma)/(mu*sigma*(beta*gamma0 + beta*mu + beta*nu - beta*gamma0*p + beta*p*rho))

%E
Emat= -((mu + nu + rho)*(b*mu^4 + 2*b*mu^3*rho + b*nu^3*rho + b*mu^3*sigma + b*nu^3*sigma + b*gamma0^2*mu^2 + b*gamma0^2*rho^2 + 3*b*mu^2*nu^2 + b*mu^2*rho^2 + b*nu^2*rho^2 + 2*b*gamma0*mu^3 + b*mu*nu^3 + 3*b*mu^3*nu - b*beta*gamma0^2*sigma - b*beta*mu^2*sigma + 2*b*gamma0*mu*nu^2 + 4*b*gamma0*mu^2*nu + b*gamma0^2*mu*nu - b*beta*nu^2*sigma + 2*b*gamma0*mu*rho^2 + 4*b*gamma0*mu^2*rho + 2*b*gamma0^2*mu*rho + 2*b*gamma0*nu*rho^2 + 2*b*gamma0*nu^2*rho + b*gamma0^2*nu*rho + 2*b*gamma0*mu^2*sigma + b*gamma0^2*mu*sigma + 2*b*gamma0*nu^2*sigma + b*gamma0^2*nu*sigma + b*gamma0^2*rho*sigma + 2*b*mu*nu*rho^2 + 4*b*mu*nu^2*rho + 5*b*mu^2*nu*rho + 3*b*mu*nu^2*sigma + 3*b*mu^2*nu*sigma + b*mu^2*rho*sigma + b*nu^2*rho*sigma - 2*b*beta*gamma0*mu*sigma - 2*b*beta*gamma0*nu*sigma - 2*b*beta*mu*nu*sigma + 6*b*gamma0*mu*nu*rho + 4*b*gamma0*mu*nu*sigma + 2*b*gamma0*mu*rho*sigma + 2*b*gamma0*nu*rho*sigma + 2*b*mu*nu*rho*sigma + b*beta*gamma0^2*p*sigma + b*beta*gamma0*mu*p*sigma + b*beta*gamma0*nu*p*sigma - b*beta*gamma0*p*rho*sigma - b*beta*mu*p*rho*sigma - b*beta*nu*p*rho*sigma))/(sigma*(- beta*gamma0^2*mu^2*p + beta*gamma0^2*mu^2 - beta*gamma0^2*mu*nu*p + beta*gamma0^2*mu*nu + sigma*beta*gamma0^2*mu*p^2 - beta*gamma0^2*mu*p*rho - 2*sigma*beta*gamma0^2*mu*p + beta*gamma0^2*mu*rho + sigma*beta*gamma0^2*mu - beta*gamma0*mu^3*p + 2*beta*gamma0*mu^3 - 2*beta*gamma0*mu^2*nu*p + 4*beta*gamma0*mu^2*nu - 2*sigma*beta*gamma0*mu^2*p + 2*beta*gamma0*mu^2*rho + 2*sigma*beta*gamma0*mu^2 - beta*gamma0*mu*nu^2*p + 2*beta*gamma0*mu*nu^2 - 2*sigma*beta*gamma0*mu*nu*p + 2*beta*gamma0*mu*nu*rho + 2*sigma*beta*gamma0*mu*nu - 2*sigma*beta*gamma0*mu*p^2*rho + beta*gamma0*mu*p*rho^2 + 2*sigma*beta*gamma0*mu*p*rho + beta*mu^4 + 3*beta*mu^3*nu + beta*mu^3*p*rho + beta*mu^3*rho + sigma*beta*mu^3 + 3*beta*mu^2*nu^2 + 2*beta*mu^2*nu*p*rho + 2*beta*mu^2*nu*rho + 2*sigma*beta*mu^2*nu + beta*mu^2*p*rho^2 + 2*sigma*beta*mu^2*p*rho + beta*mu*nu^3 + beta*mu*nu^2*p*rho + beta*mu*nu^2*rho + sigma*beta*mu*nu^2 + beta*mu*nu*p*rho^2 + 2*sigma*beta*mu*nu*p*rho + sigma*beta*mu*p^2*rho^2))
%IA
IAmat = (b*gamma0^2*mu^2*p - b*gamma0^2*mu^2 + b*gamma0^2*mu*nu*p - b*gamma0^2*mu*nu + 2*b*gamma0^2*mu*p*rho + b*sigma*gamma0^2*mu*p - 2*b*gamma0^2*mu*rho - b*sigma*gamma0^2*mu + b*gamma0^2*nu*p*rho + b*sigma*gamma0^2*nu*p - b*gamma0^2*nu*rho - b*sigma*gamma0^2*nu + b*sigma*beta*gamma0^2*p^2 + b*gamma0^2*p*rho^2 + b*sigma*gamma0^2*p*rho - 2*b*sigma*beta*gamma0^2*p - b*gamma0^2*rho^2 - b*sigma*gamma0^2*rho + b*sigma*beta*gamma0^2 + 2*b*gamma0*mu^3*p - 2*b*gamma0*mu^3 + 4*b*gamma0*mu^2*nu*p - 4*b*gamma0*mu^2*nu + 4*b*gamma0*mu^2*p*rho + 2*b*sigma*gamma0*mu^2*p - 4*b*gamma0*mu^2*rho - 2*b*sigma*gamma0*mu^2 + 2*b*gamma0*mu*nu^2*p - 2*b*gamma0*mu*nu^2 + 6*b*gamma0*mu*nu*p*rho + 4*b*sigma*gamma0*mu*nu*p - 6*b*gamma0*mu*nu*rho - 4*b*sigma*gamma0*mu*nu + b*sigma*beta*gamma0*mu*p^2 + 2*b*gamma0*mu*p*rho^2 + 2*b*sigma*gamma0*mu*p*rho - 3*b*sigma*beta*gamma0*mu*p - 2*b*gamma0*mu*rho^2 - 2*b*sigma*gamma0*mu*rho + 2*b*sigma*beta*gamma0*mu + 2*b*gamma0*nu^2*p*rho + 2*b*sigma*gamma0*nu^2*p - 2*b*gamma0*nu^2*rho - 2*b*sigma*gamma0*nu^2 + b*sigma*beta*gamma0*nu*p^2 + 2*b*gamma0*nu*p*rho^2 + 2*b*sigma*gamma0*nu*p*rho - 3*b*sigma*beta*gamma0*nu*p - 2*b*gamma0*nu*rho^2 - 2*b*sigma*gamma0*nu*rho + 2*b*sigma*beta*gamma0*nu - b*sigma*beta*gamma0*p^2*rho + b*sigma*beta*gamma0*p*rho + b*mu^4*p - b*mu^4 + 3*b*mu^3*nu*p - 3*b*mu^3*nu + 2*b*mu^3*p*rho + b*sigma*mu^3*p - 2*b*mu^3*rho - b*sigma*mu^3 + 3*b*mu^2*nu^2*p - 3*b*mu^2*nu^2 + 5*b*mu^2*nu*p*rho + 3*b*sigma*mu^2*nu*p - 5*b*mu^2*nu*rho - 3*b*sigma*mu^2*nu + b*mu^2*p*rho^2 + b*sigma*mu^2*p*rho - b*sigma*beta*mu^2*p - b*mu^2*rho^2 - b*sigma*mu^2*rho + b*sigma*beta*mu^2 + b*mu*nu^3*p - b*mu*nu^3 + 4*b*mu*nu^2*p*rho + 3*b*sigma*mu*nu^2*p - 4*b*mu*nu^2*rho - 3*b*sigma*mu*nu^2 + 2*b*mu*nu*p*rho^2 + 2*b*sigma*mu*nu*p*rho - 2*b*sigma*beta*mu*nu*p - 2*b*mu*nu*rho^2 - 2*b*sigma*mu*nu*rho + 2*b*sigma*beta*mu*nu - b*sigma*beta*mu*p^2*rho + b*sigma*beta*mu*p*rho + b*nu^3*p*rho + b*sigma*nu^3*p - b*nu^3*rho - b*sigma*nu^3 + b*nu^2*p*rho^2 + b*sigma*nu^2*p*rho - b*sigma*beta*nu^2*p - b*nu^2*rho^2 - b*sigma*nu^2*rho + b*sigma*beta*nu^2 - b*sigma*beta*nu*p^2*rho + b*sigma*beta*nu*p*rho)/(- beta*gamma0^2*mu^2*p + beta*gamma0^2*mu^2 - beta*gamma0^2*mu*nu*p + beta*gamma0^2*mu*nu + sigma*beta*gamma0^2*mu*p^2 - beta*gamma0^2*mu*p*rho - 2*sigma*beta*gamma0^2*mu*p + beta*gamma0^2*mu*rho + sigma*beta*gamma0^2*mu - beta*gamma0*mu^3*p + 2*beta*gamma0*mu^3 - 2*beta*gamma0*mu^2*nu*p + 4*beta*gamma0*mu^2*nu - 2*sigma*beta*gamma0*mu^2*p + 2*beta*gamma0*mu^2*rho + 2*sigma*beta*gamma0*mu^2 - beta*gamma0*mu*nu^2*p + 2*beta*gamma0*mu*nu^2 - 2*sigma*beta*gamma0*mu*nu*p + 2*beta*gamma0*mu*nu*rho + 2*sigma*beta*gamma0*mu*nu - 2*sigma*beta*gamma0*mu*p^2*rho + beta*gamma0*mu*p*rho^2 + 2*sigma*beta*gamma0*mu*p*rho + beta*mu^4 + 3*beta*mu^3*nu + beta*mu^3*p*rho + beta*mu^3*rho + sigma*beta*mu^3 + 3*beta*mu^2*nu^2 + 2*beta*mu^2*nu*p*rho + 2*beta*mu^2*nu*rho + 2*sigma*beta*mu^2*nu + beta*mu^2*p*rho^2 + 2*sigma*beta*mu^2*p*rho + beta*mu*nu^3 + beta*mu*nu^2*p*rho + beta*mu*nu^2*rho + sigma*beta*mu*nu^2 + beta*mu*nu*p*rho^2 + 2*sigma*beta*mu*nu*p*rho + sigma*beta*mu*p^2*rho^2)
%IS
ISmat = -(p*(mu + nu + rho)*(b*mu^3 + b*mu*rho^2 + 2*b*mu^2*rho + b*nu*rho^2 + b*nu^2*rho + b*mu^2*sigma + b*nu^2*sigma + b*gamma0*mu^2 + b*gamma0*rho^2 + b*mu*nu^2 + 2*b*mu^2*nu - b*beta*mu*sigma + b*gamma0*mu*nu - b*beta*nu*sigma + 2*b*gamma0*mu*rho + b*gamma0*nu*rho + b*gamma0*mu*sigma + b*gamma0*nu*sigma + b*gamma0*rho*sigma + 3*b*mu*nu*rho + 2*b*mu*nu*sigma + b*mu*rho*sigma + b*nu*rho*sigma - b*beta*gamma0*sigma + b*beta*gamma0*p*sigma - b*beta*p*rho*sigma))/(- beta*gamma0^2*mu^2*p + beta*gamma0^2*mu^2 - beta*gamma0^2*mu*nu*p + beta*gamma0^2*mu*nu + sigma*beta*gamma0^2*mu*p^2 - beta*gamma0^2*mu*p*rho - 2*sigma*beta*gamma0^2*mu*p + beta*gamma0^2*mu*rho + sigma*beta*gamma0^2*mu - beta*gamma0*mu^3*p + 2*beta*gamma0*mu^3 - 2*beta*gamma0*mu^2*nu*p + 4*beta*gamma0*mu^2*nu - 2*sigma*beta*gamma0*mu^2*p + 2*beta*gamma0*mu^2*rho + 2*sigma*beta*gamma0*mu^2 - beta*gamma0*mu*nu^2*p + 2*beta*gamma0*mu*nu^2 - 2*sigma*beta*gamma0*mu*nu*p + 2*beta*gamma0*mu*nu*rho + 2*sigma*beta*gamma0*mu*nu - 2*sigma*beta*gamma0*mu*p^2*rho + beta*gamma0*mu*p*rho^2 + 2*sigma*beta*gamma0*mu*p*rho + beta*mu^4 + 3*beta*mu^3*nu + beta*mu^3*p*rho + beta*mu^3*rho + sigma*beta*mu^3 + 3*beta*mu^2*nu^2 + 2*beta*mu^2*nu*p*rho + 2*beta*mu^2*nu*rho + 2*sigma*beta*mu^2*nu + beta*mu^2*p*rho^2 + 2*sigma*beta*mu^2*p*rho + beta*mu*nu^3 + beta*mu*nu^2*p*rho + beta*mu*nu^2*rho + sigma*beta*mu*nu^2 + beta*mu*nu*p*rho^2 + 2*sigma*beta*mu*nu*p*rho + sigma*beta*mu*p^2*rho^2)

%sigma=1/10;beta=2;gamma0=3/10; mu=1/35;nu=1/2;p=0.8;rho=0;b=3;
lambdamat = beta*simplify(sol.IA+sol.IS)*mu/b;
lambdamat = -(gamma0*mu^2 + gamma0*rho^2 + mu*nu^2 + 2*mu^2*nu + mu*rho^2 + 2*mu^2*rho + nu*rho^2 + nu^2*rho + mu^2*sigma + nu^2*sigma + mu^3 - beta*gamma0*sigma - beta*mu*sigma + gamma0*mu*nu - beta*nu*sigma + 2*gamma0*mu*rho + gamma0*nu*rho + gamma0*mu*sigma + gamma0*nu*sigma + gamma0*rho*sigma + 3*mu*nu*rho + 2*mu*nu*sigma + mu*rho*sigma + nu*rho*sigma + beta*gamma0*p*sigma - beta*p*rho*sigma)/(gamma0*mu + gamma0*nu + gamma0*rho + gamma0*sigma + 2*mu*nu + mu*rho + nu*rho + mu*sigma + nu*sigma + mu^2 + nu^2 - gamma0*p*sigma + p*rho*sigma);

Rp = (beta*sigma*(gamma0*(1-p) + mu + nu + p*rho))./((mu + sigma+rho).*(gamma0 + mu + nu).*(mu + nu + rho)); %ok
S  = b/(mu*Rp); %ok
E = lambdamat*S/(sigma+mu+rho);
E = beta*(Rp-1)/(beta+(sigma+rho+mu)*Rp)*b/(mu*Rp)
lambda = (Rp-1).*beta.*(sigma+rho+mu)./(beta+(sigma+rho+mu).*Rp);

%rho'
alphamat1 = -(2*gamma0*mu + gamma0*nu + gamma0*sigma + 3*mu*nu + mu*sigma + nu*sigma + (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2)^(1/2) + 2*mu^2 + nu^2 - beta*p*sigma)/(2*(gamma0 + mu + nu))
souslaracine = (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2)
alpha1 = -((2*mu + nu + sigma)*(gamma0 + mu + nu) - beta*p*sigma + sqrt(souslaracine))/(2*(gamma0 + mu + nu));
alpha1 = (beta*p*sigma - sqrt(souslaracine))/(2*(gamma0+mu+nu)) - (2*mu+nu+sigma)/2;
alphamat2 = -(2*gamma0*mu + gamma0*nu + gamma0*sigma + 3*mu*nu + mu*sigma + nu*sigma - (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2)^(1/2) + 2*mu^2 + nu^2 - beta*p*sigma)/(2*(gamma0 + mu + nu))
alpha2 = (beta*p*sigma + sqrt(souslaracine))/(2*(gamma0+mu+nu)) - (2*mu+nu+sigma)/2;

P=1-1./Rp;

%utility function
syms c rho
U = rho*(P-c);
dU = diff(U,rho);
dU = 1 - rho*(((gamma0 + mu + nu)*(mu + nu + rho))/(beta*sigma*(mu + nu + p*rho - gamma0*(p - 1))) + ((gamma0 + mu + nu)*(mu + rho + sigma))/(beta*sigma*(mu + nu + p*rho - gamma0*(p - 1))) - (p*(gamma0 + mu + nu)*(mu + nu + rho)*(mu + rho + sigma))/(beta*sigma*(mu + nu + p*rho - gamma0*(p - 1))^2)) - ((gamma0 + mu + nu)*(mu + nu + rho)*(mu + rho + sigma))/(beta*sigma*(mu + nu + p*rho - gamma0*(p - 1))) - c
solve(dU==0,rho)

%%
clear all
[beta1,beta2,gamma1,gamma2,nu1,nu2,sigma1,mu,rho] = random_parameters(true, true);
p=0.5; b=2; beta1=beta1*100;
tspan=1:1000;
options = odeset('AbsTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys1] = ode45(@(t,Y) ODE_SEIISv2(t,Y,b,beta1,nu1,eps1,sigma1,gamma1,rho,mu),tspan,[1,1,1,1], options);
Ys1(end,:)

Rp = (beta1*sigma1*(gamma1*(1-eps1) + mu + nu1 + eps1*rho))./((mu + sigma1+rho).*(gamma1 + mu + nu1).*(mu + nu1 + rho)); %ok
S  = b/(mu*max(Rp,1))

%% numerically
clear all;
sigma=1/10;beta=2;gamma0=3/10; mu=1/35;nu=1/2;p=0.8;rho=0;b=0.5;
R0 = (beta*sigma*(gamma0*(1-p) + mu + nu ))./((mu + sigma).*(gamma0 + mu + nu).*(mu + nu)); %ok


souslaracine = (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2);
alpha = (beta*p*sigma + sqrt(souslaracine))/(2*(gamma0+mu+nu)) - (2*mu+nu+sigma)/2;

[U0,c0] = U_SEIISv2(p, beta, sigma, gamma0, mu, nu, 0, 0); %c2 = c0
[U1,c1] = U_SEIISv2(p, beta, sigma, gamma0, mu, nu, alpha, 0); %c1=c1

options = optimset('Display','off'); %options for minsearch
i=1; 
vecC =(c1-(c0-c1)/2):(c0-c1)/1000:(c0+(c0-c1)/2);
vecRhomax1 = zeros(1,length(vecC));
for c=vecC    
    fun = @(rho) -U_SEIISv2(p, beta, sigma, gamma0, mu, nu, rho, c);
    vecRhomax1(i) = min(max(fminsearch(fun,0,options),0),alpha);
    i=i+1; %c
end

%plot(vecC,vecRhomax1)
plot_paper1_procedure2(vecRhomax1, c1, c0, alpha, vecC)

%ratio tests/infections
clear rho
Rp = (beta*sigma*(gamma0*(1-p) + mu + nu + p*vecRhomax1))./((mu + sigma+vecRhomax1).*(gamma0 + mu + nu).*(mu + nu + vecRhomax1)); %ok
lambda = (Rp-1).*beta.*(sigma+vecRhomax1+mu)./(beta+(sigma+vecRhomax1+mu).*Rp);
S  = b./(mu*Rp); %ok
E  = lambda.*S./(sigma+mu+vecRhomax1);
IA =  (1-p)*sigma.*E./(nu+vecRhomax1+mu);
newdiag = vecRhomax1.*(S+E+IA);
newinf  = lambda.*S;
ratio = newdiag./newinf.*(vecC>=c1);

plot(vecC,log(ratio))
plot_paper1_procedure2(log(ratio), c1, c0, max(log(ratio)), vecC)

%% verification de la version 2.2 où P=(E+IA)/N
clear all
[beta,beta2,gamma0,gamma2,nu,nu2,sigma,mu,rho] = random_parameters(true, true);
p=0.2;beta=100*beta;
souslaracine = (gamma0^2*nu^2 - 2*gamma0^2*nu*sigma - 4*beta*gamma0^2*p*sigma + gamma0^2*sigma^2 + 4*beta*gamma0^2*sigma + 2*gamma0*mu*nu^2 - 4*gamma0*mu*nu*sigma - 8*beta*gamma0*mu*p*sigma + 2*gamma0*mu*sigma^2 + 8*beta*gamma0*mu*sigma + 2*gamma0*nu^3 - 4*gamma0*nu^2*sigma - 6*beta*gamma0*nu*p*sigma + 2*gamma0*nu*sigma^2 + 8*beta*gamma0*nu*sigma - 2*beta*gamma0*p*sigma^2 + mu^2*nu^2 - 2*mu^2*nu*sigma - 4*beta*mu^2*p*sigma + mu^2*sigma^2 + 4*beta*mu^2*sigma + 2*mu*nu^3 - 4*mu*nu^2*sigma - 6*beta*mu*nu*p*sigma + 2*mu*nu*sigma^2 + 8*beta*mu*nu*sigma - 2*beta*mu*p*sigma^2 + nu^4 - 2*nu^3*sigma - 2*beta*nu^2*p*sigma + nu^2*sigma^2 + 4*beta*nu^2*sigma - 2*beta*nu*p*sigma^2 + beta^2*p^2*sigma^2);
alpha = (beta*p*sigma + sqrt(souslaracine))/(2*(gamma0+mu+nu)) - (2*mu+nu+sigma)/2;
Rp = (beta*sigma*(gamma0*(1-p) + mu + nu + p*rho))./((mu + sigma+rho).*(gamma0 + mu + nu).*(mu + nu + rho));

[U,dU,P] = U_SEIISv2(p, beta, sigma, gamma0, mu, nu, rho, 0, alpha)

b=2;
Emat= -((mu + nu + rho)*(b*mu^4 + 2*b*mu^3*rho + b*nu^3*rho + b*mu^3*sigma + b*nu^3*sigma + b*gamma0^2*mu^2 + b*gamma0^2*rho^2 + 3*b*mu^2*nu^2 + b*mu^2*rho^2 + b*nu^2*rho^2 + 2*b*gamma0*mu^3 + b*mu*nu^3 + 3*b*mu^3*nu - b*beta*gamma0^2*sigma - b*beta*mu^2*sigma + 2*b*gamma0*mu*nu^2 + 4*b*gamma0*mu^2*nu + b*gamma0^2*mu*nu - b*beta*nu^2*sigma + 2*b*gamma0*mu*rho^2 + 4*b*gamma0*mu^2*rho + 2*b*gamma0^2*mu*rho + 2*b*gamma0*nu*rho^2 + 2*b*gamma0*nu^2*rho + b*gamma0^2*nu*rho + 2*b*gamma0*mu^2*sigma + b*gamma0^2*mu*sigma + 2*b*gamma0*nu^2*sigma + b*gamma0^2*nu*sigma + b*gamma0^2*rho*sigma + 2*b*mu*nu*rho^2 + 4*b*mu*nu^2*rho + 5*b*mu^2*nu*rho + 3*b*mu*nu^2*sigma + 3*b*mu^2*nu*sigma + b*mu^2*rho*sigma + b*nu^2*rho*sigma - 2*b*beta*gamma0*mu*sigma - 2*b*beta*gamma0*nu*sigma - 2*b*beta*mu*nu*sigma + 6*b*gamma0*mu*nu*rho + 4*b*gamma0*mu*nu*sigma + 2*b*gamma0*mu*rho*sigma + 2*b*gamma0*nu*rho*sigma + 2*b*mu*nu*rho*sigma + b*beta*gamma0^2*p*sigma + b*beta*gamma0*mu*p*sigma + b*beta*gamma0*nu*p*sigma - b*beta*gamma0*p*rho*sigma - b*beta*mu*p*rho*sigma - b*beta*nu*p*rho*sigma))/(sigma*(- beta*gamma0^2*mu^2*p + beta*gamma0^2*mu^2 - beta*gamma0^2*mu*nu*p + beta*gamma0^2*mu*nu + sigma*beta*gamma0^2*mu*p^2 - beta*gamma0^2*mu*p*rho - 2*sigma*beta*gamma0^2*mu*p + beta*gamma0^2*mu*rho + sigma*beta*gamma0^2*mu - beta*gamma0*mu^3*p + 2*beta*gamma0*mu^3 - 2*beta*gamma0*mu^2*nu*p + 4*beta*gamma0*mu^2*nu - 2*sigma*beta*gamma0*mu^2*p + 2*beta*gamma0*mu^2*rho + 2*sigma*beta*gamma0*mu^2 - beta*gamma0*mu*nu^2*p + 2*beta*gamma0*mu*nu^2 - 2*sigma*beta*gamma0*mu*nu*p + 2*beta*gamma0*mu*nu*rho + 2*sigma*beta*gamma0*mu*nu - 2*sigma*beta*gamma0*mu*p^2*rho + beta*gamma0*mu*p*rho^2 + 2*sigma*beta*gamma0*mu*p*rho + beta*mu^4 + 3*beta*mu^3*nu + beta*mu^3*p*rho + beta*mu^3*rho + sigma*beta*mu^3 + 3*beta*mu^2*nu^2 + 2*beta*mu^2*nu*p*rho + 2*beta*mu^2*nu*rho + 2*sigma*beta*mu^2*nu + beta*mu^2*p*rho^2 + 2*sigma*beta*mu^2*p*rho + beta*mu*nu^3 + beta*mu*nu^2*p*rho + beta*mu*nu^2*rho + sigma*beta*mu*nu^2 + beta*mu*nu*p*rho^2 + 2*sigma*beta*mu*nu*p*rho + sigma*beta*mu*p^2*rho^2))
IAmat = (b*gamma0^2*mu^2*p - b*gamma0^2*mu^2 + b*gamma0^2*mu*nu*p - b*gamma0^2*mu*nu + 2*b*gamma0^2*mu*p*rho + b*sigma*gamma0^2*mu*p - 2*b*gamma0^2*mu*rho - b*sigma*gamma0^2*mu + b*gamma0^2*nu*p*rho + b*sigma*gamma0^2*nu*p - b*gamma0^2*nu*rho - b*sigma*gamma0^2*nu + b*sigma*beta*gamma0^2*p^2 + b*gamma0^2*p*rho^2 + b*sigma*gamma0^2*p*rho - 2*b*sigma*beta*gamma0^2*p - b*gamma0^2*rho^2 - b*sigma*gamma0^2*rho + b*sigma*beta*gamma0^2 + 2*b*gamma0*mu^3*p - 2*b*gamma0*mu^3 + 4*b*gamma0*mu^2*nu*p - 4*b*gamma0*mu^2*nu + 4*b*gamma0*mu^2*p*rho + 2*b*sigma*gamma0*mu^2*p - 4*b*gamma0*mu^2*rho - 2*b*sigma*gamma0*mu^2 + 2*b*gamma0*mu*nu^2*p - 2*b*gamma0*mu*nu^2 + 6*b*gamma0*mu*nu*p*rho + 4*b*sigma*gamma0*mu*nu*p - 6*b*gamma0*mu*nu*rho - 4*b*sigma*gamma0*mu*nu + b*sigma*beta*gamma0*mu*p^2 + 2*b*gamma0*mu*p*rho^2 + 2*b*sigma*gamma0*mu*p*rho - 3*b*sigma*beta*gamma0*mu*p - 2*b*gamma0*mu*rho^2 - 2*b*sigma*gamma0*mu*rho + 2*b*sigma*beta*gamma0*mu + 2*b*gamma0*nu^2*p*rho + 2*b*sigma*gamma0*nu^2*p - 2*b*gamma0*nu^2*rho - 2*b*sigma*gamma0*nu^2 + b*sigma*beta*gamma0*nu*p^2 + 2*b*gamma0*nu*p*rho^2 + 2*b*sigma*gamma0*nu*p*rho - 3*b*sigma*beta*gamma0*nu*p - 2*b*gamma0*nu*rho^2 - 2*b*sigma*gamma0*nu*rho + 2*b*sigma*beta*gamma0*nu - b*sigma*beta*gamma0*p^2*rho + b*sigma*beta*gamma0*p*rho + b*mu^4*p - b*mu^4 + 3*b*mu^3*nu*p - 3*b*mu^3*nu + 2*b*mu^3*p*rho + b*sigma*mu^3*p - 2*b*mu^3*rho - b*sigma*mu^3 + 3*b*mu^2*nu^2*p - 3*b*mu^2*nu^2 + 5*b*mu^2*nu*p*rho + 3*b*sigma*mu^2*nu*p - 5*b*mu^2*nu*rho - 3*b*sigma*mu^2*nu + b*mu^2*p*rho^2 + b*sigma*mu^2*p*rho - b*sigma*beta*mu^2*p - b*mu^2*rho^2 - b*sigma*mu^2*rho + b*sigma*beta*mu^2 + b*mu*nu^3*p - b*mu*nu^3 + 4*b*mu*nu^2*p*rho + 3*b*sigma*mu*nu^2*p - 4*b*mu*nu^2*rho - 3*b*sigma*mu*nu^2 + 2*b*mu*nu*p*rho^2 + 2*b*sigma*mu*nu*p*rho - 2*b*sigma*beta*mu*nu*p - 2*b*mu*nu*rho^2 - 2*b*sigma*mu*nu*rho + 2*b*sigma*beta*mu*nu - b*sigma*beta*mu*p^2*rho + b*sigma*beta*mu*p*rho + b*nu^3*p*rho + b*sigma*nu^3*p - b*nu^3*rho - b*sigma*nu^3 + b*nu^2*p*rho^2 + b*sigma*nu^2*p*rho - b*sigma*beta*nu^2*p - b*nu^2*rho^2 - b*sigma*nu^2*rho + b*sigma*beta*nu^2 - b*sigma*beta*nu*p^2*rho + b*sigma*beta*nu*p*rho)/(- beta*gamma0^2*mu^2*p + beta*gamma0^2*mu^2 - beta*gamma0^2*mu*nu*p + beta*gamma0^2*mu*nu + sigma*beta*gamma0^2*mu*p^2 - beta*gamma0^2*mu*p*rho - 2*sigma*beta*gamma0^2*mu*p + beta*gamma0^2*mu*rho + sigma*beta*gamma0^2*mu - beta*gamma0*mu^3*p + 2*beta*gamma0*mu^3 - 2*beta*gamma0*mu^2*nu*p + 4*beta*gamma0*mu^2*nu - 2*sigma*beta*gamma0*mu^2*p + 2*beta*gamma0*mu^2*rho + 2*sigma*beta*gamma0*mu^2 - beta*gamma0*mu*nu^2*p + 2*beta*gamma0*mu*nu^2 - 2*sigma*beta*gamma0*mu*nu*p + 2*beta*gamma0*mu*nu*rho + 2*sigma*beta*gamma0*mu*nu - 2*sigma*beta*gamma0*mu*p^2*rho + beta*gamma0*mu*p*rho^2 + 2*sigma*beta*gamma0*mu*p*rho + beta*mu^4 + 3*beta*mu^3*nu + beta*mu^3*p*rho + beta*mu^3*rho + sigma*beta*mu^3 + 3*beta*mu^2*nu^2 + 2*beta*mu^2*nu*p*rho + 2*beta*mu^2*nu*rho + 2*sigma*beta*mu^2*nu + beta*mu^2*p*rho^2 + 2*sigma*beta*mu^2*p*rho + beta*mu*nu^3 + beta*mu*nu^2*p*rho + beta*mu*nu^2*rho + sigma*beta*mu*nu^2 + beta*mu*nu*p*rho^2 + 2*sigma*beta*mu*nu*p*rho + sigma*beta*mu*p^2*rho^2)
(Emat+IAmat)*mu/b
