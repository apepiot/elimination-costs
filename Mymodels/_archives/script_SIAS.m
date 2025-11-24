%script for the SIAS model
[beta,beta2,gamma0,ksi0,s,s2,b,mu] = random_parameters('true', 'true');
omega=0.4;
rho=0;

% Parametres du systeme d'ODE 
tspan = 0:1:5000;
Y0 = [2; 0.1; 0.1];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_SIAS(t,Y,b,beta,gamma,omega,ksi0,s,rho,mu,'frequency'),tspan,Y0, options);
sum(Ys(end,:))

R = beta/(gamma+mu);
b/mu*(1-1/R)
sum(Ys(2:3,end))
N=b/mu;
S=Ys(end,1);I=Ys(end,2);A=Ys(end,3);
I+A
(gamma*I+(ksi0+rho)*A)/(beta*S/N-mu)
(beta*S*(I+A)/N+mu*S-b)/(beta*S/N-mu)

ksi=ksi0+s*rho;
R = beta*(gamma*(1-omega)+ksi*omega+mu)/((gamma+mu)*(ksi+mu));
b/mu/R
Ys(end,1)
(I+A)
sum(Ys(end,2:3))
b/mu*(1-1/R)

b/mu*(1-1/R)*(mu-(1-omega)*beta/R+ksi)/(ksi+mu)
b/mu*(1-1/R)*(1-1/R*(1-omega)*beta/(ksi+mu))
b/mu*(1-1/R)*(1-1/R*omega*beta/(gamma+mu))
Ys(end,3)
R = beta*(gamma0*(1-omega)+ksi0*omega+mu)/((gamma0+mu)*(ksi0+mu))

I

%%
clear all;
syms beta gamma ksi s rho mu S I A b
N=b/mu;
lambda = beta*(I+A)/N;

eqn1= 0== b - lambda*S + gamma*I + (ksi)*A - mu*S; %S
eqn2= omega*lambda*S - (gamma + mu)*I; %I
eqn3 = (1-omega)*lambda*S - (ksi + mu)*A; %A

sol=solve([eqn1,eqn2,eqn3],[S,I,A])

%% utility
clear all;
syms beta gamma0 ksi rho mu S I A b ksi0 omega s
ksi = ksi0+s*rho;
gamma = gamma0+s*rho;
R = beta*(gamma*(1-omega)+ksi*omega+mu)/((gamma+mu)*(ksi+mu));
P = 1-1/R;

solve(R==1,rho)
-(gamma0/2 - beta/2 + ksi0/2 + mu - (2*beta*gamma0 - 2*beta*ksi0 - 2*gamma0*ksi0 + beta^2 + gamma0^2 + ksi0^2 - 4*beta*gamma0*omega + 4*beta*ksi0*omega)^(1/2)/2)/s
%-(gamma0/2 - beta/2 + ksi0/2 + mu + (2*beta*gamma0 - 2*beta*ksi0 - 2*gamma0*ksi0 + beta^2 + gamma0^2 + ksi0^2 - 4*beta*gamma0*omega + 4*beta*ksi0*omega)^(1/2)/2)/s
%% utility function
clear all;
%[beta,beta2,gamma0,ksi0,s,s2,b,mu] = random_parameters('true', 'true');
%s=1;
%omega=0.4;
%rho=0:0.01:20;
ksi0=0; gamma0=12/1;
s=1;rho=0;
omega=0.5; mu=1/35;
%R0 = beta.*(gamma0.*(1-omega)+ksi0*omega+mu)./((gamma0+mu).*(ksi0+mu))
ro = 3;
beta = (ro*(gamma0 + mu)*(ksi0 + mu))/(mu + ksi0*omega - gamma0*(omega - 1));
b=1;

%%
ksi = ksi0+s*rho;
gamma = gamma0+s*rho;

R = beta.*(gamma.*(1-omega)+ksi*omega+mu)./((gamma+mu).*(ksi+mu));
P = 1-1./R;

% U=max(rho.*P,0);
% figure;
% plot(rho,U)
%%
clear rho
alpha = -(gamma0/2 - beta/2 + ksi0/2 + mu - (2*beta*gamma0 - 2*beta*ksi0 - 2*gamma0*ksi0 + beta^2 + gamma0^2 + ksi0^2 - 4*beta*gamma0*omega + 4*beta*ksi0*omega)^(1/2)/2)/s
i=0;
[U1,P1,c1] = U_SIAS(alpha,beta,gamma0,ksi0,omega,s,b,mu,0);
[U0,c0]= U_SIAS(0,beta,gamma0,ksi0,omega,s,b,mu,0); %c0 = P(0)
rhohat=[]; vecC=(c1-(c0-c1)):(c0-c1)/100:c0+(c0-c1);
options = optimset('Display','off'); %options for minsearch
for c=vecC
    i=i+1;
    fun = @(rho) -U_SIAS(rho,beta,gamma0,ksi0,omega,s,b,mu,c); 
    rhohat(i) = min(max(fminsearch(fun,0,options),0),alpha);
end
figure
plot(vecC,rhohat)
hold on;

plot_paper1_procedure2(rhohat, c1, c0, alpha, vecC)
title(['SIAS : ','$\gamma(0)=$',num2str(gamma0),' $\xi(0)=$',num2str(ksi0),...
    ' $ \omega=$',num2str(omega),' $ \mathtt R(0)=$',num2str(ro)], 'Interpreter', 'latex')
