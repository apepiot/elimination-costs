% Parameters
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters('true', 'true');
s1=1;
s2=1;

gamma1p = gamma1+s1.*rho;
gamma2p = gamma2+s2.*rho;
R1p = beta1/(gamma1p+mu);
R2p = beta2/(gamma2p+mu);

% Populations at the equilibrium
I12 = -(R2p.*b.*beta1.*beta2.^2 - b.*(R1p.^2.*R2p.^2.*beta1.^2.*beta2.^2.*mu.^2 + 2.*R1p.^2.*R2p.^2.*beta1.^2.*beta2.*gamma1.*mu.^2 + R1p.^2.*R2p.^2.*beta1.^2.*gamma1.^2.*mu.^2 - 2.*R1p.^2.*R2p.^2.*beta1.*beta2.^3.*mu.*rho - 2.*R1p.^2.*R2p.^2.*beta1.*beta2.^2.*gamma1.*mu.*rho + R1p.^2.*R2p.^2.*beta2.^4.*rho.^2 + 2.*R1p.^2.*R2p.*beta1.*beta2.^3.*mu.*rho + 2.*R1p.^2.*R2p.*beta1.*beta2.^2.*gamma1.*mu.*rho - 2.*R1p.^2.*R2p.*beta2.^4.*rho.^2 + R1p.^2.*beta2.^4.*rho.^2 + 2.*R1p.*R2p.^2.*beta1.^2.*beta2.^3.*mu + 4.*R1p.*R2p.^2.*beta1.^2.*beta2.^2.*gamma1.*mu - 2.*R1p.*R2p.^2.*beta1.^2.*beta2.^2.*mu.^2 + 2.*R1p.*R2p.^2.*beta1.^2.*beta2.^2.*mu.*rho + 2.*R1p.*R2p.^2.*beta1.^2.*beta2.*gamma1.^2.*mu - 4.*R1p.*R2p.^2.*beta1.^2.*beta2.*gamma1.*mu.^2 + 2.*R1p.*R2p.^2.*beta1.^2.*beta2.*gamma1.*mu.*rho - 2.*R1p.*R2p.^2.*beta1.^2.*beta2.*gamma2.*mu.*rho - 2.*R1p.*R2p.^2.*beta1.^2.*beta2.*mu.^2.*rho - 2.*R1p.*R2p.^2.*beta1.^2.*beta2.*mu.*rho.^2 - 2.*R1p.*R2p.^2.*beta1.^2.*gamma1.^2.*mu.^2 - 2.*R1p.*R2p.^2.*beta1.^2.*gamma1.*gamma2.*mu.*rho - 2.*R1p.*R2p.^2.*beta1.^2.*gamma1.*mu.^2.*rho - 2.*R1p.*R2p.^2.*beta1.^2.*gamma1.*mu.*rho.^2 - 2.*R1p.*R2p.^2.*beta1.*beta2.^4.*rho - 2.*R1p.*R2p.^2.*beta1.*beta2.^3.*gamma1.*rho + 2.*R1p.*R2p.^2.*beta1.*beta2.^3.*mu.*rho + 2.*R1p.*R2p.^2.*beta1.*beta2.^3.*rho.^2 + 2.*R1p.*R2p.^2.*beta1.*beta2.^2.*gamma1.*mu.*rho - 2.*R1p.*R2p.^2.*beta1.*beta2.^2.*gamma2.*rho.^2 - 2.*R1p.*R2p.^2.*beta1.*beta2.^2.*mu.*rho.^2 - 2.*R1p.*R2p.^2.*beta1.*beta2.^2.*rho.^3 + 2.*R1p.*R2p.*beta1.*beta2.^4.*rho + 2.*R1p.*R2p.*beta1.*beta2.^3.*gamma1.*rho - 2.*R1p.*R2p.*beta1.*beta2.^3.*mu.*rho - 2.*R1p.*R2p.*beta1.*beta2.^3.*rho.^2 - 2.*R1p.*R2p.*beta1.*beta2.^2.*gamma1.*mu.*rho + 2.*R1p.*R2p.*beta1.*beta2.^2.*gamma2.*rho.^2 + 2.*R1p.*R2p.*beta1.*beta2.^2.*mu.*rho.^2 + 2.*R1p.*R2p.*beta1.*beta2.^2.*rho.^3 + R2p.^2.*beta1.^2.*beta2.^4 + 2.*R2p.^2.*beta1.^2.*beta2.^3.*gamma1 - 2.*R2p.^2.*beta1.^2.*beta2.^3.*mu - 2.*R2p.^2.*beta1.^2.*beta2.^3.*rho + R2p.^2.*beta1.^2.*beta2.^2.*gamma1.^2 - 4.*R2p.^2.*beta1.^2.*beta2.^2.*gamma1.*mu - 2.*R2p.^2.*beta1.^2.*beta2.^2.*gamma1.*rho + 2.*R2p.^2.*beta1.^2.*beta2.^2.*gamma2.*rho + R2p.^2.*beta1.^2.*beta2.^2.*mu.^2 + 3.*R2p.^2.*beta1.^2.*beta2.^2.*rho.^2 - 2.*R2p.^2.*beta1.^2.*beta2.*gamma1.^2.*mu + 2.*R2p.^2.*beta1.^2.*beta2.*gamma1.*gamma2.*rho + 2.*R2p.^2.*beta1.^2.*beta2.*gamma1.*mu.^2 + 2.*R2p.^2.*beta1.^2.*beta2.*gamma1.*rho.^2 + 2.*R2p.^2.*beta1.^2.*beta2.*gamma2.*mu.*rho - 2.*R2p.^2.*beta1.^2.*beta2.*gamma2.*rho.^2 + 2.*R2p.^2.*beta1.^2.*beta2.*mu.^2.*rho - 2.*R2p.^2.*beta1.^2.*beta2.*rho.^3 + R2p.^2.*beta1.^2.*gamma1.^2.*mu.^2 + 2.*R2p.^2.*beta1.^2.*gamma1.*gamma2.*mu.*rho + 2.*R2p.^2.*beta1.^2.*gamma1.*mu.^2.*rho + 2.*R2p.^2.*beta1.^2.*gamma1.*mu.*rho.^2 + R2p.^2.*beta1.^2.*gamma2.^2.*rho.^2 + 2.*R2p.^2.*beta1.^2.*gamma2.*mu.*rho.^2 + 2.*R2p.^2.*beta1.^2.*gamma2.*rho.^3 + R2p.^2.*beta1.^2.*mu.^2.*rho.^2 + 2.*R2p.^2.*beta1.^2.*mu.*rho.^3 + R2p.^2.*beta1.^2.*rho.^4).^(1./2) + R1p.*b.*beta2.^2.*rho + R2p.*b.*beta1.*rho.^2 + R2p.*b.*beta1.*beta2.*gamma1 - R2p.*b.*beta1.*beta2.*mu - R2p.*b.*beta1.*beta2.*rho - R2p.*b.*beta1.*gamma1.*mu + R2p.*b.*beta1.*gamma2.*rho - 2.*R1p.*b.*beta2.*mu.*rho + R2p.*b.*beta1.*mu.*rho - R1p.*R2p.*b.*beta2.^2.*rho + 2.*R1p.^2.*b.*beta2.*mu.*rho - 2.*R1p.^2.*R2p.*b.*beta2.*mu.*rho + R1p.*R2p.*b.*beta1.*beta2.*mu + R1p.*R2p.*b.*beta1.*gamma1.*mu + 2.*R1p.*R2p.*b.*beta2.*mu.*rho)./(2.*(R1p.^2.*beta2.*mu.*rho.^2 - R1p.^2.*R2p.*beta2.*mu.*rho.^2 + R1p.*R2p.*beta1.*beta2.*mu.*rho + R1p.*R2p.*beta1.*gamma1.*mu.*rho));
I1 = (b/beta1.*(R1p-1)+R1p./beta1.*s1.*s2.*rho.*I12) - I12;
S = (b/mu.*(mu+gamma2p./R1p)+rho.*I12)./(mu.*(R1p-1)+beta2+rho.*mu/b.*R1p.*I12);
IR2 = b/mu.*(1-1./R2p-1./R1p) + S - I12;
I2 = b./(mu.*R1p) - S;
R1 = b/mu - S - I1 - I2 - I12 - IR2;

%Numerically
tspan = 0:1:100000;
Y0 = [15,10,2,3,4,0];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_V4(t,Y,b,beta1,beta2,s1,s2,rho,gamma1,gamma2,mu,'frequency'),tspan,Y0, options);
T = Ys(end,:)
T - [S,I1,I2,I12,IR2,R1]


%Conditions a verifier 

%I1+I12 = b/beta1*(R1p-1)+R1p/beta1*s1*s2*I12
cond1 = I1 + I12 - (b/beta1.*(R1p-1)+R1p./beta1.*s1.*s2.*rho.*I12); %0

%S+I2 = b/(mu*R1p)
cond2 = S+I2 - b./(mu.*R1p); %0

%S+I1+R1 = b/(mu*R2p)
cond3 = S+I1+R1 - b./(mu.*R2p);%0

%I2 + I12 + IR2 = b/mu*(1-1/R2p);
cond4 = I2 + I12 + IR2 - b/mu.*(1-1./R2p); %0


%lambda1 = beta1.*mu/b.*(b/beta1.*(R1p-1) + R1p/beta1.*s1.*s2.*rho.*I12); 
lambda1 = beta1*mu/b.*(I2+I12);
lambda2 = beta2.*(1-1/R2p); %beta2.*(IR2 + I2 + I12)/N;
gamma1p = gamma1+s1.*rho;
gamma2p = gamma2+s2.*rho;
gamma1t = gamma1p-s1.*s2.*rho;
gamma2t = gamma2p-s1.*s2.*rho;
gamma12t = s1.*s2.*rho;
gamma12 = gamma1p+gamma2p-s1.*s2.*rho;

eqn1 = b - (lambda1 + lambda2).*S + gamma2p.*I2 + gamma12t.*I12- mu.*S==0; %S
eqn2 = lambda1.*S - lambda2.*I1 + gamma2t.*I12 - (gamma1p + mu).*I1 == 0; %I1
eqn3 = lambda2.*S - lambda1.*I2 - (gamma2p + mu).*I2 == 0; %I2
eqn4 = lambda2.*I1 + lambda1.*I2 - (gamma12 + mu).*I12 ==0;%I12
eqn5 = gamma1t.*I12 + lambda2.*R1  - (gamma2p + mu).*IR2==0; %IR2
eqn6 = gamma1p.*I1 - lambda2.*R1 + gamma2p.*IR2 - mu.*R1==0; % R1