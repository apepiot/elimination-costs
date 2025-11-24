%% SIRxSIS - strategy 1 - U = rho1(P1+c1) + rho2(P2+c2)
clear all;close all;
%% Parametres
% Populations initiales
S0      = 50;
I10     = 5;
I20     = 3;
I30     = 2;    %I_12 (coinfection)
I40     = 0;
R10     = 0.;

% Definition des parametres 
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true)
rho1=rand(1,1)*rho;rho2=rand(1,1)*rho;
%% Equilibrium (numerically)
% Parametres du systeme d'ODE 
tspan = 0:1:10000;
Y0 = [S0; I10; I20; I30; I40; R10];
options = odeset('RelTol',1e-5,'Stats','on');%,'OutputFcn',@odeplot);
[ts,Ys] = ode45(@(t,Y) ODE_V2(t,Y,b,beta1,beta2,s1,s2,gamma1,gamma2,mu,rho1,rho2,'frequency'),tspan,Y0, options);
T=Ys(end,:)

%% Equilibrium (theoretically)
gamma1p = gamma1+s1*rho1;gamma2p = gamma2+s2*rho2;
R1p = beta1/(mu+gamma1p);R2p = beta2/(gamma2p+mu);
lambda1=mu*(R1p-1);lambda2=beta2*(1-1/R2p);

% Only disease one : ES1
SES1 = b/(mu*R1p); I1ES1 = b/beta1*(R1p-1); RES1 = gamma1p/mu*b/beta1*(R1p-1);

%Only disease two : ES2
SES2 = b/(mu*R2p);I2ES2=b/mu*(1-1/R2p);

%Coinfection equilibrium : ES12
SES12   = b/mu*(mu*(R1p-1) + gamma2p+mu)/(mu*(R1p-1)+beta2)/R1p;
I1ES12  = b*(R1p-1)/(beta2+gamma1p)*(gamma2p/beta1+mu/b*SES12);
b/mu*lambda1/R1p*((lambda1+gamma2p+mu)/(lambda1+beta2)+ gamma2p/(gamma1p+mu))/(lambda2+gamma1p+mu+gamma2p)
I2ES12  = b/(mu*R1p)-SES12;
I12ES12 = b/beta1*(R1p-1)-I1ES12;
b/mu*lambda1/R1p*((lambda2+gamma1p+mu)/(gamma1p+mu) - (lambda1+mu+gamma2p)/(lambda1+beta2))/(beta2+gamma1p)
IR2ES12 = b/mu*(1-1/R2p)-I2ES12-I12ES12;
R1ES12  = b/(mu*R2p)-SES12 -I1ES12;
N = b/mu;

%% Endemic prevalence
R10 = beta1/(gamma1+mu);R20 = beta2/(gamma2+mu);
alpha1=beta1/s1*(1-1/R10); alpha2=beta2/s2*(1-1/R20);
P1 = mu/beta1*(R1p-1)*(rho1<alpha1) + 0*(rho1>=alpha1);
P2 = (1-1/R2p)*(rho2<alpha2) + 0*(rho2>=alpha2);

%% Utility function (theoretically)
c1=0;c2=0;
U = rho1.*(P1+c1) + rho2.*(P2+c2);

%% argmaxU
rho1max = 0*(c1<=mu/beta1*(1-R10))+...
    (beta1*(sqrt(R10*mu/(mu-beta1*c1))-1)/R10/s1)*(c1>mu/beta1*(1-R10) & c1<mu/beta1*(1-1/R10))+...
    alpha1*(c1>=mu/beta1*(1-1/R10));
rho2max = 0*(c2<=1/R20-1)+...
    (beta2/(2*s2)*(1-1/R20+c2))*(c2>1/R20-1 & c2<1-1/R20)+...
    alpha2*(c2>=1-1/R20);

%% 2D Utility function
clear all;close all;
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true)
c1=-0.5;c2=-0.1;

R10 = beta1/(gamma1+mu);R20 = beta2/(gamma2+mu);
alpha1=beta1/s1*(1-1/R10); alpha2=beta2/s2*(1-1/R20);

rho1 = 0:alpha1/50:alpha1*2;
rho2 = 0:alpha2/50:alpha2*2;
[RHO1,RHO2] = meshgrid(rho1,rho2);

gamma1p = gamma1+s1*RHO1;gamma2p = gamma2+s2*RHO2;
R1p = beta1./(mu+gamma1p);R2p = beta2./(gamma2p+mu);
P1 = mu/beta1.*(R1p-1).*(RHO1<alpha1) + 0.*(RHO1>=alpha1);
P2 = (1-1./R2p)*(RHO2<alpha2) + 0.*(RHO2>=alpha2);
U = max(RHO1.*(P1+c1) + RHO2.*(P2+c2),0);

%plot
surf(RHO1,RHO2,U)
xlabel('$\rho_1$','Interpreter','latex')
ylabel('$\rho_2$','Interpreter','latex')
zlabel('$U(\rho_1,\rho_2)$','Interpreter','latex')
title([{'SIR $\times$ SIS : $U(\rho_1,\rho_2)=\rho_1(\Pi_1(\rho_1)+c_1)+\rho_2(\Pi_2(\rho_2)+c_2)$'},...
      {['$\beta_1=$',num2str(beta1),' $\beta_2=$',num2str(beta2),...
      ' $\gamma_1(0)=$',num2str(round(gamma1,2)),' $\gamma_2(0)=$',num2str(round(gamma2,2)),...
      ' $\mu=$', num2str(round(mu,2)),...
      ' $\mathtt{R}_1(0)$=', num2str(round(R10,2)),' $\mathtt{R}_2(0)$=', num2str(round(R20,2)),...
      ' $s_1$=', num2str(s1),' $s_2$=', num2str(s2)]},...
    {[' $\rho_1\prime$=', num2str(round(alpha1,2)),' $\rho_2\prime$=', num2str(round(alpha2,2)),...
    ', $c_1=$', num2str(c1), ', $c_2=$', num2str(c2),...
    ', $c_{11}=$', num2str(round(mu/beta1*(1-R10),2)),...
    ', $c_{12}=$', num2str(round(mu/beta1*(1-1/R10),2)),...
     ', $\frac{\mu}{\beta_1}=$', num2str(round(mu/beta1,2)),...
    ', $c_{21}=$', num2str(round(1/R20-1,2)), ', $c_{22}=$', num2str(round(1-1/R20,2)) ]}],...
    'Interpreter','latex')

rho1max = 0*(c1<=mu/beta1*(1-R10))+...
    (beta1*(sqrt(R10*mu/(mu-beta1*c1))-1)/R10/s1)*(c1>mu/beta1*(1-R10) & c1<mu/beta1*(1-1/R10))+...
    alpha1*(c1>=mu/beta1*(1-1/R10));
rho2max = 0*(c2<=1/R20-1)+...
    (beta2/(2*s2)*(1-1/R20+c2))*(c2>1/R20-1 & c2<1-1/R20)+...
    alpha2*(c2>=1-1/R20);
text(rho1max,-alpha2/100,0,'$\hat \rho_1$','Interpreter','latex')
text(-alpha1/100,rho2max,0,'$\hat \rho_2$','Interpreter','latex')

%% 15/04/21 - R(hatrho) en fonction de c - plot
% on definit U(rho1,rho2) = rho1*(Pi1+c1) + rho2*(Pi2+c2)
clear all;close all; Tmax=10000;t=1;
vecCond = [];
cond=1;
while(t<Tmax & cond)
    t=t+1;
    [beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);

    c = -2:0.01:2;
    c1=c;c2=c;

    % disease 1 : SIR
    R10 = beta1/(gamma1+mu); alpha1 = beta1/s1*(1-1/R10);
    rhohat1 = beta1./(R10*s1).*(sqrt(R10*mu./(mu-beta1*c1))-1).*(c1>mu/beta1.*(1-R10) & c1<mu/beta1.*(1-1./R10))+...
        alpha1.*(c1>=mu/beta1.*(1-1./R10)) + 0.*(c1<=mu/beta1.*(1-R10));
    R1hat = beta1./(gamma1+s1.*rhohat1+mu);

    %disease 2 : SIS
    R20 = beta2/(gamma2+mu); alpha2=beta2/s2*(1-1/R20);
    rhohat2 = 0.*(c2<=1-1./R20) + alpha2.*(c2>=1./R20-1) + beta2/s2.*(1-1./R20+c2).*(c2>1-1./R20 & c2<1./R20-1); %un truc qui va pas evec rho2hat
    R2hat = beta2./(gamma2+s2.*rhohat2 + mu);

    cond1 = alpha1<alpha2;
    cond2 = 1/R20-1 < mu/beta1*(1-R10);
    cond3 = R10<R20;
    
    if(cond3)
        figure(1) %rhohats
        plot(c,rhohat1,c,rhohat2)
        legend('$\hat{\rho}_1$','$\hat{\rho}_2$','Interpreter','latex')

        figure(2) %Rrhohats
        plot(c,R1hat,c,R2hat)
        legend('$\mathtt R_1(\hat{\rho}_1)$','$\mathtt R_2(\hat{\rho}_2)$','Interpreter','latex')
        cond=0;
    end
    vecCond = [vecCond;cond1,cond2,cond3];
end

%% on compare 1/R2-1 et mu/beta1*(1-R1)
clear all; close all;
syms beta1 beta2 gamma1 gamma2 mu
R10 = beta1/(gamma1+mu);R20 = beta2/(gamma2+mu);
res = 1/R20-1 - mu/beta1*(1-R10);
solve(res==0,beta2)


%% Plot des differentes surfaces 080721 3D
close all; clear all;
c11=-0.3;
c12=0.2;
c21=-0.25;
c22=0.25;

c1 = -0.8:0.005:.7; c2=-0.75:0.005:0.75;
[C1,C2] = meshgrid(c1,c2);
Z = 0*C1;
J=[1,1,0];B=[0,0,1];
% kB = (C1-c11)/(c12-c11).*(C1<c12 & C1>=c11) + (C1>=c12); %linear shading
% kJ = (C2-c21)/(c22-c21).*(C2<c22 & C2>=c21) + (C2>=c22); %linear shading

d1 = 1.5; %=mu/beta;
a1 = (d1-c12)/(c12-c11);
b1 = d1-c11;
kJ = a1.*(b1./(d1-C1)-1).*(C1<c12 & C1>=c11) + (C1>=c12); %jaune pour C1
kB = (C2-c21)/(c22-c21).*(C2<c22 & C2>=c21) + (C2>=c22);
 
col(:,:,1) = 1-kB; % red
col(:,:,2) = 1-kB; % green
col(:,:,3) = 1-kJ; % blue
% s=surface(C1,C2,Z,col);
% s.EdgeColor = 'none';
% xlim([-0.8,0.7])
% ylim([-0.75,0.75])

figure
h1 = axes ;
image(-0.8:0.01:.7,-0.75:0.01:0.75,col)
set(h1,'YDir','normal'); %to make y axe sorted
hold on
%hatched region
[a] = fill([c12,c12,-0.8,-0.8,0.7,0.7],[-0.75,c22,c22,0.75,0.75,-0.75],'w') %color has no purpose here
hPatch1 = findobj(a, 'Type', 'patch');
hh1 = hatchfill(hPatch1, 'single', -45, 10); 
set(a, 'LineStyle', 'none') %remove outline
set(hh1, 'Color', 0.4*[1,1,1]) %set color of the hatched region

%Remove axes tic
set(gca,'YTick',[]); %which will get rid of all the markings for the y axis
set(gca,'XTick',[]); %which will get rid of all the markings for the x axis


colorbar
cb1 = 1-([0:0.01:1]'*[1,1,0]) ;
colormap(cb1);
cb11 = colorbar;
cb11.Location = 'westoutside';
cb11.Label.String = '$\hat\rho_2$';
cb11.Label.Interpreter = 'latex';
cb11.XTickLabel = [{'0'},{''},{''},{''},{''},{'$\rho_2\prime$'}];
set(cb11,'TickLabelInterpreter','latex')

%newcolorbar
hold on

cb22 = newcolorbar('southoutside');
plot(0,0)
cb2 = 1-([0:0.01:1]'*[0,0,1]) ;
colormap(gca,cb2);
cb22.Label.String = '$\hat\rho_1$';
cb22.Label.Interpreter = 'latex';
cb22.XTickLabel = [{'0'},{''},{''},{''},{''},{'$\rho_1\prime$'}];
set(cb22,'TickLabelInterpreter','latex')


hold on


plot([-0.8,0.7],[-0.25,-0.25], 'Color',0.85*[1,1,1],'LineStyle','--')%c21
plot([-0.8,0.7],[0.25,0.25], 'Color',0.85*[1,1,1],'LineStyle','--')%c22
plot([-0.3,-0.3],[-0.75,0.75], 'Color',0.85*[1,1,1],'LineStyle','--')%c11
plot([0.2,0.2],[-0.75,0.75], 'Color',0.85*[1,1,1],'LineStyle','--')%c12

text(-0.55, -0.5,'$(0,0)$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(0.45, 0.5,'$(\rho_1\prime,\rho_2\prime)$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','Color','w')
text(-0.55, 0.5,'$(0,\rho_2\prime)$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','Color','w')
text(0.45, -0.5,'$(\rho_1\prime,0)$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')

text(-0.3, -0.75,'$c_{11}$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(0.2, -0.75,'$c_{12}$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(-0.8, -0.25,'$c_{21}$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(-0.8, 0.25,'$c_{22}$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')

xlim([-0.8,0.7])
ylim([-0.75,0.75])
axis off


% rhochap in fct of c (same parameters)
% kJ = a1.*(b1./(d1-c1)-1).*(c1<c12 & c1>=c11) + (c1>=c12); %jaune pour C1
% kB = (c2-c21)/(c22-c21).*(c2<c22 & c2>=c21) + (c2>=c22);
% figure
% plot(c1,kJ, 'Color', [1,1,0])
% hold on
% plot(c2,kB, 'Color', [0,0,1])

%% 3D graph with c, rho1, rho2
close all ; clear all;
[beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(true, true);
R10=beta1/(gamma1+mu); R20=beta2/(gamma2+mu); alpha1=beta1/s1*(1-1/R10);alpha2=beta2/s2*(1-1/R20);

vecC = -2:0.01:2;
C = meshgrid(vecC);
c11 = mu/beta1*(1-R10); c12 = mu/beta1*(1-1/R10);
c21 = 1/R20-1; c22=-c21;
% rhohat1 = beta1./(R10*s1).*(sqrt(R10*mu./(mu-beta1*vecC))-1).*(vecC>c11 & vecC<c12)+...
%     alpha1.*(vecC>=c12) + 0.*(vecC<=c11);
% rhohat2 = 0.*(vecC<=c21) + alpha2.*(vecC>=c22) + beta2/s2.*(1-1./R20+vecC).*(vecC>c21 & vecC<c22);
RHO1 = beta1./(R10*s1).*(sqrt(R10*mu./(mu-beta1*C))-1).*(C>c11 & C<c12)+...
    alpha1.*(C>=c12) + 0.*(C<=c11);
RHO2 = 0.*(C<=c21) + alpha2.*(C>=c22) + beta2/(2*s2).*(1-1./R20+C).*(C>c21 & C<c22);

%color %ne fonctionne pas
d1 = mu/beta1;
a1 = (d1-c12)/(c12-c11);
b1 = d1-c11;
kJ = RHO1/alpha1; %jaune pour C1
kB = RHO2/alpha2;

col(:,:,1) = 1-kB; % red
col(:,:,2) = 1-kB; % green
col(:,:,3) = 1-kJ; % blue

%plot
surf(C,RHO1,RHO2)
ylabel('$\hat \rho_1$','Interpreter','latex')
zlabel('$\hat \rho_2$','Interpreter','latex')
xlabel('$c$','Interpreter','latex')



