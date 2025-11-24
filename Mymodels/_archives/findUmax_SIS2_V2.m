%clear all
close all

% 
R1 = beta1./(gamma1+mu);
R2 = beta2./(gamma2+mu);
alpha1 = beta1./s1*(1-1./R1);
alpha2 = beta2./s2*(1-1./R2);

% rho's sampling

maxRho2 = (beta1+beta2-2*mu-gamma1/R2-gamma2/R1)/2/(s2/R1+s2*gamma1/beta2);
maxRho1 = (beta1+beta2-2*mu-gamma1/R2-gamma2/R1)/2/(s1/R2+s1*gamma2/beta1);


rho1 = 0:(alpha1/1000):1.1*alpha1;
%rho2 = [0:(maxRho2/100):0.98*maxRho2,0.98*maxRho2:(maxRho2/1000):1.1*maxRho2];
rho2 = 0:(alpha2/1000):1.1*alpha2;
[RHO1,RHO2] = meshgrid(rho1,rho2);

% gammap's & Rp's
gamma1p = gamma1+s1.*RHO1;
gamma2p = gamma2+s2.*RHO2;
R1p = beta1./(gamma1p+mu);
R2p = beta2./(gamma2p+mu);

% prevalences
P1 = 1-1./R1p;
P2 = 1-1./R2p;
P12 = (1 - (gamma1p./R2p + gamma2p./R1p + mu)./(beta1+beta2-mu)).*(RHO1<alpha1 & RHO2<alpha2); %COINFECTION EQU
P = P12.*(RHO1<alpha1 & RHO2<alpha2) + P1.*(RHO1<alpha1 & RHO2>=alpha2) + P2.*(RHO1>=alpha1 & RHO2<alpha2);

% utilities
U1 = RHO1.*P1;
U2 = RHO2.*P2;
U12 = max((RHO1 + RHO2).*P12,0);
U = U12.*(RHO1<=alpha1 & RHO2<=alpha2) + U1.*(RHO1<=alpha1 & RHO2>alpha2) + U2.*(RHO1>alpha1 & RHO2<alpha2);

%
[maxU, c1] = max(U(:));
[imax,jmax] = find(ismember(U, max(U(:))));

U(imax,jmax);
rho1max = rho1(jmax); %???
rho2max = rho2(imax);


% plot
if 0 
figure(1)
surf(RHO1,RHO2,U)
colormap(jet)
title([{'Utility of the SISxSIS model function of \rho_1 and \rho_2'},...
        {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
        num2str(gamma2), ' s1=', num2str(s1),' s2=', num2str(s2), '\mu=', num2str(mu), ...
        ' R_0^1=' num2str(round(R1,2)), ' R_0^2=' num2str(round(R2,2)),...
        ' (\rho_1,\rho_2)_{max}=', '(',num2str(round(rho1max,2)),',',num2str(round(rho2max,2)),')',...
        ' \rho_1''=',num2str(round(alpha1,2)),' \rho_2''=',num2str(round(alpha2,2))]}])
xlabel("Voluntary-testing rate \rho_1","fontweight","bold")
ylabel("Voluntary-testing rate \rho_2","fontweight","bold")
lim = caxis;

figure(2)
surf(RHO1,RHO2,U12)
caxis(lim);
colormap(jet)
title([{'Utility U_{12} of the SISxSIS model function of \rho_1 and \rho_2'},...
        {['\beta_1=',num2str(beta1), ' \beta_2=',num2str(beta2),' \gamma_1=',num2str(gamma1),' \gamma_2=',...
        num2str(gamma2), ' s1=', num2str(s1),' s2=', num2str(s2), '\mu=', num2str(mu), ...
        ' R_0^1=' num2str(round(R1,2)), ' R_0^2=' num2str(round(R2,2)),...
        ' (\rho_1,\rho_2)_{max}=', '(',num2str(round(rho1max,2)),',',num2str(round(rho2max,2)),')',...
        ' \rho_1''=',num2str(alpha1),' \rho_2''=',num2str(alpha2)]}])
xlabel("Voluntary-testing rate \rho_1","fontweight","bold")
ylabel("Voluntary-testing rate \rho_2","fontweight","bold")
zlabel('U_{12}=(\rho_1+\rho_2)\Pi_{12}')
end
%% Max of U12 (numerically)
[maxU12, c1] = max(U12(:));
[imax,jmax] = find(ismember(U12, max(U12(:))));

U12(imax,jmax);
rho1max = rho1(jmax);%???
rho2max = rho2(imax);

%max en rho2 quand rho1=rho1' (a verifier)
rho21 = (beta1+beta2-2*mu - (gamma1+s1*alpha1)/R2-gamma2)/(2*s2*((gamma1+s1*alpha1)/beta2+1)) - alpha1/2;
Ucase1 = (alpha1+alpha2)^2/4/beta2; %case simple

%max en rho1 quand rho2=rho2' (a verifier)
%rho12 = (beta1+beta2-2*mu - (gamma2+s2*alpha2)/R1-gamma1)/(2*s1*((gamma2+s2*alpha2)/beta1+1)) - alpha2/2;
%Ucase2 = (alpha1+alpha2).^2/(4*beta1); %case simple
rho12 = alpha2;
Ucase2 = alpha2*(1-1/R1);

%max en rho2 quand rho1=0
rho20 = (beta1+beta2-2*mu-gamma1/R2-gamma2/R1)/2/(s2/R1+s2*gamma1/beta2);
Ucase3 = (((beta1*beta2 - gamma1.*gamma2).^2)./(4*beta1*beta2*gamma1))*(alpha2>=rho20); %case simple

%max en rho1 quand rho2=0
rho10 = (beta1+beta2-2*mu-gamma1/R2-gamma2/R1)/2/(s1/R2+s1*gamma2/beta1);
Ucase4 = ((beta1*beta2 - gamma1.*gamma2).^2./(4*beta1*beta2*gamma2))*(alpha1>=rho10); %case simple



%% Conditions pour le max de U12

% if machin est le max verifier que...
% if on a cette conditions, verifier que...
% a la maniere de la verification des conditions 
%x1 = beta2*(-alpha1/(sqrt(beta2*gamma1))+1);
%x2 = beta2*(alpha1/(sqrt(beta2*gamma1))+1);
%x3 = beta1*(-alpha2/(sqrt(beta1*gamma2))+1);
%x4 = beta1*(alpha2/(sqrt(beta1*gamma2))+1);
x5 = beta1*((gamma1-beta2)/sqrt(gamma1*beta1)+1);
x6 = beta1*((-gamma1+beta2)/sqrt(gamma1*beta1)+1);
%x7 = beta2*((gamma2-beta1)/sqrt(beta2*gamma2)+1);
%x8 = beta2*((-gamma2+beta1)/sqrt(beta2*gamma2)+1);
x11 = beta1 + alpha2*(1-2*beta2/beta1*(1-sqrt(1-beta1/beta2)));

%conditions have to be fulfilled if u1,u2 or u3 is the max
cond11 = (alpha1<rho10 & alpha2>rho20 & (gamma2<=min(x5,x6) | gamma2>max(x5,x6)));
cond12 = (alpha2<rho20 & (beta2<beta1 | (gamma1>x11 & beta2>=beta1)));
cond1 = (cond11 | cond12);
cond2 = (alpha2<rho20 & (gamma1<x11 & beta2>=beta1));
%cond22 = (alpha1>=rho10 & alpha2<rho20);
%cond2 = (cond21 | cond22);
cond3  = (alpha1<rho10 & alpha2>rho20 & gamma2>min(x5,x6) & gamma2<max(x5,x6));
cond4 = 0;




