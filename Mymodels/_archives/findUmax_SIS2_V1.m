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


rho1 = 0:(alpha1/500):1.1*alpha1;
%rho2 = [0:(maxRho2/100):0.98*maxRho2,0.98*maxRho2:(maxRho2/1000):1.1*maxRho2];
rho2 = 0:(alpha2/500):1.1*alpha2;
[RHO1,RHO2] = meshgrid(rho1,rho2);

% gammap's & Rp's
gamma1p = gamma1+s1.*RHO1;
gamma2p = gamma2+s2.*RHO2;
R1p = beta1./(gamma1p+mu);
R2p = beta2./(gamma2p+mu);

% prevalences
P1 = 1-1./R1p;
P2 = 1-1./R2p;
P12 = 1 - (gamma1p./R2p + gamma2p./R1p + mu)./(beta1+beta2-mu); %COINFECTION EQU
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

U12(imax,jmax)
rho1max = rho1(jmax);%???
rho2max = rho2(imax);

%max en rho2 quand rho1=rho1' (a verifier)
rho21 = (beta1+beta2-2*mu - (gamma1+s1*alpha1)/R2-gamma2)/(2*s2*((gamma1+s1*alpha1)/beta2+1)) - alpha1/2;
Ucase1 = (alpha1+alpha2)^2/4/beta2; %case simple

%max en rho1 quand rho2=rho2' (a verifier)
rho12 = (beta1+beta2-2*mu - (gamma2+s2*alpha2)/R1-gamma1)/(2*s1*((gamma2+s2*alpha2)/beta1+1)) - alpha2/2;
Ucase2 = (alpha1+alpha2).^2/(4*beta1); %case simple

%max en rho2 quand rho1=0
rho20 = (beta1+beta2-2*mu-gamma1/R2-gamma2/R1)/2/(s2/R1+s2*gamma1/beta2);
Ucase3 = ((beta1*beta2 - gamma1.*gamma2).^2)./(4*beta1*beta2*gamma1); %case simple

%max en rho1 quand rho2=0
rho10 = (beta1+beta2-2*mu-gamma1/R2-gamma2/R1)/2/(s1/R2+s1*gamma2/beta1);
Ucase4 = (beta1*beta2 - gamma1.*gamma2).^2./(4*beta1*beta2*gamma2); %case simple


%% Conditions pour le max de U12

% if machin est le max verifier que...
% if on a cette conditions, verifier que...
% a la maniere de la verification des conditions 
x1 = beta2*(-alpha1/(sqrt(beta2*gamma1))+1);
x2 = beta2*(alpha1/(sqrt(beta2*gamma1))+1);
x3 = beta1*(-alpha2/(sqrt(beta1*gamma2))+1);
x4 = beta1*(alpha2/(sqrt(beta1*gamma2))+1);
x5 = beta1*((gamma1-beta2)/sqrt(gamma1*beta1)+1);
x6 = beta1*((-gamma1+beta2)/sqrt(gamma1*beta1)+1);
x7 = beta2*((gamma2-beta1)/sqrt(beta2*gamma2)+1);
x8 = beta2*((-gamma2+beta1)/sqrt(beta2*gamma2)+1);

CondU1 = abs(Ucase1-maxU12)<0.1 & (beta1>=beta2) & (gamma1<x3 | gamma1>x4) & (gamma2<min(x5,x6) | gamma2>max(x5,x6));
CondU2 = abs(Ucase2-maxU12)<0.1 & (beta1<=beta2) & (gamma2<x1 | gamma2>x2) & (gamma1<min(x7,x8) | gamma1>max(x7,x8));
CondU3 = abs(Ucase3-maxU12)<0.1 & gamma2>=max([x1,x5,gamma1]) & gamma2<=min(x2,x6);
CondU4 = abs(Ucase4-maxU12)<0.1 & gamma1>=max([x7,x3,gamma2]) & gamma1<=min(x8,x4);


%Condtot = CondU1 | CondU2 | CondU3 | CondU4


%case where rho1'>rho2' (toute comparaison a u2 ne sert a rien)
condU12_11 = gamma2<x5 & gamma1<min(beta2,x3);
condU12_12 = gamma2>x6 & gamma1>x3;
condU12_13 = gamma2<x6 & gamma1<x3;
CondU12_1 = abs(Ucase1-maxU12)<0.05 & (condU12_11 | condU12_12 |condU12_13);
CondU12_3 = abs(Ucase3-maxU12)<0.05 & gamma2>=max([x5,gamma1]) & gamma2<=x6;
CondU12_4 = abs(Ucase4-maxU12)<0.05 & gamma1>=max([x3,gamma2]);


Condtot = CondU12_1 | CondU12_3 | CondU12_4;


%% Max of U (numerically)
[maxU, c2] = max(U(:));
[imax,jmax] = find(ismember(U, max(U(:))));

U(imax,jmax)
rho1max = rho1(jmax);%???
rho2max = rho2(imax);

Ucase2tilde = alpha2*(1-gamma1/beta1);
x9 = beta2/gamma2*(beta1-2*alpha2*(1-sqrt(1-beta1/beta2)));

CondU_4     = abs(Ucase4-maxU)<0.05 & (alpha1>rho10) & (alpha2<rho20) & (gamma1>max(0,x9));
condU_31    = (gamma2>x5 & gamma2<x6) & x5<x6;
condU_32    = (gamma2>x6 & gamma2<x5) & x5>x6;
CondU_3     = abs(Ucase3-maxU)<0.05 & (alpha1<rho10 & alpha2>rho20) & (condU_31 | condU_32);
condU2t1    = (alpha1>rho10 & alpha2<rho20); %& gamma1<x9 & x9>0; %par rapport à u4
condU2t2    = (alpha1<rho10 & alpha2<rho20); %& ??? %par rapport a u1
CondU_2t    = abs(Ucase2tilde-maxU)<0.5 & (condU2t1 | condU2t2);
condU_11    = alpha1<rho10 & alpha2>rho20 & (gamma2<x5 | gamma2>x6) & x5<x6; %par rapport a u3
condU_12    = alpha1<rho10 & alpha2>rho20 & (gamma2>x5 | gamma2<x6) & x5>x6; %par rapport a u3
condU_13    = alpha1<rho10 & alpha2<rho20; % & 
CondU_1     = abs(Ucase1-maxU)<0.05 & (condU_11 | condU_12 | condU_13);

u1max  = abs(Ucase1-maxU)<0.005;
u2tmax = abs(Ucase2-maxU)<0.005;
u3max  = abs(Ucase3-maxU)<0.005;
u4max  = abs(Ucase4-maxU)<0.005;
Condtot = CondU_1 | CondU_2t | CondU_3 | CondU_4


%% 05/10 conditions selon le shema du 02/10
x10 = beta1 - alpha2*(beta2/beta1*(2-sqrt(1-beta1/beta2))-1);
x11 = beta1 - alpha2*(beta2/beta1*(2+sqrt(1-beta1/beta2))-1);

condU_11    = alpha1<rho10 & alpha2>rho20 & (gamma2<min(x5,x6) | gamma2>max(x5,x6)); %par rapport a u3
condU_12    = alpha1<rho10 & alpha2<rho20 & beta1<=beta2 & (gamma1<min(x10,x11) | gamma1>min(x10,x11)) ; %par rapport a u3
condU_13    = alpha1<rho10 & alpha2<rho20 & beta1>=beta2; % & 
CondU_1     = abs(Ucase1-maxU)<0.05 & (condU_11 | condU_12 | condU_13);

condU2t1    = (alpha1>rho10 & alpha2<rho20); 
condU2t2    = (alpha1<rho10 & alpha2<rho20) & beta1<=beta2 & (gamma1>min(x10,x11) | gamma1<max(x10,x11)); 
CondU_2t    = abs(Ucase2tilde-maxU)<0.5 & (condU2t1 | condU2t2);

CondU_3     = abs(Ucase3-maxU)<0.05 & (alpha1<rho10 & alpha2>rho20) & (gamma2>min(x5,x6) | gamma2<max(x5,x6));

Condtot = CondU_1 | CondU_2t | CondU_3 












