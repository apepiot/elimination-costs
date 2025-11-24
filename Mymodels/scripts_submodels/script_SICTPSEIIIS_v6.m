%Looking for the reproduction number of HIV and of syphilis in the
%2-disease model SICTPxSEIIIS
n=1;p=1;m=0;
[nbCompartments,M,B,tabComp] = createODEsystem_v5(n,p,m);

%% Converting the matrix product to ODE system
[X,dX,eqn,F] = matToODE_v2(nbCompartments,M,B);
eqn.'
tabComp.X = X;

%% Write the ODE system in a text file
pathW = 'C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ODEsystems\';
%fileID = fopen([pathW,'ODE_SICTPSEIIIS_v6.txt'],'w');
for k=1:size(tabComp,1)
    fprintf(fileID,'%12s%1s\r\n ',eqn(k),';');    
end
fclose(fileID);

%%
%clear all;
syms mu b;
syms Lambdah thetah sigmah ph zetah
syms Lambdas sigmas taus thetas gamma1s gamma3s nus
syms rho_h rho_s
syms eta_s_prep eta_h_prep
syms eta_s_art
syms rho_hs 
syms VTunderART %boolean, 1 if VT can be practiced under ART, 0 otherwise

syms Y [35,1]
syms dY [35,1]

%Y([2:7:30,3:7:31,5:7:33,6:7:34])=0


dY(1) = Y(29)*(gamma3s + rho_s + rho_hs) - Y(1)*(Lambdah + Lambdas + mu) - b*(ph - 1) + Y(8)*(rho_s + rho_hs) + Y(15)*(rho_s + rho_hs) + Y(22)*(rho_s + rho_hs);
dY(2) = Lambdah*Y(1) - Y(2)*(Lambdas + mu + rho_h + rho_hs + sigmah) + Y(9)*rho_s + Y(16)*rho_s + Y(23)*rho_s + Y(30)*(gamma3s + rho_s);
dY(3) = Y(10)*rho_s + Y(17)*rho_s + Y(24)*rho_s - Y(3)*(Lambdas + mu + rho_h + rho_hs + thetah) + Y(2)*sigmah + Y(31)*(gamma3s + rho_s);
dY(4) = Y(11)*eta_s_prep + Y(18)*eta_s_prep + Y(25)*eta_s_prep + b*ph - Y(4)*(Lambdas + mu - Lambdah*(zetah - 1)) + Y(32)*(eta_s_prep + gamma3s);
dY(5) = Y(12)*eta_s_prep + Y(19)*eta_s_prep + Y(26)*eta_s_prep - Y(5)*(Lambdas + eta_h_prep + mu + sigmah) + Y(33)*(eta_s_prep + gamma3s) - Lambdah*Y(4)*(zetah - 1);
dY(6) = Y(13)*eta_s_prep + Y(20)*eta_s_prep + Y(27)*eta_s_prep + Y(5)*sigmah - Y(6)*(Lambdas + eta_h_prep + mu + thetah) + Y(34)*(eta_s_prep + gamma3s);
dY(7) = Y(14)*(eta_s_art + rho_s + VTunderART*rho_hs) + Y(21)*(eta_s_art + rho_s + VTunderART*rho_hs) + Y(28)*(eta_s_art + rho_s + VTunderART*rho_hs) + Y(5)*eta_h_prep + Y(9)*rho_hs + Y(10)*rho_hs + Y(16)*rho_hs + Y(17)*rho_hs + Y(23)*rho_hs + Y(24)*rho_hs + Y(30)*rho_hs + Y(31)*rho_hs + Y(3)*(rho_h + rho_hs + thetah) + Y(35)*(eta_s_art + gamma3s + rho_s + VTunderART*rho_hs) - Y(7)*(Lambdas + mu) + Y(6)*(eta_h_prep + thetah) + Y(2)*(rho_h + rho_hs);
dY(8) = Lambdas*Y(1) - Y(8)*(Lambdah + mu + rho_s + rho_hs + sigmas);
dY(9) = Lambdah*Y(8) + Lambdas*Y(2) - Y(9)*(mu + rho_h + rho_s + rho_hs + sigmah + sigmas);
dY(10) = Lambdas*Y(3) + Y(9)*sigmah - Y(10)*(mu + rho_h + rho_s + rho_hs + sigmas + thetah);
dY(11) = Lambdas*Y(4) - Y(11)*(eta_s_prep + mu + sigmas - Lambdah*(zetah - 1));
dY(12) = Lambdas*Y(5) - Y(12)*(eta_h_prep + eta_s_prep + mu + sigmah + sigmas) - Lambdah*Y(11)*(zetah - 1);
dY(13) = Lambdas*Y(6) + Y(12)*sigmah - Y(13)*(eta_h_prep + eta_s_prep + mu + sigmas + thetah);
dY(14) = Lambdas*Y(7) + Y(12)*eta_h_prep + Y(9)*rho_h + Y(13)*(eta_h_prep + thetah) - Y(14)*(eta_s_art + mu + rho_s + sigmas + VTunderART*rho_hs) + Y(10)*(rho_h + thetah);
dY(15) = Y(8)*sigmas - Y(15)*(Lambdah + mu + rho_s + rho_hs + taus);
dY(16) = Lambdah*Y(15) + Y(9)*sigmas - Y(16)*(mu + rho_h + rho_s + rho_hs + sigmah + taus);
dY(17) = Y(16)*sigmah + Y(10)*sigmas - Y(17)*(mu + rho_h + rho_s + rho_hs + taus + thetah);
dY(18) = Y(11)*sigmas - Y(18)*(eta_s_prep + mu + taus - Lambdah*(zetah - 1));
dY(19) = Y(12)*sigmas - Y(19)*(eta_h_prep + eta_s_prep + mu + sigmah + taus) - Lambdah*Y(18)*(zetah - 1);
dY(20) = Y(19)*sigmah + Y(13)*sigmas - Y(20)*(eta_h_prep + eta_s_prep + mu + taus + thetah);
dY(21) = Y(19)*eta_h_prep + Y(16)*rho_h + Y(14)*sigmas + Y(20)*(eta_h_prep + thetah) - Y(21)*(eta_s_art + mu + rho_s + taus + VTunderART*rho_hs) + Y(17)*(rho_h + thetah);
dY(22) = Y(15)*taus - Y(22)*(Lambdah + mu + rho_s + rho_hs + thetas);
dY(23) = Lambdah*Y(22) + Y(16)*taus - Y(23)*(mu + rho_h + rho_s + rho_hs + sigmah + thetas);
dY(24) = Y(23)*sigmah + Y(17)*taus - Y(24)*(mu + rho_h + rho_s + rho_hs + thetah + thetas);
dY(25) = Y(18)*taus - Y(25)*(eta_s_prep + mu + thetas - Lambdah*(zetah - 1));
dY(26) = Y(19)*taus - Y(26)*(eta_h_prep + eta_s_prep + mu + sigmah + thetas) - Lambdah*Y(25)*(zetah - 1);
dY(27) = Y(26)*sigmah + Y(20)*taus - Y(27)*(eta_h_prep + eta_s_prep + mu + thetah + thetas);
dY(28) = Y(26)*eta_h_prep + Y(23)*rho_h + Y(21)*taus + Y(27)*(eta_h_prep + thetah) - Y(28)*(eta_s_art + mu + rho_s + thetas + VTunderART*rho_hs) + Y(24)*(rho_h + thetah);
dY(29) = Y(22)*thetas - Y(29)*(Lambdah + gamma3s + mu + rho_s + rho_hs);
dY(30) = Lambdah*Y(29) + Y(23)*thetas - Y(30)*(gamma3s + mu + rho_h + rho_s + rho_hs + sigmah);
dY(31) = Y(30)*sigmah + Y(24)*thetas - Y(31)*(gamma3s + mu + rho_h + rho_s + rho_hs + thetah);
dY(32) = Y(25)*thetas - Y(32)*(eta_s_prep + gamma3s + mu - Lambdah*(zetah - 1));
dY(33) = Y(26)*thetas - Y(33)*(eta_h_prep + eta_s_prep + gamma3s + mu + sigmah) - Lambdah*Y(32)*(zetah - 1);
dY(34) = Y(33)*sigmah - Y(34)*(eta_h_prep + eta_s_prep + gamma3s + mu + thetah) + Y(27)*thetas;
dY(35) = Y(33)*eta_h_prep + Y(30)*rho_h + Y(28)*thetas - Y(35)*(eta_s_art + gamma3s + mu + rho_s + VTunderART*rho_hs) + Y(34)*(eta_h_prep + thetah) + Y(31)*(rho_h + thetah);



%% Endemic equilibrium
sol.syst = solve([dY==0],Y);

%% DFS
syms Y [35,1]
syms dY [35,1]

%at the DFS
Y([2:7:30,3:7:31,5:7:33,6:7:34])=0;
Y([8:35])=0;
Lambdah=0;
Lambdas=0;

eqn= dY==0;
sol = solve(eqn,[Y1,Y4,Y7]);

%% Recherche de R, pour HIV et syphilis
% Calculer M en mettant les taux rho_hg, rho_hg à 0 dans createODEsystem_v5
syms Y [35,1]
syms b mu

if (0)
    N = b/mu;
    syms betas betaIh betaCh
    %Lambdas = betas*sum(Y(15:35))/N;
    %Lambdah = betaIh*sum(Y([2:7:30,5:7:33]))/N + betaCh*sum(Y([3:7:31,6:7:34]))/N;
    F = zeros(35,1);  %new infections
    F = F + Lambdah*(tabComp(:,:).HIV=="I").*[0;Y(1:end-1)];
    F = F + Lambdah*(tabComp(:,:).HIV=="Ip").*[0;Y(1:end-1)];
    F = F + Lambdas*(tabComp(:,:).syph=="E").*[zeros(7,1);Y(1:end-7)];
    
    V =  -(M*Y - F);    %other rates
    Fsorted = [F([2,3,5:35]);F(1);F(4)];
    Vsorted = [V([2,3,5:35]);V(1);V(4)];
    Ysorted = [Y([2,3,5:35]);Y(1);Y(4)];
end
Ysorted = [Y([2,3,5:35]);Y(1);Y(4)];
[Fsorted,Vsorted] = FVevaluatedInY();

syms dFsorted [35,35];
syms dVsorted [35,35];
for i=1:35
    for j=1:35
        dFsorted(i,j) = diff(Fsorted(i),Ysorted(j));
        dVsorted(i,j) = diff(Vsorted(i),Ysorted(j));
    end
end


syms betaIh betaCh betas

[eigvalues] = Rp_SICTPSEIIIS_v5(betaIh,betaCh,thetah,sigmah,zetah,ph,...
                                      betas, sigmas, taus, thetas, gamma3s,...
                                      rho_h,rho_s,rho_hs,...
                                      eta_h_prep,eta_s_prep,eta_s_art,...
                                      VTunderART,mu,b);

%eig(dV)
%det(dV)
Vmoins1 = dV^(-1);
%eig(F*V) 


%% Test
%clear all;
syms mu b;
syms betas
syms Lambdah thetah sigmah ph zetah
syms Lambdas sigmas taus thetas gamma1s gamma3s nus
syms rho_h rho_s
syms eta_s_prep eta_h_prep
syms eta_s_art
syms rho_hs 
syms VTunderART %boolean, 1 if VT can be practiced under ART, 0 otherwise

%%
b=10;
pHIV=0.5; 

[paramTab,mu,vecAlphas] = sampleParameters_v3_extent(1,1,1,1,b,pHIV);
betaIh = paramTab{3}.betaI;
betaCh = paramTab{3}.betaC;
thetah = paramTab{3}.theta0;
sigmah = paramTab{3}.sigma;
zetah  = paramTab{3}.zeta;
betas  = paramTab{4}.beta;
sigmas = paramTab{4}.sigma;
taus   = paramTab{4}.tau;
thetas = paramTab{4}.theta;
gamma3s = paramTab{4}.gamma30;
rho_h   = 1.05;
rho_s   = 0.05;
rho_hs  = 0.1;
eta_h_prep = 0;
eta_s_prep = 1;
eta_s_art  = 1;
VTunderART = 1;

[eigvalues,F,V] = Rp_SICTPSEIIIS_v5(betaIh,betaCh,thetah,sigmah,zetah,pHIV,...
                              betas, sigmas, taus, thetas, gamma3s,...
                              rho_h,rho_s,rho_hs,...
                              eta_h_prep,eta_s_prep,eta_s_art,...
                              VTunderART,mu,b);
disp('----------------------------')                         
double(unique(abs(eigvalues)))
%vpa(max(abs(eigvalues),5))
disp(['R_s=',num2str(paramTab{4}.R0), ' (calculé avec le modele SEIIIS)'])
disp(['R_s=',num2str(double(eigS)), ' (calculé avec le modele SICTPxSEIIIS)'])

disp(['R_h=',num2str(paramTab{3}.R_prep_0), ' (calculé avec le modele SICTP)'])
disp(['R_h=',num2str(double(eigH)), ' (calculé avec le modele SICTPxSEIIIS)'])




%% Test du systeme d'ODE
Y0 = ones(35,1);
tspan = [1,1000];
res = ode45(@(t,Y) ODE_SICTPSEIIIS_v5(t,Y,betaIh,betaCh,thetah,sigmah,zetah,pHIV,...
                                         betas, sigmas, taus, thetas, gamma3s,...
                                         rho_h,rho_s,rho_hs,...
                                         eta_h_prep,eta_s_prep,eta_s_art,...
                                         VTunderART,mu,b),tspan, Y0);

%sum(res.y(:,end)) - b/mu

%syphilis
Ps = sum(res.y(8:35,end))/(b/mu)
%R0=1/(1-Ps)
%1-1./paramTab{4}.R0

%sictp/hiv
Ph = sum(res.y([2:7:30,5:7:33,3:7:31,6:7:34],end))/(b/mu)
%[Rp,Lambdap,alpha,Ptot,Pun] = Rp_SICTP(betaIh,betaCh,thetah,0,sigmah,zetah,eta_h_prep,pHIV,mu,b,rho_h);
%Pun
disp('----------------------------')                         


b=10;
pHIV=0.48; 

% [paramTab,mu,vecAlphas] = sampleParameters_v3_extent(1,1,1,1,b,pHIV);
% betaIh = paramTab{3}.betaI;
% betaCh = paramTab{3}.betaC;
% thetah = paramTab{3}.theta0;
% sigmah = paramTab{3}.sigma;
 zetah  = 1%paramTab{3}.zeta;%0.5174
% betas  = paramTab{4}.beta;
% sigmas = paramTab{4}.sigma;
% taus   = paramTab{4}.tau;
% thetas = paramTab{4}.theta;
% gamma3s = paramTab{4}.gamma30;
% nus = paramTab{4}.nu
rho_h   = 0.1%.05;
rho_s   = 0%0.%05;
rho_hs  = 0.1;
eta_h_prep = 0;
eta_s_prep = 0;
eta_s_art  = 00;
VTunderART = 1;

[Rp_s,~,~] = Rp_SEIIIS_v4(betas,sigmas,taus,nus,0,thetas,gamma3s,mu,b,rho_s);
[R_prep,~,~,Ptot_prep,Pun_prep] = Rp_SICTP(betaIh,betaCh,...
    thetah,0,sigmah,zetah,eta_h_prep,...
    pHIV,mu,b,rho_h);

[eigvalues,F,V] = Rp_SICTPSEIIIS_v5(betaIh,betaCh,thetah,sigmah,zetah,pHIV,...
                              betas, sigmas, taus, thetas, gamma3s,...
                              rho_h,rho_s,rho_hs,...
                              eta_h_prep,eta_s_prep,eta_s_art,...
                              VTunderART,mu,b);
disp('----------------------------')                         
double(unique(abs(eigvalues)))
%vpa(max(abs(eigvalues),5))
disp(['R_s=',num2str(Rp_s), ' (calculé avec le modele SEIIIS)'])
%disp(['R_s=',num2str(double(eigS)), ' (calculé avec le modele SICTPxSEIIIS)'])

disp(['R_h=',num2str(R_prep), ' (calculé avec le modele SICTP)'])
%disp(['R_h=',num2str(double(eigH)), ' (calculé avec le modele SICTPxSEIIIS)'])




%% Test du systeme d'ODE
Y0 = ones(35,1);
tspan = [1,1000];
res = ode45(@(t,Y) ODE_SICTPSEIIIS_v5(t,Y,betaIh,betaCh,thetah,sigmah,zetah,pHIV,...
                                         betas, sigmas, taus, thetas, gamma3s,...
                                         rho_h,rho_s,rho_hs,...
                                         eta_h_prep,eta_s_prep,eta_s_art,...
                                         VTunderART,mu,b),tspan, Y0);

%sum(res.y(:,end)) - b/mu

%syphilis
Ps = sum(res.y(8:35,end))/(b/mu)
%R0=1/(1-Ps)
%1-1./paramTab{4}.R0

%sictp/hiv
Ph = sum(res.y([2:7:30,5:7:33,3:7:31,6:7:34],end))/(b/mu)
%[Rp,Lambdap,alpha,Ptot,Pun] = Rp_SICTP(betaIh,betaCh,thetah,0,sigmah,zetah,eta_h_prep,pHIV,mu,b,rho_h);
%Pun
disp('----------------------------')                         




