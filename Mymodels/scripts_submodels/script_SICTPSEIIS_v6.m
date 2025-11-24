% looking for the reproduction number of SEIIS when there is prep
% there is not script_SICTPSEIIS_v1
% i chose script_SICTPSEIIS_v2 to match the v2 of the global model




%Looking for the reproduction number of HIV and of syphilis in the
%2-disease model SICTPxSEIIIS
n=1;p=0;m=1;
[nbCompartments,M,B,tabComp] = createODEsystem_v5(n,p,m);


%% Converting the matrix product to ODE system
[X,dX,eqn,F] = matToODE_v2(nbCompartments,M,B);
eqn.'
tabComp.X = X;

%% Write the ODE system in a text file
pathW = 'C:\Users\Moi\Documents\IPLESP\These\Codes\multi-voluntary-testing\Mymodels\ODEsystems\';
fileID = fopen([pathW,'ODE_SICTPSEIIS_v6.txt'],'w');
for k=1:size(tabComp,1)
    fprintf(fileID,'%12s%1s\r\n ',eqn(k),';');    
end
fclose(fileID);


%% Test du systeme d'ODE
b=10;
pHIV=0.5; 

[paramTab,mu,~] = sampleParameters_v3_extent(1,1,1,1,b,pHIV);
betaIh = paramTab{3}.betaI;
betaCh = paramTab{3}.betaC;
thetah = paramTab{3}.theta0;
sigmah = paramTab{3}.sigma;
zetah  = paramTab{3}.zeta;
betaCt = paramTab{1}.beta;
sigmaCt= paramTab{1}.sigma;
nuCt   = paramTab{1}.nu;
epsCt  = paramTab{1}.eps;
gammaCt = paramTab{1}.gamma;
rho_h   = 0.1;
rho_c   = 0.1;
rho_hc  = 0.;
eta_h_prep = 1;
eta_c_prep = 1;
eta_c_art  = 1;
VTunderART = 1;

Y0 = ones(28,1);
tspan = [1,1000];
res = ode45(@(t,Y) ODE_SICTPSEIIS_v5(t,Y,betaIh,betaCh,thetah,sigmah,zetah,ph,...
                                         betaCt, gammaCt, nuCt, epsCt, sigmaCt,...
                                         rho_h,rho_c,rho_hc,...
                                         eta_h_prep,eta_c_prep,eta_c_art,...
                                         VTunderART,mu,b),tspan, Y0);
                                     
Ph = sum(res.y([2:7:23,5:7:26,3:7:24,6:7:27],end))/(b/mu);
Pc = sum(res.y([8:28],end))/(b/mu);
1/(1-Pc)

[Rp_c,~,~] = Rp_SEIIS_v4(betaCt,nuCt,epsCt,sigmaCt,gammaCt,mu,b,rho_c);
[R_prep,~,~,Ptot_prep,Pun_prep] = Rp_SICTP(betaIh,betaCh,...
    thetah,0,sigmah,zetah,eta_h_prep,pHIV,mu,b,rho_h);
disp(['R_h=',num2str(R_prep), ' (calculé avec le modele SICTP)'])
disp(['R_c=',num2str(Rp_c), ' (calculé avec le modele SEIIS)'])

                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


