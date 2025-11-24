%script to plot rhohat for the SIRxSIS2 model 
clear all
pasC=0.001;
% parameters
%Paris
mu = 1/35;
%PHIV = 16.1/100*9/100;
%PHIV=10/100;
R1 = 2.5;
gamma1 = 1/7;
beta1 = R1*(gamma1+mu);
%beta1 = mu./(mu./(gamma1+mu)-PHIV);
%R1 = beta1/(gamma1+mu);
s1=1;
alpha1=beta1/s1*(1-1/R1);

%Pch = 5.3/100;
%R2 = 1/(1-Pch);
R2 = 3;
gamma2 = 1/3;
beta2=R2*(gamma2+mu);
s2=1;

%Ps = 6.6/100;
%R3 = 1/(1-Ps);
R3 = 2.8;
gamma3=1/3;
beta3=R3*(gamma3+mu);
s3=1;

b=10;

%%a supp
mu = 1/35; R1 = 2.5; gamma1 = 1/7;
beta1=R1*(gamma1+mu);
s1=1;
alpha1=beta1/s1*(1-1/R1);

R2 = 3;
gamma2 = 1/3;
beta2=R2*(gamma2+mu);
s2=1;
R3 = 2.8;
gamma3=1/3;
beta3=R3*(gamma3+mu);
s3=1;
b=10;
% fin de a supp

alpha1=beta1/s1*(1-1/R1);
alpha2=beta2/s2*(1-1/R2);
alpha3=beta3/s3*(1-1/R3);
maxalpha=max([alpha1,alpha2,alpha3]);

%% c thresholds
c11 = -mu/beta1*(1-1/R1);
%c01 = mu/beta1*(1-R1);
c10 = mu/beta1*(R1-1);
c22 = 1/R2-1;
c02 = 1-1/R2;
c33 = 1/R3-1;
c03 = 1-1/R3;
[U120,c012] = U12_SIRSIS7(0,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,0); %c120 = P12(0)
[U130,c013] = U12_SIRSIS7(0,beta1,beta3,gamma1,gamma3,s1,s3,b,mu,0); %c130 = P13(0)
[U230,c023] = U12_SIS2(0,beta3,beta2,gamma3,gamma2,s3,s2,b,mu,0); %need to be checked
[U120,P120,c112] = U12_SIRSIS7(alpha1,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,0);
[U120,P120,c212] = U12_SIRSIS7(alpha2,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,0);
[U130,P130,c113] = U12_SIRSIS7(alpha1,beta1,beta3,gamma1,gamma3,s1,s3,b,mu,0);
[U130,P130,c313] = U12_SIRSIS7(alpha3,beta1,beta3,gamma1,gamma3,s1,s3,b,mu,0);
[U230,P230,c223] = U12_SIS2(alpha2,beta2,beta3,gamma2,gamma3,s2,s3,b,mu,0);
[U230,P230,c323] = U12_SIS2(alpha3,beta2,beta3,gamma2,gamma3,s2,s3,b,mu,0);
[U1230,c0123,dU123] = U123_SIRSIS2(0,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,0);
[U1230,P1230,c1123] = U123_SIRSIS2(alpha1,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,0);
[U1230,P1230,c2123] = U123_SIRSIS2(alpha2,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,0);
[U1230,P1230,c3123] = U123_SIRSIS2(alpha3,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,0);

clear U120 U130 U230 P120 P130 P230 P1230 U1230

%% argmax U
vecC=-1:pasC:1; l=length(vecC);
vecRhohat=zeros(1,l);
vecRhohat1=zeros(1,l);vecRhohat2=zeros(1,l);vecRhohat3=zeros(1,l);
vecRhohat12=zeros(1,l);vecRhohat23=zeros(1,l);vecRhohat13=zeros(1,l);
vecRhohat123=zeros(1,l); 
vecRhohatSIRSIS12=zeros(1,l); vecRhohatSIRSIS13=zeros(1,l);vecRhohatSISSIS23=zeros(1,l);
vecRhotothat=zeros(1,l);
i=1;
for c=vecC
    alpha1 = beta1/s1*(1-1/R1);
    alpha2 = beta2/s2*(1-1/R2);
    alpha3 = beta3/s3*(1-1/R3);
    maxalpha=max([alpha1,alpha2,alpha3]);
    minalpha=min([alpha1,alpha2,alpha3]);
    rho1hat = beta1/R1/s1*(sqrt(R1*mu/(mu+beta1*c))-1);
    %     if (imag(rho1hat))
    %         rho1hat=0;
    %     end
    rho1hat = rho1hat*(c>c11 & c<c10) + alpha1*(c<=c11);
    rho2hat = min(max(alpha2/2-c*beta2/(2*s2),0),alpha2);
    rho3hat = min(max(alpha3/2-c*beta3/(2*s3),0),alpha3);
    
    options = optimset('Display','off');
    fun12 = @(rho) -U12_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    rho12hat = min(max(fminsearch(fun12,0,options),0),max(alpha1,alpha2)); %argmax U12
    fun13 = @(rho) -U12_SIRSIS7(rho,beta1,beta3,gamma1,gamma3,s1,s3,b,mu,c);
    rho13hat = min(max(fminsearch(fun13,0,options),0),max(alpha1,alpha3));
    fun23 = @(rho) -U12_SIS2(rho,beta3,beta2,gamma3,gamma2,s3,s2,b,mu,c);
    rho23hat = min(max(fminsearch(fun23,0,options),0),max(alpha2,alpha3));
    fun123 = @(rho) -U123_SIRSIS2(rho,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c);
    rho123hat = min(max(fminsearch(fun123,0,options),0),max([alpha1,alpha2,alpha3]));
    
    %argmaxU in total
    list_rhohat = max([rho1hat,rho2hat,rho3hat,rho12hat,rho13hat,rho23hat,rho123hat],0);
    U_list = U_SIRSIS2(list_rhohat,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c);
    [Umax,imax] = max(U_list);
    rhohat = list_rhohat(imax);
    
    %argmax U_SIRxSIS 12
    list_rhohat12 = max([rho1hat,rho2hat,rho12hat],0);
    U_list12 = U_SIRSIS7(list_rhohat12,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,c);
    [Umax12,imax12] = max(U_list12);
    rhohat_SIRSIS12 = list_rhohat12(imax12);
    
    %argmax U_SIRxSIS 13
    list_rhohat13 = max([rho1hat,rho3hat,rho13hat],0);
    U_list13 = U_SIRSIS7(list_rhohat13,beta1,beta3,gamma1,gamma3,s1,s3,b,mu,c);
    [Umax13,imax13] = max(U_list13);
    rhohat_SIRSIS13 = list_rhohat13(imax13);
       
    %argmax U_SISxSIS 23
    list_rhohat23 = max([rho2hat,rho3hat,rho23hat],0);
    U_list23 = U_SIS2(list_rhohat23,beta3,beta2,gamma3,gamma2,s3,s2,b,mu,c);
    [Umax23,imax23] = max(U_list23);
    rhohat_SISSIS23 = list_rhohat23(imax23);
    
    vecRhohat1(i) = rho1hat;
    vecRhohat2(i) = rho2hat;
    vecRhohat3(i) = rho3hat;
    vecRhohat12(i) = rho12hat;
    vecRhohat23(i) = rho23hat;
    vecRhohat13(i) = rho13hat;
    vecRhohat123(i) = rho123hat;
    vecRhohat(i) = rhohat;
    
    vecRhohatSIRSIS12(i) = rhohat_SIRSIS12;
    vecRhohatSIRSIS13(i) = rhohat_SIRSIS13;
    vecRhohatSISSIS23(i) = rhohat_SISSIS23;
    
    %for comparison
    fun = @(rho) -U_SIRSIS2(rho,beta1,beta2,beta3,gamma1,gamma2,gamma3,s1,s2,s3,b,mu,c);
    rhoTOThat = min(fminsearch(fun,0,options),max([alpha1,alpha2,alpha3])); %if there are two peaks, is it gonna take the global max ?
    vecRhotothat(i) = max(rhoTOThat,0);
    i=i+1;
end
%% c thresholds
vecAlpha = [alpha1,alpha2,alpha3];
[alphai,ii]=max(vecAlpha);
[alphak,kk]=min(vecAlpha);
jj = 6-ii-kk;
alphaj = vecAlpha(jj);

%cii
vecCii = [c11,c22,c33]; listCii = {' $c_1^1=$',' $c_2^2=$',' $c_3^3=$'}; 
cii = vecCii(ii);
%cjij
vecCjij(:,:,1) = [0,c112,c113;c112,0,0;c113,0,0];
vecCjij(:,:,2) = [0,c212,0;c212,0,c223;0,c223,0];
vecCjij(:,:,3) = [0,0,c313;0,0,c323;c313,c323,0];
cjij = vecCjij(ii,jj,jj);
%listCjij...
%ck123
vecCk123 = [c1123,c2123,c3123];
ck123 = vecCk123(kk);

%cswitch
%Il doit y avoir 2 switches dans ce modèle. Un entre rhoihat et rhoijhat et
%un entre rhoijhat et rhoijkhat.
tabRhohats = [vecRhohat1;vecRhohat2;vecRhohat3;vecRhohat12;vecRhohat13;vecRhohat23];
vecRhohati = tabRhohats(ii,:);
vecRhohatij = tabRhohats(ii+jj+1,:); %voir les notes du 23/02/22 (fin de page du verso)

cs123toij = max(vecC(vecRhohat==vecRhohatij & vecRhohatij~=0));%ok (affiner le c)
csijtoi = max(vecC(vecRhohat==vecRhohati & vecRhohati~=0));%ok (affiner le c)
cs123toi = max(vecC(vecRhohat==vecRhohatij & vecRhohatij~=0));%ok (affiner le c)
cs3dis = max(vecC(vecRhohat==maxalpha)); %when 3 diseases are eliminated
cs2dis = max(vecC(vecRhohat<maxalpha & vecRhohat>=alphak & vecRhohat>=alphaj))+pasC/2; %when 2 diseases are eliminated
cs1dis = max(vecC(vecRhohat>=minalpha))+pasC/2; %when only one disease is eliminated

%%
figure
p_rho1 = plot(vecC,vecRhohat1,':')
hold on
p_rho2  = plot(vecC,vecRhohat2,':')
p_rho3  = plot(vecC,vecRhohat3,':')
p_rho12 = plot(vecC,vecRhohatSIRSIS12,'--')
p_rho23 = plot(vecC,vecRhohatSISSIS23,'--')
p_rho13 = plot(vecC,vecRhohatSIRSIS13,'--')
%plot(vecC,vecRhohat12,'--',vecC,vecRhohat13,'--',vecC,vecRhohat23,'--')
%plot(vecC,vecRhohat123)
plot(vecC,vecRhohat,'Linewidth',1.5)
%legend('$\hat\rho_1$','$\hat\rho_2$','$\hat\rho_3$','$\hat\rho_{12}$','$\hat\rho_{23}$','$\hat\rho_{13}$','$\hat\rho_{123}$','$\hat\rho$','Interpreter','latex')

title([{'Utility U of the SIR$\times$SIS$^2$ model'},...
    {['$\beta_1=$',num2str(round(beta1,2)), ' $\beta_2=$',num2str(round(beta2,2)),' $\beta_3=$',num2str(round(beta3,2)),' $\gamma_1(0)^{-1}=$',num2str(round(1/gamma1,2)),' $\gamma_2(0)^{-1}=$',...
    num2str(round(1/gamma2,2)),' $\gamma_3(0)^{-1}=$',num2str(round(1/gamma3,2)), ' $s_1=$', num2str(s1),' $s_2=$', num2str(s2),' $s_3=$', num2str(s3), ' $\mu^{-1}=$', num2str(1/mu)]}, ...
    {[' $\mathtt R_1(0)=$' num2str(round(R1,2)), ' $\mathtt R_2(0)=$' num2str(round(R2,2)),' $\mathtt R_3(0)=$' num2str(round(R3,2)),...
    ' $\rho_1\prime=$',num2str(round(alpha1,2)),' $\rho_2\prime=$',num2str(round(alpha2,2)),' $\rho_3\prime=$',num2str(round(alpha3,2))]},...
    {['$i=$',num2str(ii),' $j=$',num2str(jj),' $k=$',num2str(kk),' $c_i^i=$',num2str(round(cii,2)),' $c_j^{ij}=$',num2str(round(cjij,2)),' $c_k^{123}=$',num2str(round(ck123,2)),...
    ' $c_0^{123}=$',num2str(round(c0123,2)),' $c_{\hat\rho=\hat\rho_i}=$',num2str(round(csijtoi,2)),' $c_{\hat\rho=\hat\rho_{ij}}=$',num2str(round(cs123toij,2))]}],...
    'Interpreter','latex')
%     {[' $c_0^1=$',num2str(round(c01,2)),' $c_1^1=$',num2str(round(c11,2)),' $c_0^2=$',num2str(round(c02,2)),' $c_2^2=$',num2str(round(c22,2)),...
%     ' $c_0^3=$',num2str(round(c03,2)),' $c_3^3=$',num2str(round(c33,2))]},...
%     {[' $c_0^{12}=$',num2str(round(c012,2)),' $c_1^{12}=$',num2str(round(c112,2)), ' $c_2^{12}=$',num2str(round(c212,2)),...
%     ' $c_0^{23}=$',num2str(round(c023,2)),' $c_2^{23}=$',num2str(round(c223,2)), ' $c_3^{23}=$',num2str(round(c323,2)),...
%     ' $c_0^{13}=$',num2str(round(c013,2)),' $c_1^{13}=$',num2str(round(c113,2)), ' $c_3^{13}=$',num2str(round(c313,2))]}
xlabel("$c$","fontweight","bold",'Interpreter','latex')
ylabel("$\hat \rho$","fontweight","bold",'Interpreter','latex')
limy=1.2*maxalpha;
ylim([0,limy])


%% ajout de 3 zones sur le graphe: epidemie when rho=0, 1 disease eliminated, 2 diseases eliminated, 3 diseases eliminated
hold on

% 0, 1, 2, 3 diseases eliminated zones
if(size(cs1dis,2)~=0 && size(cs2dis,2)~=0 && cs1dis>cs2dis && cs2dis>cs3dis)
    a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
    a=area([cs3dis cs3dis cs2dis cs2dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.3;  
    a=area([cs2dis cs2dis cs1dis cs1dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.15;  
    a=area([cs1dis cs1dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
    legend('$\hat\rho_1$','$\hat\rho_2$','$\hat\rho_3$','$\hat\rho_{1\times 2}$',...
    '$\hat\rho_{2\times 3}$','$\hat\rho_{1 \times 3}$','$\hat\rho$',...
    '1,2 \& 3 eliminated', [num2str(jj),' \& ', num2str(kk), ' eliminated'],...
    ['only ', num2str(kk), ' eliminated'],...
    'no disease eliminated',...
    'Interpreter','latex')

% 0, 3 diseases eliminated zones
elseif(size(cs1dis,2)==0 && size(cs2dis,2)==0)
    a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
    a=area([cs3dis cs3dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
    legend('$\hat\rho_1$','$\hat\rho_2$','$\hat\rho_3$','$\hat\rho_{1\times 2}$',...
    '$\hat\rho_{2\times 3}$','$\hat\rho_{1 \times 3}$','$\hat\rho$',...
    '1,2 \& 3 eliminated', 'no disease eliminated',...
    'Interpreter','latex')

% 0, 3 diseases eliminated
elseif(size(cs1dis,2)~=0 && size(cs2dis,2)==0 && cs3dis>=cs1dis)
    a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
    a=area([cs3dis cs3dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
    legend('$\hat\rho_1$','$\hat\rho_2$','$\hat\rho_3$','$\hat\rho_{1\times 2}$',...
    '$\hat\rho_{2\times 3}$','$\hat\rho_{1 \times 3}$','$\hat\rho$',...
    '1,2 \& 3 eliminated', 'no disease eliminated',...
    'Interpreter','latex')

% 0, 2, 3 diseases eliminated
elseif(size(cs1dis,2)==0 && size(cs2dis,2)~=0 && cs3dis<cs2dis || size(cs1dis,2)~=0 && size(cs2dis,2)~=0 && cs1dis==cs2dis && cs3dis<cs2dis)
    a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
    a=area([cs3dis cs3dis cs2dis cs2dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.3;  
    a=area([cs2dis cs2dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
    legend('$\hat\rho_1$','$\hat\rho_2$','$\hat\rho_3$','$\hat\rho_{1\times 2}$',...
    '$\hat\rho_{2\times 3}$','$\hat\rho_{1 \times 3}$','$\hat\rho$',...
    '1,2 \& 3 eliminated', [num2str(jj),' \& ', num2str(kk), ' eliminated'], 'no disease eliminated',...
    'Interpreter','latex')

% 0, 1, 3 diseases eliminated
elseif(size(cs1dis,2)~=0 && size(cs2dis,2)==0 && cs3dis<cs1dis)
    a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
    a=area([cs3dis cs3dis cs1dis cs1dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.15;  
    a=area([cs1dis cs1dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
    legend('$\hat\rho_1$','$\hat\rho_2$','$\hat\rho_3$','$\hat\rho_{1\times 2}$',...
    '$\hat\rho_{2\times 3}$','$\hat\rho_{1 \times 3}$','$\hat\rho$',...
    '1,2 \& 3 eliminated', ['only ', num2str(kk),' eliminated'], 'no disease eliminated',...
    'Interpreter','latex')
end

