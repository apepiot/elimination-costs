%script to plot rhohat for the SIRxSIS2 model 
clear all; close all
pasC = 0.01; plotSIRSIS=true;plotSIS2=false;
% parameters
mu = 1/35;

R1 = 2;
gamma1 = 365/14;
beta1 = R1*(gamma1+mu);
P1 = mu/beta1*(R1-1);
s1 = 1;

R2 = 2.5;
gamma2 = 12/1; 
beta2 = R2*(gamma2+mu);
s2=1;
P2=1-1/R2;

R3 = 3.1;
gamma3 = 365/21; %10 days
beta3=R3*(gamma3+mu);
s3=1;

b=10;

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
vecC=-1.5:pasC:1; l=length(vecC);
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
    fun23 = @(rho) -U12_SIS2(rho,beta3,beta2,gamma3,gamma2,s3,s2,b,mu,c); %
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


%% plot SIRxSIS (e.g. HIV, chlamydia)
if(plotSIRSIS)
    close all
    figure(1)
    p_rho1 = plot(vecC,vecRhohat1,':','Linewidth',2)
    hold on
    p_rho2  = plot(vecC,vecRhohat2,':','Linewidth',2)
    p_rho3  = plot(vecC,vecRhohat3,':','Linewidth',2)
    p_rho12 = plot(vecC,vecRhohatSIRSIS12,'--','Linewidth',2)
    p_rho23 = plot(vecC,vecRhohatSISSIS23,'--','Linewidth',2)
    p_rho13 = plot(vecC,vecRhohatSIRSIS13,'--','Linewidth',2)
    p_rho = plot(vecC,vecRhohat,'Linewidth',1.5,'Linewidth',2.5)

    title('Fréquence de dépistage volontaire atteinte en fonction du cout c')
    xlabel("$c$","fontweight","bold",'Interpreter','latex')
    ylabel("frequence de depistage volontaire","fontweight","bold",'Interpreter','latex')
    limy=1.2*maxalpha;
    ylim([0,limy])

    % plot de R
    R1rho1 = beta1./(gamma1+mu+s1*vecRhohat1);R2rho2 = beta2./(gamma2+mu+s2*vecRhohat2);
    R1rho12 = beta1./(gamma1+mu+s1*vecRhohatSIRSIS12);R2rho12 = beta2./(gamma2+mu+s2*vecRhohatSIRSIS12);
    figure(2)
    p_R1_1 = plot(vecC,R1rho1,':','Linewidth',1.5);
    hold on
    p_R2_2 = plot(vecC,R2rho2,':','Linewidth',1.5);
    p_R1_12 = plot(vecC,R1rho12,'--','Linewidth',1.5,'Color', 'b');
    p_R2_12 = plot(vecC,R2rho12,'--','Linewidth',1.5, 'Color', 'r');
    title(['Nombres de reproduction effectifs en tenant compte'],...
        ['de la fréquence de dépistage volontaire atteinte en fonction du coût c'])
    xlabel("$c$","fontweight","bold",'Interpreter','latex')
    ylabel("$R(\hat \gamma)$","fontweight","bold",'Interpreter','latex')
    limyR=1.05*max(R1,R2);
    ylim([1,limyR])

    %plot de la prevalence
    P1_1 = mu/beta1.*(R1rho1-1); P2_2 = 1-1./R2rho2;
    P1_12 = max(mu/beta1.*(R1rho12-1),0); P2_12 = max(1-1./R2rho12,0);
    figure(3)
    p_P1_1 = plot(vecC,P1_1,':','Linewidth',1.5)
    hold on
    p_P2_2 = plot(vecC,P2_2,':','Linewidth',1.5)
    p_P1_12  = plot(vecC,P1_12,'--','Linewidth',1.5,'Color', 'b')
    p_P2_12  = plot(vecC,P2_12,'--','Linewidth',1.5, 'Color', 'r')
    title(['Prévalences en tenant compte de la fréquence '],...
        ['de dépistage volontaire atteinte en fonction du coût c'])
    xlabel("$c$","fontweight","bold",'Interpreter','latex')
    ylabel("$Prevalence(\hat \gamma)$","fontweight","bold",'Interpreter','latex')
    limyP=1.05*max([P1_1,P2_2]);
    ylim([0,limyP])
    %

    % SIRxSIS : 2 diseases 3 zones at max /rhohat
    figure(1)
    delete(p_rho3);delete(p_rho);delete(p_rho13);delete(p_rho23)
    if(0) %if 0,1,2 (code a modifier en consequences)
        a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.4;  
        a=area([cs3dis cs3dis cs2dis cs2dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.15;  
        a=area([cs2dis cs2dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
        legend('freq. de dep. vol. VIH ','freq. de dep. vol. chlamydia ',...
        'freq. de dep. vol. VIH \& chlamydia',...
        'elimination de VIH \& chlamydia',...
        'seulement VIH est elimine', 'HIV \& chlamydia persistent',...
        'Interpreter','latex')  

        figure(2) %R's
        a=area([vecC(1) vecC(1) cs3dis cs3dis],[1,limyR,limyR,1 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.4;  
        %uistack(a,'bottom') %met derriere
        a=area([cs3dis cs3dis cs2dis cs2dis],[1,limyR,limyR,1 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.15;  
        %uistack(a,'bottom')
        a=area([cs2dis cs2dis vecC(end) vecC(end)],[1,limyR,limyR,1],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
        %uistack(a,'bottom')
        legend('$R_{vih}$','$R_{ch}$','$R_{vih}$','$R_{ch}$',...
        'elimination de HIV \& chlamydia',...
        'seulement chlamydia est eliminee', 'HIV \& chlamydia persistent',...
        'Interpreter','latex')

        figure(3) %prevalence
        a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limyP,limyP,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.4;  
        %uistack(a,'bottom') %met derriere
        a=area([cs3dis cs3dis cs2dis cs2dis],[0,limyP,limyP,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.15;  
        %uistack(a,'bottom')
        a=area([cs2dis cs2dis vecC(end) vecC(end)],[0,limyP,limyP,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
        %uistack(a,'bottom')
        legend('Prevalence VIH','Prevalence chlamydia','Prevalence VIH ','Prevalence chlamydia',...
        'elimination de HIV \& chlamydia',...
        'seulement chlamydia est eliminee', 'HIV \& chlamydia persistent',...
        'Interpreter','latex')
    end
    if(0) %0, 2 hiv, chlamydia
        a=area([vecC(1) vecC(1) cs2dis cs2dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
        a=area([cs2dis cs2dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
        legend('$\hat\gamma_{hiv}$ (specifique VIH)','$\hat\gamma_{ch}$ (specifique chlamydia)','$\hat\gamma_{hiv\times ch}$ (combine HIV et chlamydia)',...
        'elimination of HIV \& chlamydia',...
        'HIV \& chlamydia persistent',...
        'Interpreter','latex')
    end
    if(1) %0, 2 chlamydia,hiv
        %c2323 = -0.176;
        a=area([vecC(1) vecC(1) c2323 c2323],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
        a=area([c2323 c2323 vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
        legend('freq. de dep. vol. VIH ','freq. de dep. vol. chlamydia ',...
        'freq. de dep. vol. VIH \& chlamydia',...
        'elimination de VIH \& chlamydia',...
        'HIV \& chlamydia persistent',...
        'Interpreter','latex')
    end
end %end plot SIRxSIS

%% plot SISxSIS (e.g. syphilis, chlamydia)
if(plotSIS2)
    close all
    figure(4)
    p_rho1 = plot(vecC,vecRhohat1,':','Linewidth',2)
    hold on
    p_rho2  = plot(vecC,vecRhohat2,':','Linewidth',2)
    p_rho3  = plot(vecC,vecRhohat3,':','Linewidth',2)
    p_rho12 = plot(vecC,vecRhohatSIRSIS12,'--','Linewidth',2)
    p_rho23 = plot(vecC,vecRhohatSISSIS23,'--','Linewidth',2)
    p_rho13 = plot(vecC,vecRhohatSIRSIS13,'--','Linewidth',2)
    p_rho = plot(vecC,vecRhohat,'Linewidth',1.5,'Linewidth',2.5)

    title('Fréquence de dépistage volontaire atteinte en fonction du cout c')
    xlabel("$c$","fontweight","bold",'Interpreter','latex')
    ylabel("frequence de depistage volontaire","fontweight","bold",'Interpreter','latex')
    limy=1.2*maxalpha;
    ylim([0,limy])

    % plot de R
    R3rho3 = beta3./(gamma3+mu+s3*vecRhohat3);R2rho2 = beta2./(gamma2+mu+s2*vecRhohat2);
    R3rho23 = beta3./(gamma3+mu+s3*vecRhohatSISSIS23);R2rho23 = beta2./(gamma2+mu+s2*vecRhohatSISSIS23);
    figure(5)
    p_R3_3 = plot(vecC,R3rho3,':','Linewidth',1.5);
    hold on
    p_R2_2 = plot(vecC,R2rho2,':','Linewidth',1.5);
    p_R3_23 = plot(vecC,R3rho23,'--','Linewidth',1.5,'Color', 'b');
    p_R2_23 = plot(vecC,R2rho23,'--','Linewidth',1.5, 'Color', 'r');
    title(['Nombres de reproduction effectifs en tenant compte'],...
        ['de la fréquence de dépistage volontaire atteinte en fonction du coût c'])
    xlabel("$c$","fontweight","bold",'Interpreter','latex')
    ylabel("$R(\hat \gamma)$","fontweight","bold",'Interpreter','latex')
    limyR=1.05*max(R3,R2);
    ylim([1,limyR])

    %plot de la prevalence
    P3_3 = 1-1./R3rho3; P2_2 = 1-1./R2rho2;
    P2_23 = max(1-1./R2rho23,0); P3_23 = max(1-1./R3rho23,0);
    figure(6)
    p_P2_2 = plot(vecC,P2_2,':','Linewidth',1.5)
    hold on
    p_P3_3 = plot(vecC,P3_3,':','Linewidth',1.5)
    p_P2_23  = plot(vecC,P2_23,'--','Linewidth',1.5,'Color', 'b')
    p_P3_23  = plot(vecC,P3_23,'--','Linewidth',1.5, 'Color', 'r')
    title(['Prévalences en tenant compte de la fréquence '],...
        ['de dépistage volontaire atteinte en fonction du coût c'])
    xlabel("$c$","fontweight","bold",'Interpreter','latex')
    ylabel("$Prevalence(\hat \gamma)$","fontweight","bold",'Interpreter','latex')
    limyP=1.05*max([P2_2,P3_3]);
    ylim([0,limyP])
    %

    %% SISxSIS : 2 diseases 3 zones at max /rhohat
    if(1) %if 0,1,2 (code a modifier en consequences)
        figure(4)
        delete(p_rho1);delete(p_rho);delete(p_rho13);delete(p_rho12)
        a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [78/255, 178/255, 55/255];a.FaceAlpha = 0.4;  
        a=area([cs3dis cs3dis cs2dis cs2dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [78/255, 178/255, 55/255];a.FaceAlpha = 0.15;  
        a=area([cs2dis cs2dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [78/255, 178/255, 55/255];a.FaceAlpha = 0.05;  
       if(alpha2>alpha3)
            legend('freq. de dep. vol. IST 1 ','freq. de dep. vol. IST 2 ',...
            'freq. de dep. vol. comb. IST 1 \& IST 2',...
            'elimination de IST 1 \& IST 2',...
            'seulement IST 2 est eliminee', 'IST 1 \& IST 2 persistent',...
            'Interpreter','latex')
        else
            legend('freq. de dep. vol. IST 1 ','freq. de dep. vol. IST 2 ',...
            'freq. de dep. vol. comb. IST 1 \& IST 2',...
            'elimination de IST 1 \& IST 2',...
            'seulement IST 1 est eliminee', 'IST 1 \& IST 2 persistent',...
            'Interpreter','latex')
        end
        
       
        figure(5) %R's
        a=area([vecC(1) vecC(1) cs3dis cs3dis],[1,limyR,limyR,1 ],'LineStyle','none'); a(1).FaceColor = [78/255, 178/255, 55/255];a.FaceAlpha = 0.4;  
        %uistack(a,'bottom') %met derriere
        a=area([cs3dis cs3dis cs2dis cs2dis],[1,limyR,limyR,1 ],'LineStyle','none'); a(1).FaceColor = [78/255, 178/255, 55/255];a.FaceAlpha = 0.15;  
        %uistack(a,'bottom')
        a=area([cs2dis cs2dis vecC(end) vecC(end)],[1,limyR,limyR,1],'LineStyle','none'); a(1).FaceColor = [78/255, 178/255, 55/255];a.FaceAlpha = 0.05;  
        %uistack(a,'bottom')
       if(alpha2>alpha3)
            legend('$R_{IST 1}$','$R_{IST 2}$','$R_{IST 1}$','$R_{IST 2}$',...
            'elimination de IST 1 \& IST 2',...
            'seulement IST 2 est eliminee', 'IST 1 \& IST 2 persistent',...
            'Interpreter','latex')
        else
            legend('$R_{IST 1}$','$R_{IST 2}$','$R_{IST 1}$','$R_{IST 2}$',...
            'elimination de IST 1 \& IST 2',...
            'seulement IST 1 est eliminee', 'IST 1 \& IST 2 persistent',...
            'Interpreter','latex')
        end

        figure(6) %prevalence
        a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limyP,limyP,0 ],'LineStyle','none'); a(1).FaceColor = [78/255, 178/255, 55/255];a.FaceAlpha = 0.4;  
        %uistack(a,'bottom') %met derriere
        a=area([cs3dis cs3dis cs2dis cs2dis],[0,limyP,limyP,0 ],'LineStyle','none'); a(1).FaceColor = [78/255, 178/255, 55/255];a.FaceAlpha = 0.15;  
        %uistack(a,'bottom')
        a=area([cs2dis cs2dis vecC(end) vecC(end)],[0,limyP,limyP,0 ],'LineStyle','none'); a(1).FaceColor = [78/255, 178/255, 55/255];a.FaceAlpha = 0.05;  
        %uistack(a,'bottom')        
        if(alpha2>alpha3)
            legend('Prevalence IST 1','Prevalence IST 2','Prevalence IST 1 ','Prevalence IST 2',...
            'elimination de IST 1 \& IST 2',...
            'seulement IST 2 est eliminee', 'IST 1 \& IST 2 persistent',...
            'Interpreter','latex')
        else
            legend('Prevalence IST 1','Prevalence IST 2','Prevalence IST 1 ','Prevalence IST 2',...
            'elimination de IST 1 \& IST 2',...
            'seulement IST 1 est elimine', 'IST 1 \& IST 2 persistent',...
            'Interpreter','latex')
        end
    end
    if(0) %0, 2 hiv, chlamydia
        delete(p_rho3);delete(p_rho);delete(p_rho13);delete(p_rho23)
        a=area([vecC(1) vecC(1) cs2dis cs2dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
        a=area([cs2dis cs2dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
        legend('$\hat\gamma_{hiv}$ (specifique VIH)','$\hat\gamma_{ch}$ (specifique chlamydia)','$\hat\gamma_{hiv\times ch}$ (combine HIV et chlamydia)',...
        'elimination of HIV \& chlamydia',...
        'HIV \& chlamydia persistent',...
        'Interpreter','latex')
    end
end %end plot SISxSIS

%% SIRxSIS² : ajout de 3 zones sur le graphe: epidemie when rho=0, 1 disease eliminated, 2 diseases eliminated, 3 diseases eliminated
% hold on
% 
% % 0, 1, 2, 3 diseases eliminated zones
% if(size(cs1dis,2)~=0 && size(cs2dis,2)~=0 && cs1dis>cs2dis && cs2dis>cs3dis)
%     a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
%     a=area([cs3dis cs3dis cs2dis cs2dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.3;  
%     a=area([cs2dis cs2dis cs1dis cs1dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.15;  
%     a=area([cs1dis cs1dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
%     legend('$\hat\gamma_{hiv}$','$\hat\gamma_{ch}$','$\hat\gamma_{s}$','$\hat\gamma_{hiv\times ch}$',...
%     '$\hat\gamma_{ch\times s}$','$\hat\gamma_{hiv \times s}$','$\hat\gamma$',...
%     'hiv,ch \& s eliminated', [num2str(jj),' \& ', num2str(kk), ' eliminated'],...
%     ['only ', num2str(kk), ' eliminated'],...
%     'no disease eliminated',...
%     'Interpreter','latex')
% 
% % 0, 3 diseases eliminated zones
% elseif(size(cs1dis,2)==0 && size(cs2dis,2)==0)
%     a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
%     a=area([cs3dis cs3dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
%     legend('$\hat\rho_1$','$\hat\rho_2$','$\hat\rho_3$','$\hat\rho_{1\times 2}$',...
%     '$\hat\rho_{2\times 3}$','$\hat\rho_{1 \times 3}$','$\hat\rho$',...
%     '1,2 \& 3 eliminated', 'no disease eliminated',...
%     'Interpreter','latex')
% 
% % 0, 3 diseases eliminated
% elseif(size(cs1dis,2)~=0 && size(cs2dis,2)==0 && cs3dis>=cs1dis)
%     a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
%     a=area([cs3dis cs3dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
%     legend('$\hat\rho_1$','$\hat\rho_2$','$\hat\rho_3$','$\hat\rho_{1\times 2}$',...
%     '$\hat\rho_{2\times 3}$','$\hat\rho_{1 \times 3}$','$\hat\rho$',...
%     '1,2 \& 3 eliminated', 'no disease eliminated',...
%     'Interpreter','latex')
% 
% % 0, 2, 3 diseases eliminated
% elseif(size(cs1dis,2)==0 && size(cs2dis,2)~=0 && cs3dis<cs2dis || size(cs1dis,2)~=0 && size(cs2dis,2)~=0 && cs1dis==cs2dis && cs3dis<cs2dis)
%     a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
%     a=area([cs3dis cs3dis cs2dis cs2dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.3;  
%     a=area([cs2dis cs2dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
%     legend('$\hat\rho_1$','$\hat\rho_2$','$\hat\rho_3$','$\hat\rho_{1\times 2}$',...
%     '$\hat\rho_{2\times 3}$','$\hat\rho_{1 \times 3}$','$\hat\rho$',...
%     '1,2 \& 3 eliminated', [num2str(jj),' \& ', num2str(kk), ' eliminated'], 'no disease eliminated',...
%     'Interpreter','latex')
% 
% % 0, 1, 3 diseases eliminated
% elseif(size(cs1dis,2)~=0 && size(cs2dis,2)==0 && cs3dis<cs1dis)
%     a=area([vecC(1) vecC(1) cs3dis cs3dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.5;  
%     a=area([cs3dis cs3dis cs1dis cs1dis],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.15;  
%     a=area([cs1dis cs1dis vecC(end) vecC(end)],[0,limy,limy,0 ],'LineStyle','none'); a(1).FaceColor = [132/255, 151/255, 176/255];a.FaceAlpha = 0.05;  
%     legend('$\hat\rho_1$','$\hat\rho_2$','$\hat\rho_3$','$\hat\rho_{1\times 2}$',...
%     '$\hat\rho_{2\times 3}$','$\hat\rho_{1 \times 3}$','$\hat\rho$',...
%     '1,2 \& 3 eliminated', ['only ', num2str(kk),' eliminated'], 'no disease eliminated',...
%     'Interpreter','latex')
% end

%%
