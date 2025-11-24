% Comparaison des resultats obtenus pour les modèles SIRxSIS targeting
% testing vs combined testing

%% Parameters
close all; 
%clear all;
%beta1=3.44;beta2=1.71;gamma1=1.32;gamma2=0.13;s1=1;s2=1;mu=1.224; b=0.4; %s1=0.67;s2=0.42
%beta1=0.42;beta2=2.18;gamma1=0.37;gamma2=1.71;s1=1;s2=1;mu=0.03; b=0.03; %s1=0.67;s2=0.42
%Paris
mu = 1/35;
PHIV = 16.1/100*9/100;
%PHIV=10/100;
gamma1 = 1/2.7;
beta1 = mu./(mu./(gamma1+mu)-PHIV);
R1 = beta1/(gamma1+mu);
s1=1;
Pch = 5.3/100;
R2 = 1/(1-Pch);
gamma2 = 1/1.5;
beta2=R2*(gamma2+mu);
s2=1;
b=10;


%%
R1     = beta1/(gamma1+mu);
R2     = beta2/(gamma2+mu);
alpha1 = beta1/s1*(1-1/R1);
alpha2 = beta2/s2*(1-1/R2);
minalpha = min(alpha1,alpha2);
maxalpha = max(alpha1,alpha2);

c11 = -mu*(1-1/R1)/beta1; %c such that argmax U1 = alpha1
c10 = mu*(R1-1)/beta1; %c such that argmax U1 = 0
c22 = 1/R2-1; %c such that argmax U2 = alpha2
c20 = -c22; %c such that argmax U2 = 0

%vectors c1,c2 for the map
lowboundC1  = c11-(c10-c11)/2;
highboundC1 = c10+(c10-c11)/2;
lowboundC2  = c22-(c20-c22)/2;
highboundC2 = c20+(c20-c22)/2;
maxhighbound = max(highboundC1,highboundC2);
minlowbound = min(lowboundC1,lowboundC2);
c1 = minlowbound:(c10-c11)/20:maxhighbound;
c2 = minlowbound:(c20-c22)/20:maxhighbound;

[C1,C2] = meshgrid(c1,c2);

%% Plot 3D : gradient du rhomax1 and rhomax2 - SIRxSIS targeting testing

%gradient colors kJ:yellow and kB:blue as a function of rhohat1 (yellow) and rhohat2
%(blue)
kJ = (1-beta1/R1/s1*(sqrt(R1.*mu./(mu+beta1.*C1))-1)./alpha1).*(C1>=c11 & C1<c10) + (C1>=c10); %jaune pour C1
kB = (C2-c22)/(c20-c22).*(C2<c20 & C2>=c22) + (C2>=c20); %hatrho2 est lineaire
 
col(:,:,1) = kB; % red
col(:,:,2) = kB; % green
col(:,:,3) = kJ; % blue

%plot
figure()
h1 = axes ;
image(c1,c2,col)
set(h1,'YDir','normal'); %to make y axe sorted
hold on

%hatched region
% [a] = fill([c11,c11,highboundC1,highboundC1,lowboundC1,lowboundC1],[highboundC2,c22,c22,lowboundC2,lowboundC2,highboundC2],'w') %color has no purpose here
% hPatch1 = findobj(a, 'Type', 'patch');
% hh1 = hatchfill(hPatch1, 'single', -45, 10); 
% set(a, 'LineStyle', 'none') %remove outline
% set(hh1, 'Color', 0.4*[1,1,1]) %set color of the hatched region

%Remove axes tic
set(gca,'YTick',[]); %which will get rid of all the markings for the y axis
set(gca,'XTick',[]); %which will get rid of all the markings for the x axis

%colorbar for C2 (blue)
colorbar
cb1 = 1-[1:-0.01:0]'*[1,1,0] ;
colormap(cb1);
cb11 = colorbar;
cb11.Location = 'westoutside';
cb11.Label.String = '$\hat\rho_2$';
cb11.Label.Interpreter = 'latex';
cb11.XTickLabel = [{'$\rho_2\prime$'},{''},{''},{''},{''},{''},{''},{''},{''},{''},{'$0$'}];
set(cb11,'TickLabelInterpreter','latex')

%newcolorbar (for C1, yellow)
hold on

cb22 = newcolorbar('southoutside');
%plot(0,0)
cb2 = 1-[1:-0.01:0]'*[0,0,1] ;
colormap(gca,cb2);
cb22.Label.String = '$\hat\rho_1$';
cb22.Label.Interpreter = 'latex';
cb22.XTickLabel = [{'$\rho_1\prime$'},{''},{''},{''},{''},{''},{''},{''},{''},{''},{'$0$'}];
set(cb22,'TickLabelInterpreter','latex')

hold on

plot([lowboundC1,highboundC1],[c20,c20], 'Color',0.85*[1,1,1],'LineStyle','--')%c20
plot([lowboundC1,highboundC1],[c22,c22], 'Color',0.85*[1,1,1],'LineStyle','--')%c22
plot([c11,c11],[lowboundC2,highboundC2], 'Color',0.85*[1,1,1],'LineStyle','--')%c11
plot([c10,c10],[lowboundC2,highboundC2], 'Color',0.85*[1,1,1],'LineStyle','--')%c10

%Zones
text((highboundC1+c10)/2,(highboundC2+c20)/2,'$(0,0)$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text((lowboundC1+c11)/2,(lowboundC2+c22)/2,'$(\rho_1\prime,\rho_2\prime)$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','Color','w')

text((highboundC1+c10)/2,(lowboundC2+c22)/2,'$(0,\rho_2\prime)$', 'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'Interpreter', 'latex','Color','w')
text((lowboundC1+c11)/2, (highboundC2+c20)/2,'$(\rho_1\prime,0)$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')

%labels sur les axes
ex = (highboundC1-lowboundC1)/25;
ey = (highboundC2-lowboundC2)/25;
text(c11, lowboundC2-ey,'$c_{1}^1$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(c10, lowboundC2-ey,'$c_{0}^1$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(lowboundC1-ex, c22,'$c_{2}^2$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(lowboundC1-ex, c20,'$c_{0}^2$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')

xlim([lowboundC1,highboundC1])
ylim([lowboundC2,highboundC2])
axis off

title([{'$\hat\rho(c)$'},...
            {['$\beta_1$=',num2str(beta1), ' $\beta_2$=',num2str(beta2),' $\gamma_1(0)$=',num2str(round(gamma1,2)),' $\gamma_2(0)$=',...
            num2str(round(gamma2,2)), ' $s_1$=', num2str(s1),' $s_2$=', num2str(s2), ' $\mu$=', num2str(mu), ' $\pi$=',  num2str(b)]},...
            {[' $\mathtt R_1(0)$=' num2str(round(R1,2)), ' $\mathtt R_2(0)=$' num2str(round(R2,2)),...
            ' $\rho_1\prime$=' num2str(round(alpha1,2)),' $\rho_2\prime$=' num2str(round(alpha2,2)),...
             ' $c_0^1$=', num2str(round(c10,2)),' $c_1^1$=', num2str(round(c11,2)),...
             ' $c_0^2$=', num2str(round(c20,2)),' $c_2^2$=', num2str(round(c22,2))...
            ]}],...
            'Interpreter','latex')
hold off
%% Plot 3D : gradient du rhomax - SIRxSIS combined testing
%C = max(C1,C2);
%C=C1+C2;
%C = min(C1,C2)
C=C1/2+C2/2;
%cost thresholds for U12
c121 = (beta1*((gamma1 + mu)/beta1 - 1)*(s2/beta2 - (mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)*((s2 - s1*s2)/beta1 + ((s2 - (beta1*mu*s1)/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))^2)*(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1))) + (s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1))) + (mu*s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/((beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1))^2*(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)))))/(beta2 + gamma1 - beta1*((gamma1 + mu)/beta1 - 1) + beta1*s2*((gamma1 + mu)/beta1 - 1)) + (mu*(s1 - s1*s2)*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)*((gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta1 + ((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)))))/(beta2 + gamma1 - beta1*((gamma1 + mu)/beta1 - 1) + beta1*s2*((gamma1 + mu)/beta1 - 1))^2 + (beta1*mu*s1*((gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta1 + ((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)))))/((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))^2*(beta2 + gamma1 - beta1*((gamma1 + mu)/beta1 - 1) + beta1*s2*((gamma1 + mu)/beta1 - 1)))))/s1 - (gamma2 + mu - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta2 + (mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)*((gamma2 + beta1*s2*((gamma1 + mu)/beta1 - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1)/beta1 + ((gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1))*(gamma2 + mu + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1) - (beta1*s2*((gamma1 + mu)/beta1 - 1))/s1))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - beta1*((gamma1 + mu)/beta1 - 1)) - 1)))))/(beta2 + gamma1 - beta1*((gamma1 + mu)/beta1 - 1) + beta1*s2*((gamma1 + mu)/beta1 - 1)) + 1;
c122 = (mu*((gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1))/beta1 + ((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))))*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))/(beta2 + gamma1 + beta2*s1*((gamma2 + mu)/beta2 - 1) - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - (gamma2 + mu - beta2*((gamma2 + mu)/beta2 - 1))/beta2 + (beta2*((gamma2 + mu)/beta2 - 1)*(s2/beta2 - (mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1)*((s2 - s1*s2)/beta1 + ((s2 - (beta1*mu*s1)/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)^2)*(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))) + (s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))) + (mu*s1*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/((beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))^2*(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2))))/(beta2 + gamma1 + beta2*s1*((gamma2 + mu)/beta2 - 1) - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) + (mu*((gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1))/beta1 + ((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))))*(s1 - s1*s2)*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1))/(beta2 + gamma1 + beta2*s1*((gamma2 + mu)/beta2 - 1) - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)^2 + (beta1*mu*s1*((gamma2 - beta2*((gamma2 + mu)/beta2 - 1) + beta2*s1*((gamma2 + mu)/beta2 - 1))/beta1 + ((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)*(gamma2 + mu + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1) - beta2*((gamma2 + mu)/beta2 - 1)))/(beta1*(beta2 + mu*(beta1/(gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2) - 1)))))/((gamma1 + mu - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2)^2*(beta2 + gamma1 + beta2*s1*((gamma2 + mu)/beta2 - 1) - (beta2*s1*((gamma2 + mu)/beta2 - 1))/s2))))/s2 + 1;
P20  = 1-1./R2; lambda10 = mu*(R1-1); i1ES120 = lambda10./(beta2+gamma1)*(gamma2/beta1+(lambda10+gamma2+mu)/R1/(lambda10+beta2)); P120 = P20 + i1ES120;
c120 = P120;

%cost thresholds (c11,c10,c22,c20) for U1 and U2 :already defined in preamble

rhomax1 = (beta1./(R1*s1).*(sqrt(R1*mu./(mu+beta1*C))-1)).*(C>c11 & C<c10) + alpha1*(C<=c11);
rhomax2 = beta2./(2*s2).*(1-1./R2-C).*(C>c22 & C<c20) + alpha2*(C<=c22);
if(alpha1<alpha2)
    %then argmaxU = rhomax12 or rhomax2 when c>c22
    MAT = (C>c121 & C<c120);
    
    %methode bourrin où on recalcule le max pour chaque C (alors qu'avec la
    %fonction max, il y a bcp de doublons)
    RHOMAX12_2 = double(MAT); %dans cette matrice, on stocke le max entre rho2 et rho12
    [ni,nj]    = size(C); 

    for i=1:ni
        for j=1:nj
            if (MAT(i,j)~=0) %i.e. C est trop grand ou trop petit, rhohat=0 ou alpha1
                cij = C(i,j);
                rhomax2ij = rhomax2(i,j);
                U2max  = rhomax2ij*(1-(gamma2+mu+s2*rhomax2ij)/beta2-cij); %maxU1, i.e. U1 en hatrho1
                fun12 = @(rho) U12_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,cij); %-U12
                rhomax12ij = fmincon(fun12,minalpha);
                U12max = -U12_SIRSIS7(rhomax12ij,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,cij);
                if(U2max<U12max && rhomax2ij~=alpha2)
                    RHOMAX12_2(i,j) = rhomax12ij;
                else
                    RHOMAX12_2(i,j) = rhomax2ij;
                end
            end 
        end
    end
    RHOMAX = 0.*(C>=c120) + RHOMAX12_2.*(C>c121 & C<c120) + rhomax2.*(C<=c121);
else %if alpha2<=alpha1
    %then argmaxU = rhomax12 or rhomax1 when c>c11
    MAT = (C>c122 & C<c120);

    RHOMAX12_1 = double(MAT); %dans cette matrice, on stocke le max entre rho2 et rho12
    [ni,nj]    = size(C); 

    for i=1:ni
        for j=1:nj
            if (MAT(i,j)~=0) %i.e. C est trop grand ou trop petit, rhohat=0 ou alpha1
                cij = C(i,j);
                rhomax1ij = rhomax1(i,j);
                U1max  = rhomax1ij*(mu/beta1*(beta1/(gamma1+s1*rhomax1ij+mu)-1)-cij); %maxU1, i.e. U1 en hatrho1
                fun12 = @(rho) U12_SIRSIS7(rho,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,cij); %-U12
                rhomax12ij = fmincon(fun12,minalpha);
                U12max = -U12_SIRSIS7(rhomax12ij,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,cij);
                if(U1max<U12max && rhomax1ij~=alpha1)
                    RHOMAX12_1(i,j) = rhomax12ij;
                else
                    RHOMAX12_1(i,j) = rhomax1ij;
                end
            end 
        end
    end
    RHOMAX = 0.*(C>=c120) + RHOMAX12_1.*(C>c122 & C<c120) + rhomax1.*(C<=c122);
end
%% green plot
% %facteur de reduction de couleur
% epsilon = 1;
% kG = epsilon*RHOMAX./maxalpha;%smooth color
% col(:,:,1) = 0; % red
% col(:,:,2) = (1-kG); % green %mettre à 1 pour du vert
% col(:,:,3) = 0; % blue
% 
% figure
% h1 = axes ;
% image(c1,c2,col)
% set(h1,'YDir','normal'); %to make y axe sorted
% hold on
% 
% %colorbar
% colorbar
% %cb1 = 1-[0.6:-0.01:0]'*[1,0,1] ; %green
% cb1 = 1-[1*epsilon:-0.01:0]'*[1,1,1] ;
% colormap(cb1);
% cb11 = colorbar;
% cb11.Location = 'westoutside';
% cb11.Label.String = '$\hat\rho_2$';
% cb11.Label.Interpreter = 'latex';
% cb11.XTickLabel = [{'$\rho_2\prime$'},{''},{''},{''},{''},{''},{''},{''},{''},{''},{'$0$'}];
% set(cb11,'TickLabelInterpreter','latex')
% 
% plot([minlowbound,maxhighbound],[minlowbound,maxhighbound],'k-')

%% manually - colorful plot

%creating colorbar map, such that there is gradient of blue (disease 2) and
%yellow (disease 1);
kJ = subdivisedColormap([[1.,1.,1.];[1,1,0]],round(min(alpha1/alpha2,1)*100,0),'lin'); %jaune pour C1
kB = subdivisedColormap([[1.,1.,1.];[0,0,1]],round(min(alpha2/alpha1,1)*100,0),'lin'); 

if(alpha1<alpha2)
    KJ =[kJ;ones(1,size(kB,1)-size(kJ,1))'*[1,1,0]];
    mapJB = 1-((1-kB) + (1-KJ));
end

figure()
s = surf(C1,C2,RHOMAX);
%s=contourf(C1,C2,RHOMAX)
%map1 = subdivisedColormap([[1,1,1];[0.,0.5,0.];[0.,0.,0.]],10); %green colorbar

colormap(mapJB)
colorbar('Direction','reverse')
%colormap(flipud(parula))
s.EdgeColor = 'none';
s.FaceColor = 'interp';
cb2 = colorbar;

cb2.Label.Interpreter = 'latex';
cb2.Label.String = '$\hat\rho$';
set(cb2,'YTick',[0,alpha1,alpha2],'YTickLabel',{'0','$\rho_1\prime$','$\rho_2\prime$'},'TickLabelInterpreter','latex')

hold on;

%plot3([minlowbound,maxhighbound],[minlowbound,maxhighbound],[maxalpha,maxalpha],'k-')

%labels sur les axes
ex = (maxhighbound-minlowbound)/25;
ey = (maxhighbound-minlowbound)/25;
text(c11, minlowbound-ey,'$c_{1}^1$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(c10, minlowbound-ey,'$c_{0}^1$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(0, minlowbound-ey,'$0$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(c120, minlowbound-ey,'$c_0^{12}$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(c121, minlowbound-ey,'$c_1^{12}$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(c122, minlowbound-ey,'$c_2^{12}$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')

text(minlowbound-ex, c22,'$c_{2}^2$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(minlowbound-ex, c20,'$c_{0}^2$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(minlowbound-ex, 0,'$0$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(minlowbound-ex, c120,'$c_0^{12}$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(minlowbound-ex, c121,'$c_1^{12}$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')
text(minlowbound-ex, c122,'$c_2^{12}$', 'HorizontalAlignment', 'center', 'Interpreter', 'latex')


xlim([minlowbound,maxhighbound])
ylim([minlowbound,maxhighbound])
axis off

title([{'$\hat\rho(c)$'},...
            {['$\beta_1$=',num2str(round(beta1,2)), ' $\beta_2$=',num2str(round(beta2,2)),' $\gamma_1(0)$=',num2str(round(gamma1,2)),' $\gamma_2(0)$=',...
            num2str(round(gamma2,2)), ' $s_1$=', num2str(s1),' $s_2$=', num2str(s2), ' $\mu$=', num2str(round(mu,2)), ' $\pi$=',  num2str(round(b,2)),...
            ' $\mathtt R_1(0)$=' num2str(round(R1,2)), ' $\mathtt R_2(0)=$' num2str(round(R2,2))]},...
            {[' $\rho_1\prime$=' num2str(round(alpha1,2)),' $\rho_2\prime$=' num2str(round(alpha2,2)),...
            ' $c_0^1$=', num2str(round(c10,2)),' $c_1^1$=', num2str(round(c11,2)),...
             ' $c_0^2$=', num2str(round(c20,2)),' $c_2^2$=', num2str(round(c22,2)),...
             ' $c_0^{12}$=', num2str(round(c120,2)),' $c_1^{12}$=', num2str(round(c121,2)),' $c_2^{12}$=', num2str(round(c122,2)),' $c_s$=', num2str(round(cswitch,2))...
            ]}],...
            'Interpreter','latex')


%text((minlowbound+c11)/2,(minlowbound+c11)/2, alpha1, '$\rho_1\prime$','FontSize',20,'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'Color','w')
text((maxhighbound+c120)/2,(maxhighbound+c120)/2, maxalpha, '$0$', 'FontSize',20,'HorizontalAlignment', 'center', 'Interpreter', 'latex')

%% for c=max(c1,c2) (case alpha1<alpha2)
%let's find the cost c such that argmaxU = rho2'
ic = sum(RHOMAX(1,:)==alpha2);
calpha2 = c1(ic);
%let's find the cost c such that argmaxU = rho1'
ic1 = sum(RHOMAX(1,:)>alpha1);
calpha1 = c1(ic1);

%Zones 
text((minlowbound+calpha2)/2,(minlowbound+calpha2)/2, alpha2, '$\rho_2\prime$','FontSize',20,'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'Color','w')
text((0+c22)/1.4,(0+c22)/1.4, alpha2, '$\hat\rho_2$','FontSize',18,'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'Color','w')
text((0+c120)/2,(0+c120)/2, alpha2, '$\hat\rho_{12}$','FontSize',18,'HorizontalAlignment', 'center', 'Interpreter', 'latex')

%zone 0
plot3(c120*ones(1,100), linspace(minlowbound,c120,100), zeros(1,100)+maxalpha*0.01,'--k')
plot3(linspace(minlowbound,c120,100), c120*ones(1,100),zeros(1,100)+maxalpha*0.01,'--k')

%zones maxalpha
c12min = min(c121,c122);
%plot3(c12min*ones(1,100), linspace(minlowbound,c12min,100), maxalpha*1.01*ones(1,100),':w')
%plot3(linspace(minlowbound,c12min,100), c12min*ones(1,100), maxalpha*1.01*ones(1,100),':w')

plot3(calpha2*ones(1,100), linspace(minlowbound,calpha2,100), maxalpha*1.01*ones(1,100),'--w')
plot3(linspace(minlowbound,calpha2,100), calpha2*ones(1,100), maxalpha*1.01*ones(1,100),'--w')

%line where rhohat=rho1'
plot3(linspace(minlowbound,calpha1,100), calpha1*ones(1,100), maxalpha*1.01*ones(1,100),'-k')
plot3(calpha1*ones(1,100), linspace(minlowbound,calpha1,100), maxalpha*1.01*ones(1,100),'-k')


%% for c=c1+c2 (case alpha1<alpha2)
%let's find the cost c such that argmaxU = rho2'
ic = sum(RHOMAX(1,:)==alpha2);
calpha2 = c1(ic)+c2(1);
icbis = sum(RHOMAX(:,1)==alpha2);
text((minlowbound+calpha2)/2,(minlowbound+calpha2)/2, alpha2, '$\rho_2\prime$','FontSize',20,'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'Color','w')
plot3([c1(ic),c2(1)],[c2(1),c2(icbis)],[alpha2,alpha2],'--w')

%zone 0
plot3([c120-c2(end),c1(end)],[c2(end),c120-c1(end)],[alpha2*0.01,alpha2*0.01],'--k')

%argmaxU=rho1'
plot3([c121-c2(end),c1(end)],[c2(end),c121-c1(end)],[alpha1,alpha1],'-k')


%% autre representations
s1=surf(C1(1:28,1:23),C2(1:28,1:23),RHOMAX(1:28,1:23));
s2=surf(C1(29:60,24:60),C2(29:60,24:60),RHOMAX(29:60,24:60));
s3=surf(C1(1:28,24:60),C2(1:28,24:60),RHOMAX(1:28,24:60));
s4=surf(C1(29:60,1:23),C2(29:60,1:23),RHOMAX(29:60,1:23));

RHOMAX2=RHOMAX;
RHOMAX2(C>c122 & C<c120)=NaN;
figure
s5 = surf(C1,C2,RHOMAX2);
