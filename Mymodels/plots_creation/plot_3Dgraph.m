%this code plots different 3D graphs with R1(0) in x-axis, R2(0) in y-axis,
%and color grading in function of how one disease *drives* the other through
%CVT

yR2 = (gamma1+mu)/(gamma2+mu)*(R10-1)+1;
close all;

titre = ['IST1 : ', num2str(round(1/gamma1,2)),' years and IST 2 : ', num2str(round(1/gamma2,2)),' years'];
fname = ['C:\Users\Moi\Documents\IPLESP\These\Graphes\strategies_comparison\',num2str(type)]; %PC

if(1)
n=10; %nb de subdvisions colormap

figure(4)
%plot :  gradient of change
%colorbars
greyred = subdivisedColormap([[0.,0.,0.];[1,1,1];[0.85,0.325,0.098]],n, 'quad'); %2^n+1
greyyellow = subdivisedColormap([[0.,0.,0.];[1,1,1];[0.929,0.8,0.125]],n, 'quad'); %2^n+1

%put 0 as a median of the colorbar
amp = maxv-minv; %amplitude des valeurs
%maxv/amplitude %pourcentage de valeurs +
%minv/amplitude %pourcentage de valeurs -
minv1 = min(min(dis2drives1));maxv1 = max(max(dis2drives1));
minv2 = min(min(dis1drives2));maxv2 = max(max(dis1drives2));
    
if(abs(minv1)<maxv1)
   %greyred2 = greyred(ceil((2^n+1)*0.5*(1-abs(minv)/(maxv))):end,:);
    greyred2 = greyred(round((2^n+1)*(minv1 - (-maxv1)) / (2 * maxv1)):end, :);
else
    greyred2 = greyred(1:round((2^n+1)*abs(maxv1 - (-minv1)) / abs(2 * minv1)),:); %ok ?
end
if(abs(minv2)<maxv2)
    %greyyellow2 = greyyellow(ceil((2^n+1)*0.5*(1-abs(minv)/(maxv))):end,:);
    greyyellow2 = greyyellow(round((2^n+1)*(minv2 - (-maxv2)) / (2 * maxv2)):end, :);
else
    %greyred2 = greyred(1:ceil((2^n+1)*0.5*(1+abs(maxv)/abs(minv))),:);
    %greyyellow2 = greyyellow(1:ceil((2^n+1)*0.5*(1+abs(maxv)/abs(minv))),:);    
    greyyellow2 = greyyellow(1:round(abs(2^n+1)*(maxv2 - (-minv2)) / abs(2 * minv2)),:); %ok ?
end

%surf1=surf(R1,R2,tot)
surf1 = surf(R1,R2,dis2drives1, 'FaceColor','interp');
surf1.EdgeColor = 'none';
colormap(greyred2)
cb = colorbar;
cb.YTick = [-2 -1 0 1 2 3 4 5 6];
cb.YTickLabel = {'-200%','-100%', 'No Change', '+100%', '+200%','+300%','+400%','+500%','+600%'};
%cb.Position = cb.Position + 1e-10;
freezeColors %adding to the second area where aplha1>alpha2
title([{'Is combined testing better than specific testing ?'},...
    {'$(c_{j1\times2}-c_{j}^j)/|c_{j}^j|$ such that $\rho_j\prime<\rho_i\prime$'},{titre}], 'Interpreter','latex')
xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')
hold on
view(2)
%fname = '/Users/amandine/Desktop/These/Graphes/strategies_comparison/SIS2'; %MAC
set(gcf,'Position',[10 10 550 500])

if saveFiles
filename = ['strat_comp_',num2str(type),'_gamma1_',replace(num2str(round(gamma1,2)),'.',','),...
    '_gamma2_',replace(num2str(round(gamma2,2)),'.',','),'-1'];
saveas(gca, fullfile(fname, filename),'jpeg');
end
%%
surf2=surf(R1,R2,dis1drives2,'FaceColor','interp')
surf2.EdgeColor = 'none';
cb2 = colorbar('southoutside');
colormap(greyyellow2)
cb2.YTick = [-2 -1 0 1 2];
cb2.YTickLabel = {'-200%','-100%', 'No Change','+100%','+200%',};
plot3(R10,yR2,maxv*ones(length(R10),1),'-k','LineWidth',1.5);
ylim([1,max(R20)]);
view(2)
set(gcf,'Position',[10 10 550 500])

if saveFiles
filename = ['strat_comp_',num2str(type),'_gamma1_',replace(num2str(round(gamma1,2)),'.',',')...
    ,'_gamma2_',replace(num2str(round(gamma2,2)),'.',','),'-2'];
saveas(gca, fullfile(fname, filename),'jpeg');
end
end
%%
figure(5)
tot2 = (ALPHA1<ALPHA2 & (C11<CEL1))*1 +...
    (ALPHA1<ALPHA2 & (C11>=CEL1))*2 +...
    (ALPHA1>=ALPHA2 & (C22<CEL1))*4 +...
    (ALPHA1>=ALPHA2 & (C22>CEL1))*3;
surf4=surf(R1,R2,tot2);
surf4.EdgeColor = 'none';
new_cb2 = [[0.85,0.325,0.098];[0.8,0.8,0.8];[0.8,0.8,0.8];[0.9290,0.6940,0.125]]; %2^n+1
colormap(new_cb2)
%text(xG1,yG1,4,'\rho\prime_1<\rho\prime_2','HorizontalAlignment','center','FontSize',20)
%text(xG2,yG2,4,'\rho\prime_1>\rho\prime_2','HorizontalAlignment','center','FontSize',20)
hold on;
plot3(R10,yR2,4*ones(length(R10),1),'-k','LineWidth',1.5);
ylim([1,max(R20)]);
title([{'Is combined testing better than specific testing ?'},{titre}], 'Interpreter','latex')
xlabel('$\mathtt R_1(0)$', 'Interpreter','latex')
ylabel('$\mathtt R_2(0)$', 'Interpreter','latex')
view(2)

if saveFiles
filename = ['strat_comp_',num2str(type),'_gamma1_',replace(num2str(round(gamma1,2)),'.',','),...
    '_gamma2_',replace(num2str(round(gamma2,2)),'.',','),'-0'];
saveas(gca, fullfile(fname, filename),'png');
close all;
end
