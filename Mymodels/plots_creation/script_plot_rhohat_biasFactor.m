% script to compare different values of the biasFactor

vecC=linspace(-0.1,0.3,100);
biasFactor=1;
[tab,tabco,tabcn,tabTimes] = findRhohat_v4(2,paramTab,mu,b,vecC,biasFactor);
plot(vecC,tab.one(1).rhohat,'b-','DisplayName','Ct')
hold on
plot(vecC,tab.one(2).rhohat,'r-','DisplayName','HIV')
plot(vecC,tab.two.rhohat,'k-','DisplayName','b=1')
hold on
linS = {'--','-.',':'};i=1;
for biasFactor=[2,5,10]
    [tab,tabco,tabcn,tabTimes] = findRhohat_v4(2,paramTab,mu,b,vecC,biasFactor);
    plot(vecC,tab.two.rhohat,'Color','k','linestyle',linS{i},'DisplayName',['b=',num2str(biasFactor)])
    i=i+1;
end
legend
ylim([0,1.3*max(paramTab{1}.alpha,paramTab{2}.alpha)])
