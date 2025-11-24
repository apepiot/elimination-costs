function plot_4dis(vecC,tab)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%4 diseases
plot(vecC,tab.four.rhohat,'-','DisplayName','CtxNgxHIVxsyph.','LineWidth',3)
hold on
%3 diseases
plot(vecC,tab.three(1).rhohat,'--','DisplayName','CtxNgxHIV','LineWidth',2.5)
plot(vecC,tab.three(2).rhohat,'--','DisplayName','CtxNgxsyph.','LineWidth',2.5)
plot(vecC,tab.three(3).rhohat,'--','DisplayName','CtxHIVxsyph.','LineWidth',2.5)
plot(vecC,tab.three(4).rhohat,'--','DisplayName','NgxHIVxsyph.','LineWidth',2.5)

%2 diseases
plot(vecC,tab.two(1).rhohat,'-.','DisplayName','CtxNg','LineWidth',2)
plot(vecC,tab.two(2).rhohat,'-.','DisplayName','CtxHIV','LineWidth',2)
plot(vecC,tab.two(3).rhohat,'-.','DisplayName','Ctxsyph.','LineWidth',2)
plot(vecC,tab.two(4).rhohat,'-.','DisplayName','NgxHIV','LineWidth',2)
plot(vecC,tab.two(5).rhohat,'-.','DisplayName','Ngxsyph.','LineWidth',2)
plot(vecC,tab.two(6).rhohat,'-.','DisplayName','HIVxsyph.','LineWidth',2)

%1 disease
plot(vecC,tab.one(1).rhohat,':','DisplayName','Ct','LineWidth',1.5)
hold on
plot(vecC,tab.one(2).rhohat,':','DisplayName','Ng','LineWidth',1.5)
plot(vecC,tab.one(3).rhohat,':','DisplayName','HIV','LineWidth',1.5)
plot(vecC,tab.one(4).rhohat,':','DisplayName','syph.','LineWidth',1.5)

legend()
xlim([vecC(1),vecC(end)])
end

