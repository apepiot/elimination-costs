function [U] = U_SEIIS_v4(param,mu,b,rho,c,f)
    alpha=param.alpha;
    [Ui,~,~] = U1_SEIISv4(param,mu,b,rho,c,f);
    U = Ui.*(rho<alpha) + (-alpha*c).*(rho>=alpha);
end