function [U] = U_SICR_v4(param, mu, b, rho, c, f)
    [Ui,~] = U1_SICR_v4(param, mu, b, rho, c, f);
    alpha = param.alpha;
    U = Ui.*(rho<alpha)+ (-alpha*c).*(rho>=alpha);
end