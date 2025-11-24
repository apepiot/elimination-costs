function [U] = U_SICTP(param,mu,b,rho,c,f)  
    alpha=param.alpha;
    [Ui,~] = U1_SICTP(param,mu,b,rho,c,f) ;    
    U = Ui.*(rho<=alpha) + (-c*alpha)*(rho>alpha);
end