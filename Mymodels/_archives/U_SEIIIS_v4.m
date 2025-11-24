function [U] = U_SEIIIS_v4(param,mu,b,rho,c,f)  
    alpha=param.alpha;
    [Ui,~] = U1_SEIIIS_v4(param,mu,b,rho,c,f) ;    
    U = Ui.*(rho<=alpha) + (-c*alpha)*(rho>alpha);
end