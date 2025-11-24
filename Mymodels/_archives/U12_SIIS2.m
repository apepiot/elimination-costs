function [U] = U12_SIIS2(b,beta1,beta2,gamma10,gamma20,eps1,eps2,sigma1,sigma2,nu1,nu2,mu,vecRho,c)
 %utility function of the SIISxSIIS model %numerically
    tspan = 0:1:500; U=zeros(1,length(vecRho)); i=1;
    for rho=vecRho
        Y0 = [99; 1; 1 ; 1 ; 1 ; 1 ; 1 ;1 ; 1];
        options = odeset('RelTol',1e-5,'Stats','off');%,'OutputFcn',@odeplot);
        [~,Ys] = ode45(@(t,Y) SIIS2(t,Y,b,beta1,beta2,gamma10,gamma20,rho,eps1,eps2,sigma1,sigma2,nu1,nu2,mu,'frequency'),tspan,Y0, options);
        N=b/mu;
        P = sum(Ys(end,2:9))/N;
        U(i) = rho*(P-c);
    end   
end