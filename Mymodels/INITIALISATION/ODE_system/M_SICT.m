function [M] = M_SICT(Lambda,theta0,sigma,p,zeta,eta,b,rho,mu)
	% dS/dt  = (1-p) \pi - \Lambda S-\mu S
	% dI/dt  = \Lambda S - (\sigma+\rho+\mu) I
	% dC/dt  = \sigma I - (\gamma(\rho)+\mu) C
	% dT/dt  = \rho I+\gamma(\rho) C + \eta I_p + (\gamma(0)+\eta) C_p - \mu T
    % with \Lambda = (\beta_I (I+I_p)+\beta_C (C+C_p))/N
    theta = theta0+rho;
        %S              %I              %C                 %T
    M = [-Lambda-mu,    0,              0,                 0;...
        Lambda,         -(sigma+rho+mu),0,                 0;...
        0,              sigma,          -(theta+mu),       0;...
        0,              rho,            theta,             -mu];
end