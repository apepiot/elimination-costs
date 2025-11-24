syms mu Lambda1 Lambda2 Lambda3 Lambda4 nu1 nu2 nu3 nu4 eps1 eps2 eps3 eps4 sigma1 sigma2 sigma3 sigma4 gamma1 gamma2 gamma3 gamma4 b rho;
syms X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16

eqn=[0 == b - X1*(Lambda1 + Lambda2 + mu) + X2*rho + X5*rho + X6*rho + X7*rho + X8*rho + X10*rho + X11*rho + X12*rho + X14*rho + X15*rho + X4*(gamma1 + nu1 + rho) + X13*(gamma2 + nu2 + rho) + X3*(nu1 + rho) + X9*(nu2 + rho);
0 == Lambda1*X1 + X10*nu2 - X2*(Lambda2 + mu + rho + sigma1) + X14*(gamma2 + nu2);
0 == X11*nu2 - X3*(Lambda2 + mu + nu1 + rho) + X15*(gamma2 + nu2) - X2*sigma1*(eps1 - 1);
0 == X12*nu2 - X4*(Lambda2 + gamma1 + mu + nu1 + rho) + X16*(gamma2 + nu2) + X2*eps1*sigma1;
0 == Lambda2*X1 + X7*nu1 - X5*(Lambda1 + mu + rho + sigma2) + X8*(gamma1 + nu1);
0 == Lambda2*X2 + Lambda1*X5 - X6*(mu + rho + sigma1 + sigma2);
0 == Lambda2*X3 - X7*(mu + nu1 + rho + sigma2) - X6*sigma1*(eps1 - 1);
0 == Lambda2*X4 - X8*(gamma1 + mu + nu1 + rho + sigma2) + X6*eps1*sigma1;
0 == X11*nu1 - X9*(Lambda1 + mu + nu2 + rho) + X12*(gamma1 + nu1) - X5*sigma2*(eps2 - 1);
0 == Lambda1*X9 - X10*(mu + nu2 + rho + sigma1) - X6*sigma2*(eps2 - 1);
0 == - X11*(mu + nu1 + nu2 + rho) - X7*sigma2*(eps2 - 1) - X10*sigma1*(eps1 - 1);
0 == X10*eps1*sigma1 - X12*(gamma1 + mu + nu1 + nu2 + rho) - X8*sigma2*(eps2 - 1);
0 == X15*nu1 - X13*(Lambda1 + gamma2 + mu + nu2 + rho) + X16*(gamma1 + nu1) + X5*eps2*sigma2;
0 == Lambda1*X13 - X14*(gamma2 + mu + nu2 + rho + sigma1) + X6*eps2*sigma2;
0 == X7*eps2*sigma2 - X15*(gamma2 + mu + nu1 + nu2 + rho) - X14*sigma1*(eps1 - 1);
0 == X8*eps2*sigma2 - X16*(gamma1 + gamma2 + mu + nu1 + nu2) + X14*eps1*sigma1]

solve(eqn, [X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16]);
 