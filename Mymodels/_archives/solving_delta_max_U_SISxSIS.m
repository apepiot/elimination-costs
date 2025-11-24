%% Solving dU/ddelta pour le modele both voluntary testing the same of the SISxSIS model

clear all;

%% symbolic

syms gamma1 gamma2 beta1 beta2 mu delta s1 s2 s

s1 = s;
s2 = s ;
gamma1p = gamma1+s1*delta;
gamma2p = gamma2+s2*delta;


R01 = beta1/(gamma1p + mu) ;
R02 = beta2/(gamma2p + mu) ;

beta_hat = beta1+beta2 -mu;

i = (1 - (gamma1p/R02 + gamma2p/R01 + mu)/beta_hat);  %I/N

dU = i - delta/beta_hat* (s1*(gamma2p+mu)/beta2 + gamma1p*s2/beta2 + s2*(gamma1p+mu)/beta1 + gamma2p*s1/beta1);

sol = solve(dU == 0, delta);

simplify(sol);
%%% du/dt=0

sol1 = sol(1);

sol1*3*s1*s2*(beta1+beta2);
    
A = 3*(1/beta1+ 1/beta2);
%B = 2*(gamma1*s2+gamma2*s1)*(1/beta1+1/beta2) + 2*s1*mu/beta2 + 2*s2*mu/beta1;
B = 2*(gamma1 + gamma2 + mu)*(1/beta2 + 1/beta1);
%C = mu - beta_hat + gamma1*gamma2/beta1 + gamma2*mu/beta1 + gamma1*gamma2/beta2 + gamma1*mu/beta2;
C = -beta1 - beta2 + 2*mu + gamma1*(gamma2+mu)/beta2 + gamma2*(gamma1+mu)/beta1;
DELTA = B^2 - 4*A*C;

simplify(DELTA) 

syms x
X = solve(A*x^2 + B*x + C == 0, x)

simplify(X(1))

%%evaluer U en sol1 ou sol2
d = X(1);
gamma1p = gamma1+s1*d;
gamma2p = gamma2+s2*d;
R01 = beta1/(gamma1p + mu) ;
R02 = beta2/(gamma2p + mu) ;


U = 1 - (gamma1p/R02 + gamma2p/R01 + mu)/beta_hat;

simplify(U)



