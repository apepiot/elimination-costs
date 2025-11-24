[beta1,beta2,gamma1,gamma2,s1,s2,b,mu] = random_parameters('true', 'true')

s1=1;s2=1;mu=0;
epsilon1 = 0:0.005:0.2;
epsilon2 = 0.01;

beta

a = gamma1/beta2 + gamma2/beta1 + 1;
c = beta1+beta2;
b = gamma1*gamma2*(1/beta1+1/beta2);
rho_hat = c*(1-sqrt((1+a)*(a+gamma1*gamma2/(beta1*beta2)))/(1+a))





