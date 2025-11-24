function [condTOT,cond1,cond2,cond3,cond4,cond5] = verifConditionsUmax(beta1,beta2,gamma1,gamma2,s1,s2,b,mu)
%fonction qui comparee les umax candidates entre eux
%   a verifier 
    R1 = beta1/(gamma1+mu);
    R2 = beta2/(gamma2+mu);
    alpha1 = beta1/s1*(1-1/R1);%rho1'
    alpha2 = beta2/s2*(1-1/R2);%rho2'
    rho10 = beta1*beta2/2/gamma2 - gamma1/2; %rho1**
    rho20 = beta1*beta2/2/gamma1 - gamma2/2; %rho2**

    Ucase1 = (alpha1+alpha2)^2/4/beta2;
    Ucase2 = (alpha1+alpha2)^2/(4*beta1);
    Ucase3 = (beta1*beta2 - gamma1.*gamma2).^2./(4*beta1*beta2*gamma1);
    Ucase4 = (beta1*beta2 - gamma1.*gamma2).^2./(4*beta1*beta2*gamma2);


    %x1 = beta2*(-alpha1/(sqrt(beta2*gamma1))+1);
    %x2 = beta2*(alpha1/(sqrt(beta2*gamma1))+1);
    %x3 = beta1*(-alpha2/(sqrt(beta1*gamma2))+1);
    %x4 = beta1*(alpha2/(sqrt(beta1*gamma2))+1);
    x5 = beta1*((gamma1-beta2)/sqrt(gamma1*beta1)+1);
    x6 = beta1*((-gamma1+beta2)/sqrt(gamma1*beta1)+1);
    x11 = beta1 + alpha2*(1-2*beta2/beta1*(1-sqrt(1-beta1/beta2)));

%     %conditions entre u1 et u2
%     cond11  = (Ucase1>Ucase2) & (beta2<beta1);
%     cond12  = (Ucase1<=Ucase2) & (beta2>=beta1);
%     cond1   = cond11 | cond12;
% 
%     %conditions entre u3 et u4
%     cond21  = (Ucase3>Ucase4) & (gamma1<gamma2);
%     cond22  = (Ucase3<=Ucase4) & (gamma1>=gamma2);
%     cond2   = cond21 | cond22;
% 
%     %conditions entre u3 et u2
%     cond311 = (Ucase3>Ucase2) & (beta2>gamma1) & gamma2>max(x1,0) & gamma2<x2;
%     cond312 = (Ucase3>Ucase2) & (beta2<gamma1) & (gamma2<=max(x1,0) | gamma2>x2);
%     cond313 = (Ucase3>Ucase2) & (beta2==gamma1) & (gamma2<x1) & x1==x2;
%     cond31  =  cond311 | cond312 | cond313;
% 
%     cond321 = (Ucase3<Ucase2) & (beta2>gamma1) & (gamma2<max(x1,0) | gamma2>x2);
%     cond322 = (Ucase3<Ucase2) & (beta2<gamma1) & gamma2>max(x1,0) & gamma2<x2;
%     %cond323 = (Ucase3<Ucase2) & (beta2==gamma1) & (gamma2>x1) & x1==x2;
%     %beta1=gamma1 implies Ucase1=Ucase2
%     cond32  = cond321 | cond322;
% 
%     cond331 = (Ucase3==Ucase2) & (gamma2==x1 | gamma2==x2);
%     cond332 = (Ucase3-Ucase2 < 10e-15) & (beta2==gamma1);
%     cond33 = cond331 | cond332;
% 
%     cond3   = cond31 | cond32 | cond33;
% 
%     %conditions entre u4 et u1
%     cond411 = (Ucase4>=Ucase1) & (beta1>=gamma2) & gamma1>=max(x3,0) & gamma2<=x4;
%     cond412 = (Ucase4>Ucase1) & (beta1<gamma2) & (gamma1<max(x3,0) | gamma2>x4);
%     cond41  =  cond411 | cond412;
% 
%     cond421 = (Ucase4<=Ucase1) & (beta1>=gamma2) & (gamma1<=max(x3,0) | gamma1>=x4);
%     cond422 = (Ucase4<Ucase1) & (beta1<gamma2) & gamma1>max(x3,0) & gamma1<=x4;
%     cond42  = cond421 | cond422;
% 
%     cond431 = (Ucase4==Ucase1) & (gamma1==x3 | gamma1==x4);
%     cond432 = (Ucase4-Ucase1 < 10e-15) & (beta1==gamma2);
%     cond43 = cond431 | cond432;
% 
%     cond4   = cond41 | cond42 |cond43;
% 
%     %conditions entre u1 et u3
%     cond511 = (Ucase1>=Ucase3) & (gamma1<=beta1) & (gamma2<=x5 | gamma2>=x6);
%     cond512 = (Ucase1>Ucase3) & (gamma1>beta1) & (gamma2>x6 & gamma2<x5);
%     cond51  = cond511 | cond512;
% 
%     cond521 = (Ucase1<=Ucase3) & (gamma1<=beta1) & (gamma2>=x5 & gamma2<=x6);
%     cond522 = (Ucase1<Ucase3) & (gamma1>beta1) & (gamma2<x6 | gamma2>x5);
%     cond52  = cond521 | cond522;
% 
%     cond5   = cond51 | cond52;
    condTOT = cond1&cond2&cond3&cond4&cond5;


    
end

