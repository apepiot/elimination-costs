%Verif solution i(t) de Wolfram Alpha
close all; clear all;
p=0.2;
mu=0.05;
beta=0.6;
gamma=0.1;
N0=15;
I0=3;

beta+(mu+gamma)*p/mu

i_t=[];
vecTime = 0:0.1:1000;
for t=vecTime
    i_t = [i_t, i_WA(t, p, mu, beta, gamma, N0, I0)];
end

figure(1)
plot(vecTime, i_t)


%Par le modele SIS directement

%load ('/Users/amandine/Desktop/These/Codes/Basiques/SIS.m')

% *Computing the ODE's system*
%MaxTime = 100;
S0=N0-I0;
[tps, pop]= ode45(@(t,y) SIS(t,y,[beta gamma p mu], 'density'),[vecTime],[S0 I0]);
S = pop(:,1); I = pop(:,2);

figure(2)
plot(vecTime, I)


figure(3)

plot(vecTime, abs(i_t' - I))

diff = abs(i_t' - I);
diff(end)