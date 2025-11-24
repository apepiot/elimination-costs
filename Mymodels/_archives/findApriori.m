function [Y] = findApriori(betaIh,betaCh,sigmah,thetah,zetah,etah,ph,...
                betaS,sigmaS,gamma1S,gamma3S,tauS,thetaS,nuS,...
                betaC,gammaC,nuC,epsC,sigmaC,...
                betaG,gammaG,nuG,epsG,sigmaG,...
                mu,b,...
                rho_h, rho_s,rho_c,rho_g)

Y=ones(560,1);           
%on peut calculer les valeurs theoriques, mais on verra plus tard
%hiv
X0 = ones(7,1); tspan=[1,100];
[res_sictp] = ode45(@(t,Y) ODE_SICTPrEP(t,Y,b,betaIh,betaCh,sigmah,...
                    thetah,ph,zetah,etah,mu,rho_h,'frequency'),...
                    tspan, X0);
popHIV = res_sictp.y(:,end);

%syph
X0 = ones(5,1); tspan=[1,100];
[res_seiiis] = ode45(@(t,Y) ODE_SEIIIS_v4(t,Y,betaS,sigmaS,tauS,...
                    thetaS,gamma1S,gamma3S,nuS,rho_s,mu,b),...
                    tspan, X0);
popS = res_seiiis.y(:,end);

%Ct
X0 = ones(4,1);
[res_seiis] = ode45(@(t,Y) ODE_SEIIS_v4(t,Y,betaC,nuC,gammaC,...
                    sigmaC,epsC,rho_c,mu,b),...
                    tspan, X0);
popCt = res_seiis.y(:,end);

%Ng
X0 = ones(4,1);
[res_seiis] = ode45(@(t,Y) ODE_SEIIS_v4(t,Y,betaG,nuG,gammaG,...
                    sigmaG,epsG,rho_g,mu,b),...
                    tspan, X0);
popNg = res_seiis.y(:,end);


Y(1:7:554) = popHIV(1)/(560/7); %S
Y(2:7:555) = popHIV(2)/(560/7); %I
Y(3:7:556) = popHIV(3)/(560/7); %C
Y(4:7:557) = popHIV(4)/(560/7); %P
Y(5:7:558) = popHIV(5)/(560/7); %Ip
Y(6:7:559) = popHIV(6)/(560/7); %Cp
Y(7:7:560) = popHIV(7)/(560/7); %T

Y = Y.*repmat(repelem(popS,7),560/7,1)'/(560/5);
%Y = repmat(popS(1)

%tabComp(tabComp.HIV=="S",:).no'


end

