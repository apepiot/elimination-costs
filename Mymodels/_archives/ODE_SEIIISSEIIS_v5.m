function [dY] = ODE_SEIIISSEIIS_v5(t,Y,...
                            betas,sigmas,gamma3s,taus,...
                            betaX,gammaX,nuX,epsX,...
                            tabComp,mu,b,rho)
    %probably the same that SEIIISSEIIS_v4
    
    N = b/mu; %endemic N
    Lambdas = betas*Y(tabComp(tabComp.syph1=="I1" | tabComp.syph1=="I2" | tabComp.syph1=="I3" ,"no"))./N;
    LambdaX = betaX*Y(tabComp(tabComp.STI1=="IA" | tabComp.STI1=="IS","no"))./N;

    dY(1) = b - Y(1)*(Lambdas + LambdaX + mu) + Y(2)*rho + Y(3)*rho + Y(4)*rho + Y(6)*rho + Y(7)*rho + Y(8)*rho + Y(9)*rho + Y(10)*rho + Y(12)*rho + Y(13)*rho + Y(14)*rho + Y(15)*rho + Y(16)*(gammaX + nuX) + Y(5)*(gamma3s + rho) + Y(11)*(nuX + rho);
    dY(2) = Lambdas*Y(1) - Y(3)*rho - Y(4)*rho - Y(5)*rho - Y(6)*rho - Y(7)*rho - Y(8)*rho - Y(9)*rho - Y(10)*rho - Y(11)*rho - Y(13)*rho - Y(14)*rho - Y(15)*rho + Y(12)*(nuX - rho) - Y(2)*(LambdaX + mu + rho + sigmas) + Y(17)*(gammaX + nuX);
    dY(3) = Y(13)*(nuX - rho) - Y(5)*rho - Y(6)*rho - Y(7)*rho - Y(8)*rho - Y(9)*rho - Y(10)*rho - Y(11)*rho - Y(12)*rho - Y(14)*rho - Y(15)*rho - Y(4)*rho - Y(2)*(rho - sigmas) - Y(3)*(LambdaX + mu + rho + taus) + Y(18)*(gammaX + nuX);
    dY(4) = Y(14)*(nuX - rho) - Y(5)*rho - Y(6)*rho - Y(7)*rho - Y(8)*rho - Y(9)*rho - Y(10)*rho - Y(11)*rho - Y(12)*rho - Y(13)*rho - Y(15)*rho - Y(2)*rho - Y(3)*(rho - taus) - Y(4)*(LambdaX + mu + rho + thetas) + Y(19)*(gammaX + nuX);
    dY(5) = Y(15)*(nuX - rho) - Y(3)*rho - Y(6)*rho - Y(7)*rho - Y(8)*rho - Y(9)*rho - Y(10)*rho - Y(11)*rho - Y(12)*rho - Y(13)*rho - Y(14)*rho - Y(2)*rho - Y(4)*(rho - thetas) - Y(5)*(LambdaX + gamma3s + mu + rho) + Y(20)*(gammaX + nuX);
    dY(6) = LambdaX*Y(1) - Y(2)*rho - Y(3)*rho - Y(4)*rho - Y(5)*rho - Y(7)*rho - Y(8)*rho - Y(9)*rho - Y(11)*rho - Y(12)*rho - Y(13)*rho - Y(14)*rho - Y(15)*rho + Y(10)*(gamma3s - rho) - Y(6)*(Lambdas + mu + rho + sigmaX);
    dY(7) = Y(6)*(Lambdas - rho) - Y(4)*rho - Y(5)*rho - Y(8)*rho - Y(9)*rho - Y(10)*rho - Y(11)*rho - Y(12)*rho - Y(13)*rho - Y(14)*rho - Y(15)*rho - Y(3)*rho + Y(2)*(LambdaX - rho) - Y(7)*(mu + rho + sigmas + sigmaX);
    dY(8) = Y(3)*(LambdaX - rho) - Y(4)*rho - Y(5)*rho - Y(6)*rho - Y(9)*rho - Y(10)*rho - Y(11)*rho - Y(12)*rho - Y(13)*rho - Y(14)*rho - Y(15)*rho - Y(2)*rho - Y(7)*(rho - sigmas) - Y(8)*(mu + rho + sigmaX + taus);
    dY(9) = Y(4)*(LambdaX - rho) - Y(3)*rho - Y(5)*rho - Y(6)*rho - Y(7)*rho - Y(10)*rho - Y(11)*rho - Y(12)*rho - Y(13)*rho - Y(14)*rho - Y(15)*rho - Y(2)*rho - Y(8)*(rho - taus) - Y(9)*(mu + rho + sigmaX + thetas);
    dY(10) = Y(5)*(LambdaX - rho) - Y(3)*rho - Y(4)*rho - Y(6)*rho - Y(7)*rho - Y(8)*rho - Y(11)*rho - Y(12)*rho - Y(13)*rho - Y(14)*rho - Y(15)*rho - Y(2)*rho - Y(9)*(rho - thetas) - Y(10)*(gamma3s + mu + rho + sigmaX);
    dY(11) = Y(15)*(gamma3s - rho) - Y(2)*rho - Y(3)*rho - Y(4)*rho - Y(5)*rho - Y(7)*rho - Y(8)*rho - Y(9)*rho - Y(10)*rho - Y(12)*rho - Y(13)*rho - Y(14)*rho - Y(6)*(rho + sigmaX*(epsX - 1)) - Y(11)*(Lambdas + mu + nuX + rho);
    dY(12) = Y(11)*(Lambdas - rho) - Y(2)*rho - Y(3)*rho - Y(4)*rho - Y(5)*rho - Y(6)*rho - Y(8)*rho - Y(9)*rho - Y(10)*rho - Y(13)*rho - Y(14)*rho - Y(15)*rho - Y(7)*(rho + sigmaX*(epsX - 1)) - Y(12)*(mu + nuX + rho + sigmas);
    dY(13) = - Y(8)*(rho + sigmaX*(epsX - 1)) - Y(2)*rho - Y(3)*rho - Y(4)*rho - Y(5)*rho - Y(6)*rho - Y(7)*rho - Y(9)*rho - Y(10)*rho - Y(11)*rho - Y(14)*rho - Y(15)*rho - Y(12)*(rho - sigmas) - Y(13)*(mu + nuX + rho + taus);
    dY(14) = - Y(9)*(rho + sigmaX*(epsX - 1)) - Y(2)*rho - Y(3)*rho - Y(4)*rho - Y(5)*rho - Y(6)*rho - Y(7)*rho - Y(8)*rho - Y(10)*rho - Y(11)*rho - Y(12)*rho - Y(15)*rho - Y(13)*(rho - taus) - Y(14)*(mu + nuX + rho + thetas);
    dY(15) = - Y(10)*(rho + sigmaX*(epsX - 1)) - Y(2)*rho - Y(3)*rho - Y(4)*rho - Y(5)*rho - Y(6)*rho - Y(7)*rho - Y(8)*rho - Y(9)*rho - Y(11)*rho - Y(12)*rho - Y(13)*rho - Y(14)*(rho - thetas) - Y(15)*(gamma3s + mu + nuX + rho);
    dY(16) = Y(20)*gamma3s - Y(16)*(Lambdas + gammaX + mu + nuX) + Y(6)*epsX*sigmaX;
    dY(17) = Lambdas*Y(16) - Y(17)*(gammaX + mu + nuX + sigmas) + Y(7)*epsX*sigmaX;
    dY(18) = Y(17)*sigmas - Y(18)*(gammaX + mu + nuX + taus) + Y(8)*epsX*sigmaX;
    dY(19) = Y(18)*taus - Y(19)*(gammaX + mu + nuX + thetas) + Y(9)*epsX*sigmaX;
    dY(20) = Y(19)*thetas - Y(20)*(gamma3s + gammaX + mu + nuX) + Y(10)*epsX*sigmaX;
    
    %         syph1    STI1    no        X    
    %     _____    ____    __    _________
    % 
    %     "S"      "S"      1    [1×1 sym]
    %     "E"      "S"      2    [1×1 sym]
    %     "I1"     "S"      3    [1×1 sym]
    %     "I2"     "S"      4    [1×1 sym]
    %     "I3"     "S"      5    [1×1 sym]
    %     "S"      "E"      6    [1×1 sym]
    %     "E"      "E"      7    [1×1 sym]
    %     "I1"     "E"      8    [1×1 sym]
    %     "I2"     "E"      9    [1×1 sym]
    %     "I3"     "E"     10    [1×1 sym]
    %     "S"      "IA"    11    [1×1 sym]
    %     "E"      "IA"    12    [1×1 sym]
    %     "I1"     "IA"    13    [1×1 sym]
    %     "I2"     "IA"    14    [1×1 sym]
    %     "I3"     "IA"    15    [1×1 sym]
    %     "S"      "IS"    16    [1×1 sym]
    %     "E"      "IS"    17    [1×1 sym]
    %     "I1"     "IS"    18    [1×1 sym]
    %     "I2"     "IS"    19    [1×1 sym]
    %     "I3"     "IS"    20    [1×1 sym]
end

