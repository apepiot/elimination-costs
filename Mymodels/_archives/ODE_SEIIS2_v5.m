function [dY] = ODE_SEIIS2_v5(t,Y,...
                            betaX,gammaX,nuX,epsX,sigmaX,...
                            betaY,gammaY,nuY,epsY,sigmaY,...
                            tabComp,mu,b,rho)
                        
    %probably the same that ODE_SEIIS2_v4
    
    N = b/mu; %endemic N
   
    LambdaX = betaX*Y(tabComp(tabComp.STI1=="IA" | tabComp.STI1=="IS","no"))./N;
    LambdaY = betaY*Y(tabComp(tabComp.STI2=="IA" | tabComp.STI2=="IS","no"))./N;  
    
    dY(1) = b - Y(1)*(LambdaX + LambdaY + mu) + Y(2)*rho + Y(5)*rho + Y(6)*rho + Y(7)*rho + Y(10)*rho + Y(11)*rho + Y(4)*(gammaX + nuX) + Y(13)*(gammaY + nuY) + Y(3)*(nuX + rho) + Y(9)*(nuY + rho);
    dY(2) = LambdaX*Y(1) - Y(3)*rho - Y(5)*rho - Y(6)*rho - Y(7)*rho - Y(9)*rho - Y(11)*rho + Y(10)*(nuY - rho) - Y(2)*(LambdaY + mu + rho + sigmaX) + Y(14)*(gammaY + nuY);
    dY(3) = Y(11)*(nuY - rho) - Y(5)*rho - Y(6)*rho - Y(7)*rho - Y(9)*rho - Y(10)*rho - Y(2)*(rho + sigmaX*(epsX - 1)) - Y(3)*(LambdaY + mu + nuX + rho) + Y(15)*(gammaY + nuY);
    dY(4) = Y(12)*nuY - Y(4)*(LambdaY + gammaX + mu + nuX) + Y(16)*(gammaY + nuY) + Y(2)*epsX*sigmaX;
    dY(5) = LambdaY*Y(1) - Y(2)*rho - Y(3)*rho - Y(6)*rho - Y(9)*rho - Y(10)*rho - Y(11)*rho + Y(7)*(nuX - rho) - Y(5)*(LambdaX + mu + rho + sigmaY) + Y(8)*(gammaX + nuX);
    dY(6) = Y(2)*(LambdaY - rho) - Y(7)*rho - Y(9)*rho - Y(10)*rho - Y(11)*rho - Y(3)*rho + Y(5)*(LambdaX - rho) - Y(6)*(mu + rho + sigmaX + sigmaY);
    dY(7) = Y(3)*(LambdaY - rho) - Y(2)*rho - Y(5)*rho - Y(9)*rho - Y(10)*rho - Y(11)*rho - Y(6)*(rho + sigmaX*(epsX - 1)) - Y(7)*(mu + nuX + rho + sigmaY);
    dY(8) = LambdaY*Y(4) - Y(8)*(gammaX + mu + nuX + sigmaY) + Y(6)*epsX*sigmaX;
    dY(9) = Y(11)*(nuX - rho) - Y(2)*rho - Y(3)*rho - Y(6)*rho - Y(7)*rho - Y(10)*rho - Y(5)*(rho + sigmaY*(epsY - 1)) - Y(9)*(LambdaX + mu + nuY + rho) + Y(12)*(gammaX + nuX);
    dY(10) = Y(9)*(LambdaX - rho) - Y(2)*rho - Y(3)*rho - Y(5)*rho - Y(7)*rho - Y(11)*rho - Y(6)*(rho + sigmaY*(epsY - 1)) - Y(10)*(mu + nuY + rho + sigmaX);
    dY(11) = - Y(7)*(rho + sigmaY*(epsY - 1)) - Y(10)*(rho + sigmaX*(epsX - 1)) - Y(2)*rho - Y(3)*rho - Y(5)*rho - Y(6)*rho - Y(9)*rho - Y(11)*(mu + nuX + nuY + rho);
    dY(12) = Y(10)*epsX*sigmaX - Y(12)*(gammaX + mu + nuX + nuY) - Y(8)*sigmaY*(epsY - 1);
    dY(13) = Y(15)*nuX - Y(13)*(LambdaX + gammaY + mu + nuY) + Y(16)*(gammaX + nuX) + Y(5)*epsY*sigmaY;
    dY(14) = LambdaX*Y(13) - Y(14)*(gammaY + mu + nuY + sigmaX) + Y(6)*epsY*sigmaY;
    dY(15) = Y(7)*epsY*sigmaY - Y(15)*(gammaY + mu + nuX + nuY) - Y(14)*sigmaX*(epsX - 1);
    dY(16) = Y(8)*epsY*sigmaY - Y(16)*(gammaX + gammaY + mu + nuX + nuY) + Y(14)*epsX*sigmaX;
    
    %         STI1    STI2    no        X    
    %     ____    ____    __    _________
    % 
    %     "S"     "S"      1    [1×1 sym]
    %     "E"     "S"      2    [1×1 sym]
    %     "IA"    "S"      3    [1×1 sym]
    %     "IS"    "S"      4    [1×1 sym]
    %     "S"     "E"      5    [1×1 sym]
    %     "E"     "E"      6    [1×1 sym]
    %     "IA"    "E"      7    [1×1 sym]
    %     "IS"    "E"      8    [1×1 sym]
    %     "S"     "IA"     9    [1×1 sym]
    %     "E"     "IA"    10    [1×1 sym]
    %     "IA"    "IA"    11    [1×1 sym]
    %     "IS"    "IA"    12    [1×1 sym]
    %     "S"     "IS"    13    [1×1 sym]
    %     "E"     "IS"    14    [1×1 sym]
    %     "IA"    "IS"    15    [1×1 sym]
    %     "IS"    "IS"    16    [1×1 sym]

end

