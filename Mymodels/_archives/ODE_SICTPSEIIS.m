function [dY] = ODE_SICTPSEIIS(t,Y,betaIh,betaCh,sigmah,gammah,zetah,etah,ph,...
                                betaX,gammaX,nuX,epsX,sigmaX,...
                                tabComp,mu,b,rho)
    %forget to include mandatory routine testing rate of stis under prep
    %N = sum(Y);
    N = b/mu; %endemic N
    dY = zeros(length(Y),1);
    
    %Lambdah = betaIh*sum(Y(table2array(tabComp(tabComp.HIV1=="I" | tabComp.HIV1=="Ip","no"))))./N +...
    %    betaCh*sum(Y(table2array(tabComp(tabComp.HIV1=="C" | tabComp.HIV1=="Cp","no"))))./N;
    %LambdaX = betaX*Y(table2array(tabComp(tabComp.STI1=="IA" | tabComp.STI1=="IS","no")))./N;
    
    Lambdah = betaIh*sum(Y([2     5     9    12    16    19    23    26]))./N +...
        betaCh*sum(Y([3     6    10    13    17    20    24    27]))./N;
    LambdaX = betaX*sum(Y([15    16    17    18    19    20    21    22    23    24    25    26    27    28]))./N;
    
    dY(1) = Y(8)*rho - Y(1)*(Lambdah + LambdaX + mu) - b*(ph - 1) + Y(22)*(gammaX + nuX) + Y(15)*(nuX + rho);
    dY(2) = Lambdah*Y(1) + Y(16)*nuX - Y(2)*(LambdaX + mu + rho + sigmah) + Y(23)*(gammaX + nuX);
    dY(3) = Y(17)*nuX + Y(2)*sigmah - Y(3)*(LambdaX + gammah + mu + rho) + Y(24)*(gammaX + nuX);
    dY(4) = Y(11)*rho + b*ph - Y(4)*(LambdaX + mu - Lambdah*(zetah - 1)) + Y(25)*(gammaX + nuX) + Y(18)*(nuX + rho);
    dY(5) = Y(19)*nuX - Y(5)*(LambdaX + etah + mu + rho + sigmah) + Y(26)*(gammaX + nuX) - Lambdah*Y(4)*(zetah - 1);
    dY(6) = Y(20)*nuX - Y(6)*(LambdaX + etah + gammah + mu + rho) + Y(5)*sigmah + Y(27)*(gammaX + nuX);
    dY(7) = Y(2)*rho + Y(9)*rho + Y(10)*rho + Y(12)*rho + Y(13)*rho + Y(14)*rho + Y(16)*rho + Y(17)*rho + Y(19)*rho + Y(20)*rho + Y(6)*(etah + gammah + rho) - Y(7)*(LambdaX + mu) + Y(28)*(gammaX + nuX) + Y(5)*(etah + rho) + Y(3)*(gammah + rho) + Y(21)*(nuX + rho);
    dY(8) = LambdaX*Y(1) - Y(8)*(Lambdah + mu + rho + sigmaX);
    dY(9) = Lambdah*Y(8) + LambdaX*Y(2) - Y(9)*(mu + rho + sigmah + sigmaX);
    dY(10) = LambdaX*Y(3) + Y(9)*sigmah - Y(10)*(gammah + mu + rho + sigmaX);
    dY(11) = LambdaX*Y(4) - Y(11)*(mu + rho + sigmaX - Lambdah*(zetah - 1));
    dY(12) = LambdaX*Y(5) - Y(12)*(etah + mu + rho + sigmah + sigmaX) - Lambdah*Y(11)*(zetah - 1);
    dY(13) = LambdaX*Y(6) + Y(12)*sigmah - Y(13)*(etah + gammah + mu + rho + sigmaX);
    dY(14) = LambdaX*Y(7) + Y(12)*etah + Y(10)*gammah - Y(14)*(mu + rho + sigmaX) + Y(13)*(etah + gammah);
    dY(15) = - Y(15)*(Lambdah + mu + nuX + rho) - Y(8)*sigmaX*(epsX - 1);
    dY(16) = Lambdah*Y(15) - Y(16)*(mu + nuX + rho + sigmah) - Y(9)*sigmaX*(epsX - 1);
    dY(17) = Y(16)*sigmah - Y(17)*(gammah + mu + nuX + rho) - Y(10)*sigmaX*(epsX - 1);
    dY(18) = - Y(18)*(mu + nuX + rho - Lambdah*(zetah - 1)) - Y(11)*sigmaX*(epsX - 1);
    dY(19) = - Y(19)*(etah + mu + nuX + rho + sigmah) - Lambdah*Y(18)*(zetah - 1) - Y(12)*sigmaX*(epsX - 1);
    dY(20) = Y(19)*sigmah - Y(20)*(etah + gammah + mu + nuX + rho) - Y(13)*sigmaX*(epsX - 1);
    dY(21) = Y(19)*etah + Y(17)*gammah - Y(21)*(mu + nuX + rho) + Y(20)*(etah + gammah) - Y(14)*sigmaX*(epsX - 1);
    dY(22) = Y(8)*epsX*sigmaX - Y(22)*(Lambdah + gammaX + mu + nuX);
    dY(23) = Lambdah*Y(22) - Y(23)*(gammaX + mu + nuX + sigmah) + Y(9)*epsX*sigmaX;
    dY(24) = Y(23)*sigmah - Y(24)*(gammah + gammaX + mu + nuX) + Y(10)*epsX*sigmaX;
    dY(25) = Y(11)*epsX*sigmaX - Y(25)*(gammaX + mu + nuX - Lambdah*(zetah - 1));
    dY(26) = Y(12)*epsX*sigmaX - Y(26)*(etah + gammaX + mu + nuX + sigmah) - Lambdah*Y(25)*(zetah - 1);
    dY(27) = Y(26)*sigmah - Y(27)*(etah + gammah + gammaX + mu + nuX) + Y(13)*epsX*sigmaX;
    dY(28) = Y(26)*etah + Y(24)*gammah - Y(28)*(gammaX + mu + nuX) + Y(27)*(etah + gammah) + Y(14)*epsX*sigmaX;
    
    %         HIV1    STI1    no        X    
    %     ____    ____    __    _________
    % 
    %     "S"     "S"      1    [1×1 sym]
    %     "I"     "S"      2    [1×1 sym]
    %     "C"     "S"      3    [1×1 sym]
    %     "P"     "S"      4    [1×1 sym]
    %     "Ip"    "S"      5    [1×1 sym]
    %     "Cp"    "S"      6    [1×1 sym]
    %     "T"     "S"      7    [1×1 sym]
    %     "S"     "E"      8    [1×1 sym]
    %     "I"     "E"      9    [1×1 sym]
    %     "C"     "E"     10    [1×1 sym]
    %     "P"     "E"     11    [1×1 sym]
    %     "Ip"    "E"     12    [1×1 sym]
    %     "Cp"    "E"     13    [1×1 sym]
    %     "T"     "E"     14    [1×1 sym]
    %     "S"     "IA"    15    [1×1 sym]
    %     "I"     "IA"    16    [1×1 sym]
    %     "C"     "IA"    17    [1×1 sym]
    %     "P"     "IA"    18    [1×1 sym]
    %     "Ip"    "IA"    19    [1×1 sym]
    %     "Cp"    "IA"    20    [1×1 sym]
    %     "T"     "IA"    21    [1×1 sym]
    %     "S"     "IS"    22    [1×1 sym]
    %     "I"     "IS"    23    [1×1 sym]
    %     "C"     "IS"    24    [1×1 sym]
    %     "P"     "IS"    25    [1×1 sym]
    %     "Ip"    "IS"    26    [1×1 sym]
    %     "Cp"    "IS"    27    [1×1 sym]
    %     "T"     "IS"    28    [1×1 sym]
end

