function [dY] = ODE_SICTPSEIIS2(t,Y,betaIh,betaCh,sigmah,gammah,zetah,etah,ph,...
                            betaX,gammaX,nuX,epsX,sigmaX,...
                            betaY,gammaY,nuY,epsY,sigmaY,...
                            tabComp,mu,b,rho)
    %forget to include mandatory routine testing rate of stis under prep
     dY = zeros(length(Y),1);
     %N = sum(Y);
     N = b/mu; %endemic N
     
%      Lambdah = betaIh*sum(Y(table2array(tabComp(tabComp.HIV1=="I" | tabComp.HIV1=="Ip","no"))))./N +...
%         betaCh*sum(Y(table2array(tabComp(tabComp.HIV1=="C" | tabComp.HIV1=="Cp","no"))))./N;
%      LambdaX = betaX*sum(Y(table2array(tabComp(tabComp.STI1=="IA" | tabComp.STI1=="IS","no"))))./N;
%      LambdaY = betaY*sum(Y(table2array(tabComp(tabComp.STI2=="IA" | tabComp.STI2=="IS","no"))))./N;
     
     Lambdah = betaIh*sum(Y([2     5     9    12    16    19    23    26    30    33    37    40    44    47    51    54    58    61    65    68    72 ...
         75    79    82    86    89    93    96   100   103   107   110]))./N +...
        betaCh*sum(Y([3     6    10    13    17    20    24    27    31    34    38    41    45    48    52    55    59    62    66    69    73 ...
        76    80    83    87    90    94    97   101   104   108   111]))./N;
     LambdaX = betaX*sum(Y([15    16    17    18    19    20    21    22    23    24    25    26    27    28    43    44    45    46    47    48    49 ...
         50    51    52    53    54    55    56    71    72    73    74    75    76    77    78    79    80    81    82    83    84 ...
         99   100   101   102   103   104   105   106   107   108   109   110   111   112]))./N;
     LambdaY = betaY*sum(Y([57    58    59    60    61    62    63    64    65    66    67    68    69    70    71    72    73    74    75    76    77 ...
         78    79    80    81    82    83    84    85    86    87    88    89    90    91    92    93    94    95    96    97    98 ...
         99   100   101   102   103   104   105   106   107   108   109   110   111   112]))./N;
     
     dY(1) = Y(8)*rho + Y(29)*rho + Y(36)*rho + Y(43)*rho + Y(64)*rho + Y(71)*rho - Y(1)*(Lambdah + LambdaX + LambdaY + mu) - b*(ph - 1) + Y(22)*(gammaX + nuX) + Y(85)*(gammaY + nuY) + Y(15)*(nuX + rho) + Y(57)*(nuY + rho);
     dY(2) = Lambdah*Y(1) - Y(2)*(LambdaX + LambdaY + mu + rho + sigmah) + Y(16)*nuX + Y(58)*nuY + Y(23)*(gammaX + nuX) + Y(86)*(gammaY + nuY);
     dY(3) = Y(17)*nuX + Y(59)*nuY + Y(2)*sigmah + Y(24)*(gammaX + nuX) + Y(87)*(gammaY + nuY) - Y(3)*(LambdaX + LambdaY + gammah + mu + rho);
     dY(4) = Y(11)*rho + Y(32)*rho + Y(39)*rho + Y(46)*rho + Y(67)*rho + Y(74)*rho + b*ph - Y(4)*(LambdaX + LambdaY + mu - Lambdah*(zetah - 1)) + Y(25)*(gammaX + nuX) + Y(88)*(gammaY + nuY) + Y(18)*(nuX + rho) + Y(60)*(nuY + rho);
     dY(5) = Y(19)*nuX + Y(61)*nuY - Y(5)*(LambdaX + LambdaY + etah + mu + rho + sigmah) + Y(26)*(gammaX + nuX) + Y(89)*(gammaY + nuY) - Lambdah*Y(4)*(zetah - 1);
     dY(6) = Y(20)*nuX + Y(62)*nuY + Y(5)*sigmah - Y(6)*(LambdaX + LambdaY + etah + gammah + mu + rho) + Y(27)*(gammaX + nuX) + Y(90)*(gammaY + nuY);
     dY(7) = Y(2)*rho - Y(7)*(LambdaX + LambdaY + mu) + Y(9)*rho + Y(10)*rho + Y(12)*rho + Y(13)*rho + Y(14)*rho + Y(16)*rho + Y(17)*rho + Y(19)*rho + Y(20)*rho + Y(30)*rho + Y(31)*rho + Y(33)*rho + Y(34)*rho + Y(35)*rho + Y(37)*rho + Y(38)*rho + Y(40)*rho + Y(41)*rho + Y(42)*rho + Y(44)*rho + Y(45)*rho + Y(47)*rho + Y(48)*rho + Y(49)*rho + Y(58)*rho + Y(59)*rho + Y(61)*rho + Y(62)*rho + Y(65)*rho + Y(66)*rho + Y(68)*rho + Y(69)*rho + Y(70)*rho + Y(72)*rho + Y(73)*rho + Y(75)*rho + Y(76)*rho + Y(77)*rho + Y(6)*(etah + gammah + rho) + Y(28)*(gammaX + nuX) + Y(91)*(gammaY + nuY) + Y(5)*(etah + rho) + Y(3)*(gammah + rho) + Y(21)*(nuX + rho) + Y(63)*(nuY + rho);
     dY(8) = LambdaX*Y(1) - Y(8)*(Lambdah + LambdaY + mu + rho + sigmaX) + Y(64)*nuY + Y(92)*(gammaY + nuY);
     dY(9) = Lambdah*Y(8) + LambdaX*Y(2) + Y(65)*nuY - Y(9)*(LambdaY + mu + rho + sigmah + sigmaX) + Y(93)*(gammaY + nuY);
     dY(10) = LambdaX*Y(3) - Y(10)*(LambdaY + gammah + mu + rho + sigmaX) + Y(66)*nuY + Y(9)*sigmah + Y(94)*(gammaY + nuY);
     dY(11) = LambdaX*Y(4) + Y(67)*nuY - Y(11)*(LambdaY + mu + rho + sigmaX - Lambdah*(zetah - 1)) + Y(95)*(gammaY + nuY);
     dY(12) = LambdaX*Y(5) + Y(68)*nuY - Y(12)*(LambdaY + etah + mu + rho + sigmah + sigmaX) + Y(96)*(gammaY + nuY) - Lambdah*Y(11)*(zetah - 1);
     dY(13) = LambdaX*Y(6) + Y(69)*nuY + Y(12)*sigmah - Y(13)*(LambdaY + etah + gammah + mu + rho + sigmaX) + Y(97)*(gammaY + nuY);
     dY(14) = LambdaX*Y(7) + Y(12)*etah + Y(10)*gammah + Y(70)*nuY - Y(14)*(LambdaY + mu + rho + sigmaX) + Y(13)*(etah + gammah) + Y(98)*(gammaY + nuY);
     dY(15) = Y(71)*nuY + Y(99)*(gammaY + nuY) - Y(15)*(Lambdah + LambdaY + mu + nuX + rho) - Y(8)*sigmaX*(epsX - 1);
     dY(16) = Lambdah*Y(15) - Y(16)*(LambdaY + mu + nuX + rho + sigmah) + Y(72)*nuY + Y(100)*(gammaY + nuY) - Y(9)*sigmaX*(epsX - 1);
     dY(17) = Y(73)*nuY - Y(17)*(LambdaY + gammah + mu + nuX + rho) + Y(16)*sigmah + Y(101)*(gammaY + nuY) - Y(10)*sigmaX*(epsX - 1);
     dY(18) = Y(74)*nuY - Y(18)*(LambdaY + mu + nuX + rho - Lambdah*(zetah - 1)) + Y(102)*(gammaY + nuY) - Y(11)*sigmaX*(epsX - 1);
     dY(19) = Y(75)*nuY - Y(19)*(LambdaY + etah + mu + nuX + rho + sigmah) + Y(103)*(gammaY + nuY) - Lambdah*Y(18)*(zetah - 1) - Y(12)*sigmaX*(epsX - 1);
     dY(20) = Y(76)*nuY + Y(19)*sigmah - Y(20)*(LambdaY + etah + gammah + mu + nuX + rho) + Y(104)*(gammaY + nuY) - Y(13)*sigmaX*(epsX - 1);
     dY(21) = Y(19)*etah + Y(17)*gammah + Y(77)*nuY - Y(21)*(LambdaY + mu + nuX + rho) + Y(20)*(etah + gammah) + Y(105)*(gammaY + nuY) - Y(14)*sigmaX*(epsX - 1);
     dY(22) = Y(78)*nuY + Y(106)*(gammaY + nuY) - Y(22)*(Lambdah + LambdaY + gammaX + mu + nuX) + Y(8)*epsX*sigmaX;
     dY(23) = Lambdah*Y(22) - Y(23)*(LambdaY + gammaX + mu + nuX + sigmah) + Y(79)*nuY + Y(107)*(gammaY + nuY) + Y(9)*epsX*sigmaX;
     dY(24) = Y(80)*nuY - Y(24)*(LambdaY + gammah + gammaX + mu + nuX) + Y(23)*sigmah + Y(108)*(gammaY + nuY) + Y(10)*epsX*sigmaX;
     dY(25) = Y(81)*nuY - Y(25)*(LambdaY + gammaX + mu + nuX - Lambdah*(zetah - 1)) + Y(109)*(gammaY + nuY) + Y(11)*epsX*sigmaX;
     dY(26) = Y(82)*nuY - Y(26)*(LambdaY + etah + gammaX + mu + nuX + sigmah) + Y(110)*(gammaY + nuY) + Y(12)*epsX*sigmaX - Lambdah*Y(25)*(zetah - 1);
     dY(27) = Y(83)*nuY + Y(26)*sigmah - Y(27)*(LambdaY + etah + gammah + gammaX + mu + nuX) + Y(111)*(gammaY + nuY) + Y(13)*epsX*sigmaX;
     dY(28) = Y(26)*etah + Y(24)*gammah + Y(84)*nuY - Y(28)*(LambdaY + gammaX + mu + nuX) + Y(27)*(etah + gammah) + Y(112)*(gammaY + nuY) + Y(14)*epsX*sigmaX;
     dY(29) = LambdaY*Y(1) - Y(29)*(Lambdah + LambdaX + mu + rho + sigmaY) + Y(43)*nuX + Y(50)*(gammaX + nuX);
     dY(30) = Lambdah*Y(29) + LambdaY*Y(2) + Y(44)*nuX - Y(30)*(LambdaX + mu + rho + sigmah + sigmaY) + Y(51)*(gammaX + nuX);
     dY(31) = LambdaY*Y(3) - Y(31)*(LambdaX + gammah + mu + rho + sigmaY) + Y(45)*nuX + Y(30)*sigmah + Y(52)*(gammaX + nuX);
     dY(32) = LambdaY*Y(4) + Y(46)*nuX - Y(32)*(LambdaX + mu + rho + sigmaY - Lambdah*(zetah - 1)) + Y(53)*(gammaX + nuX);
     dY(33) = LambdaY*Y(5) + Y(47)*nuX - Y(33)*(LambdaX + etah + mu + rho + sigmah + sigmaY) + Y(54)*(gammaX + nuX) - Lambdah*Y(32)*(zetah - 1);
     dY(34) = LambdaY*Y(6) + Y(48)*nuX + Y(33)*sigmah - Y(34)*(LambdaX + etah + gammah + mu + rho + sigmaY) + Y(55)*(gammaX + nuX);
     dY(35) = LambdaY*Y(7) + Y(33)*etah + Y(31)*gammah + Y(49)*nuX - Y(35)*(LambdaX + mu + rho + sigmaY) + Y(34)*(etah + gammah) + Y(56)*(gammaX + nuX);
     dY(36) = LambdaY*Y(8) + LambdaX*Y(29) - Y(36)*(Lambdah + mu + rho + sigmaX + sigmaY);
     dY(37) = Lambdah*Y(36) + LambdaY*Y(9) + LambdaX*Y(30) - Y(37)*(mu + rho + sigmah + sigmaX + sigmaY);
     dY(38) = LambdaY*Y(10) + LambdaX*Y(31) + Y(37)*sigmah - Y(38)*(gammah + mu + rho + sigmaX + sigmaY);
     dY(39) = LambdaY*Y(11) + LambdaX*Y(32) - Y(39)*(mu + rho + sigmaX + sigmaY - Lambdah*(zetah - 1));
     dY(40) = LambdaY*Y(12) + LambdaX*Y(33) - Y(40)*(etah + mu + rho + sigmah + sigmaX + sigmaY) - Lambdah*Y(39)*(zetah - 1);
     dY(41) = LambdaY*Y(13) + LambdaX*Y(34) + Y(40)*sigmah - Y(41)*(etah + gammah + mu + rho + sigmaX + sigmaY);
     dY(42) = LambdaY*Y(14) + LambdaX*Y(35) + Y(40)*etah + Y(38)*gammah - Y(42)*(mu + rho + sigmaX + sigmaY) + Y(41)*(etah + gammah);
     dY(43) = LambdaY*Y(15) - Y(43)*(Lambdah + mu + nuX + rho + sigmaY) - Y(36)*sigmaX*(epsX - 1);
     dY(44) = Lambdah*Y(43) + LambdaY*Y(16) - Y(44)*(mu + nuX + rho + sigmah + sigmaY) - Y(37)*sigmaX*(epsX - 1);
     dY(45) = LambdaY*Y(17) + Y(44)*sigmah - Y(45)*(gammah + mu + nuX + rho + sigmaY) - Y(38)*sigmaX*(epsX - 1);
     dY(46) = LambdaY*Y(18) - Y(46)*(mu + nuX + rho + sigmaY - Lambdah*(zetah - 1)) - Y(39)*sigmaX*(epsX - 1);
     dY(47) = LambdaY*Y(19) - Y(47)*(etah + mu + nuX + rho + sigmah + sigmaY) - Lambdah*Y(46)*(zetah - 1) - Y(40)*sigmaX*(epsX - 1);
     dY(48) = LambdaY*Y(20) + Y(47)*sigmah - Y(48)*(etah + gammah + mu + nuX + rho + sigmaY) - Y(41)*sigmaX*(epsX - 1);
     dY(49) = LambdaY*Y(21) + Y(47)*etah + Y(45)*gammah - Y(49)*(mu + nuX + rho + sigmaY) + Y(48)*(etah + gammah) - Y(42)*sigmaX*(epsX - 1);
     dY(50) = LambdaY*Y(22) - Y(50)*(Lambdah + gammaX + mu + nuX + sigmaY) + Y(36)*epsX*sigmaX;
     dY(51) = Lambdah*Y(50) + LambdaY*Y(23) - Y(51)*(gammaX + mu + nuX + sigmah + sigmaY) + Y(37)*epsX*sigmaX;
     dY(52) = LambdaY*Y(24) + Y(51)*sigmah - Y(52)*(gammah + gammaX + mu + nuX + sigmaY) + Y(38)*epsX*sigmaX;
     dY(53) = LambdaY*Y(25) - Y(53)*(gammaX + mu + nuX + sigmaY - Lambdah*(zetah - 1)) + Y(39)*epsX*sigmaX;
     dY(54) = LambdaY*Y(26) - Y(54)*(etah + gammaX + mu + nuX + sigmah + sigmaY) + Y(40)*epsX*sigmaX - Lambdah*Y(53)*(zetah - 1);
     dY(55) = LambdaY*Y(27) + Y(54)*sigmah - Y(55)*(etah + gammah + gammaX + mu + nuX + sigmaY) + Y(41)*epsX*sigmaX;
     dY(56) = LambdaY*Y(28) + Y(54)*etah + Y(52)*gammah - Y(56)*(gammaX + mu + nuX + sigmaY) + Y(55)*(etah + gammah) + Y(42)*epsX*sigmaX;
     dY(57) = Y(71)*nuX + Y(78)*(gammaX + nuX) - Y(57)*(Lambdah + LambdaX + mu + nuY + rho) - Y(29)*sigmaY*(epsY - 1);
     dY(58) = Lambdah*Y(57) - Y(58)*(LambdaX + mu + nuY + rho + sigmah) + Y(72)*nuX + Y(79)*(gammaX + nuX) - Y(30)*sigmaY*(epsY - 1);
     dY(59) = Y(73)*nuX - Y(59)*(LambdaX + gammah + mu + nuY + rho) + Y(58)*sigmah + Y(80)*(gammaX + nuX) - Y(31)*sigmaY*(epsY - 1);
     dY(60) = Y(74)*nuX - Y(60)*(LambdaX + mu + nuY + rho - Lambdah*(zetah - 1)) + Y(81)*(gammaX + nuX) - Y(32)*sigmaY*(epsY - 1);
     dY(61) = Y(75)*nuX - Y(61)*(LambdaX + etah + mu + nuY + rho + sigmah) + Y(82)*(gammaX + nuX) - Lambdah*Y(60)*(zetah - 1) - Y(33)*sigmaY*(epsY - 1);
     dY(62) = Y(76)*nuX + Y(61)*sigmah - Y(62)*(LambdaX + etah + gammah + mu + nuY + rho) + Y(83)*(gammaX + nuX) - Y(34)*sigmaY*(epsY - 1);
     dY(63) = Y(61)*etah + Y(59)*gammah + Y(77)*nuX - Y(63)*(LambdaX + mu + nuY + rho) + Y(62)*(etah + gammah) + Y(84)*(gammaX + nuX) - Y(35)*sigmaY*(epsY - 1);
     dY(64) = LambdaX*Y(57) - Y(64)*(Lambdah + mu + nuY + rho + sigmaX) - Y(36)*sigmaY*(epsY - 1);
     dY(65) = Lambdah*Y(64) + LambdaX*Y(58) - Y(65)*(mu + nuY + rho + sigmah + sigmaX) - Y(37)*sigmaY*(epsY - 1);
     dY(66) = LambdaX*Y(59) + Y(65)*sigmah - Y(66)*(gammah + mu + nuY + rho + sigmaX) - Y(38)*sigmaY*(epsY - 1);
     dY(67) = LambdaX*Y(60) - Y(67)*(mu + nuY + rho + sigmaX - Lambdah*(zetah - 1)) - Y(39)*sigmaY*(epsY - 1);
     dY(68) = LambdaX*Y(61) - Y(68)*(etah + mu + nuY + rho + sigmah + sigmaX) - Lambdah*Y(67)*(zetah - 1) - Y(40)*sigmaY*(epsY - 1);
     dY(69) = LambdaX*Y(62) + Y(68)*sigmah - Y(69)*(etah + gammah + mu + nuY + rho + sigmaX) - Y(41)*sigmaY*(epsY - 1);
     dY(70) = LambdaX*Y(63) + Y(68)*etah + Y(66)*gammah - Y(70)*(mu + nuY + rho + sigmaX) + Y(69)*(etah + gammah) - Y(42)*sigmaY*(epsY - 1);
     dY(71) = - Y(71)*(Lambdah + mu + nuX + nuY + rho) - Y(43)*sigmaY*(epsY - 1) - Y(64)*sigmaX*(epsX - 1);
     dY(72) = Lambdah*Y(71) - Y(72)*(mu + nuX + nuY + rho + sigmah) - Y(44)*sigmaY*(epsY - 1) - Y(65)*sigmaX*(epsX - 1);
     dY(73) = Y(72)*sigmah - Y(73)*(gammah + mu + nuX + nuY + rho) - Y(45)*sigmaY*(epsY - 1) - Y(66)*sigmaX*(epsX - 1);
     dY(74) = - Y(74)*(mu + nuX + nuY + rho - Lambdah*(zetah - 1)) - Y(46)*sigmaY*(epsY - 1) - Y(67)*sigmaX*(epsX - 1);
     dY(75) = - Y(75)*(etah + mu + nuX + nuY + rho + sigmah) - Lambdah*Y(74)*(zetah - 1) - Y(47)*sigmaY*(epsY - 1) - Y(68)*sigmaX*(epsX - 1);
     dY(76) = Y(75)*sigmah - Y(76)*(etah + gammah + mu + nuX + nuY + rho) - Y(48)*sigmaY*(epsY - 1) - Y(69)*sigmaX*(epsX - 1);
     dY(77) = Y(75)*etah + Y(73)*gammah - Y(77)*(mu + nuX + nuY + rho) + Y(76)*(etah + gammah) - Y(49)*sigmaY*(epsY - 1) - Y(70)*sigmaX*(epsX - 1);
     dY(78) = Y(64)*epsX*sigmaX - Y(78)*(Lambdah + gammaX + mu + nuX + nuY) - Y(50)*sigmaY*(epsY - 1);
     dY(79) = Lambdah*Y(78) - Y(79)*(gammaX + mu + nuX + nuY + sigmah) + Y(65)*epsX*sigmaX - Y(51)*sigmaY*(epsY - 1);
     dY(80) = Y(79)*sigmah - Y(80)*(gammah + gammaX + mu + nuX + nuY) + Y(66)*epsX*sigmaX - Y(52)*sigmaY*(epsY - 1);
     dY(81) = Y(67)*epsX*sigmaX - Y(81)*(gammaX + mu + nuX + nuY - Lambdah*(zetah - 1)) - Y(53)*sigmaY*(epsY - 1);
     dY(82) = Y(68)*epsX*sigmaX - Y(82)*(etah + gammaX + mu + nuX + nuY + sigmah) - Lambdah*Y(81)*(zetah - 1) - Y(54)*sigmaY*(epsY - 1);
     dY(83) = Y(82)*sigmah - Y(83)*(etah + gammah + gammaX + mu + nuX + nuY) + Y(69)*epsX*sigmaX - Y(55)*sigmaY*(epsY - 1);
     dY(84) = Y(82)*etah + Y(80)*gammah - Y(84)*(gammaX + mu + nuX + nuY) + Y(83)*(etah + gammah) + Y(70)*epsX*sigmaX - Y(56)*sigmaY*(epsY - 1);
     dY(85) = Y(99)*nuX + Y(106)*(gammaX + nuX) - Y(85)*(Lambdah + LambdaX + gammaY + mu + nuY) + Y(29)*epsY*sigmaY;
     dY(86) = Lambdah*Y(85) - Y(86)*(LambdaX + gammaY + mu + nuY + sigmah) + Y(100)*nuX + Y(107)*(gammaX + nuX) + Y(30)*epsY*sigmaY;
     dY(87) = Y(101)*nuX - Y(87)*(LambdaX + gammah + gammaY + mu + nuY) + Y(86)*sigmah + Y(108)*(gammaX + nuX) + Y(31)*epsY*sigmaY;
     dY(88) = Y(102)*nuX - Y(88)*(LambdaX + gammaY + mu + nuY - Lambdah*(zetah - 1)) + Y(109)*(gammaX + nuX) + Y(32)*epsY*sigmaY;
     dY(89) = Y(103)*nuX - Y(89)*(LambdaX + etah + gammaY + mu + nuY + sigmah) + Y(110)*(gammaX + nuX) + Y(33)*epsY*sigmaY - Lambdah*Y(88)*(zetah - 1);
     dY(90) = Y(104)*nuX + Y(89)*sigmah - Y(90)*(LambdaX + etah + gammah + gammaY + mu + nuY) + Y(111)*(gammaX + nuX) + Y(34)*epsY*sigmaY;
     dY(91) = Y(89)*etah + Y(87)*gammah + Y(105)*nuX - Y(91)*(LambdaX + gammaY + mu + nuY) + Y(90)*(etah + gammah) + Y(112)*(gammaX + nuX) + Y(35)*epsY*sigmaY;
     dY(92) = LambdaX*Y(85) - Y(92)*(Lambdah + gammaY + mu + nuY + sigmaX) + Y(36)*epsY*sigmaY;
     dY(93) = Lambdah*Y(92) + LambdaX*Y(86) - Y(93)*(gammaY + mu + nuY + sigmah + sigmaX) + Y(37)*epsY*sigmaY;
     dY(94) = LambdaX*Y(87) + Y(93)*sigmah - Y(94)*(gammah + gammaY + mu + nuY + sigmaX) + Y(38)*epsY*sigmaY;
     dY(95) = LambdaX*Y(88) - Y(95)*(gammaY + mu + nuY + sigmaX - Lambdah*(zetah - 1)) + Y(39)*epsY*sigmaY;
     dY(96) = LambdaX*Y(89) - Y(96)*(etah + gammaY + mu + nuY + sigmah + sigmaX) + Y(40)*epsY*sigmaY - Lambdah*Y(95)*(zetah - 1);
     dY(97) = LambdaX*Y(90) + Y(96)*sigmah - Y(97)*(etah + gammah + gammaY + mu + nuY + sigmaX) + Y(41)*epsY*sigmaY;
     dY(98) = LambdaX*Y(91) + Y(96)*etah + Y(94)*gammah - Y(98)*(gammaY + mu + nuY + sigmaX) + Y(97)*(etah + gammah) + Y(42)*epsY*sigmaY;
     dY(99) = Y(43)*epsY*sigmaY - Y(99)*(Lambdah + gammaY + mu + nuX + nuY) - Y(92)*sigmaX*(epsX - 1);
     dY(100) = Lambdah*Y(99) - Y(100)*(gammaY + mu + nuX + nuY + sigmah) + Y(44)*epsY*sigmaY - Y(93)*sigmaX*(epsX - 1);
     dY(101) = Y(100)*sigmah - Y(101)*(gammah + gammaY + mu + nuX + nuY) + Y(45)*epsY*sigmaY - Y(94)*sigmaX*(epsX - 1);
     dY(102) = Y(46)*epsY*sigmaY - Y(102)*(gammaY + mu + nuX + nuY - Lambdah*(zetah - 1)) - Y(95)*sigmaX*(epsX - 1);
     dY(103) = Y(47)*epsY*sigmaY - Y(103)*(etah + gammaY + mu + nuX + nuY + sigmah) - Lambdah*Y(102)*(zetah - 1) - Y(96)*sigmaX*(epsX - 1);
     dY(104) = Y(103)*sigmah - Y(104)*(etah + gammah + gammaY + mu + nuX + nuY) + Y(48)*epsY*sigmaY - Y(97)*sigmaX*(epsX - 1);
     dY(105) = Y(103)*etah + Y(101)*gammah - Y(105)*(gammaY + mu + nuX + nuY) + Y(104)*(etah + gammah) + Y(49)*epsY*sigmaY - Y(98)*sigmaX*(epsX - 1);
     dY(106) = Y(50)*epsY*sigmaY - Y(106)*(Lambdah + gammaX + gammaY + mu + nuX + nuY) + Y(92)*epsX*sigmaX;
     dY(107) = Lambdah*Y(106) - Y(107)*(gammaX + gammaY + mu + nuX + nuY + sigmah) + Y(51)*epsY*sigmaY + Y(93)*epsX*sigmaX;
     dY(108) = Y(107)*sigmah - Y(108)*(gammah + gammaX + gammaY + mu + nuX + nuY) + Y(52)*epsY*sigmaY + Y(94)*epsX*sigmaX;
     dY(109) = Y(53)*epsY*sigmaY - Y(109)*(gammaX + gammaY + mu + nuX + nuY - Lambdah*(zetah - 1)) + Y(95)*epsX*sigmaX;
     dY(110) = Y(54)*epsY*sigmaY - Y(110)*(etah + gammaX + gammaY + mu + nuX + nuY + sigmah) + Y(96)*epsX*sigmaX - Lambdah*Y(109)*(zetah - 1);
     dY(111) = Y(110)*sigmah - Y(111)*(etah + gammah + gammaX + gammaY + mu + nuX + nuY) + Y(55)*epsY*sigmaY + Y(97)*epsX*sigmaX;
     dY(112) = Y(110)*etah + Y(108)*gammah - Y(112)*(gammaX + gammaY + mu + nuX + nuY) + Y(111)*(etah + gammah) + Y(56)*epsY*sigmaY + Y(98)*epsX*sigmaX;

end

