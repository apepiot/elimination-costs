function [dY] = ODE_SICTPSEIIISSEIIS(t,Y,betaIh,betaCh,sigmah,gammah,zetah,etah,ph,...
                            betas,sigmas,gamma3s,taus,thetas,...
                            betaX,gammaX,nuX,epsX,sigmaX,...
                            tabComp,mu,b,rho)
    %forget to include mandatory routine testing rate of stis under prep
   dY = zeros(length(Y),1);
     
    %N = sum(Y);
    N = b/mu; %endemic N

%     Lambdah = betaIh*sum(Y(table2array(tabComp(tabComp.HIV1=="I" | tabComp.HIV1=="Ip","no"))))./N +...
%         betaCh*sum(Y(table2array(tabComp(tabComp.HIV1=="C" | tabComp.HIV1=="Cp","no"))))./N;
%     LambdaX = betaX*sum(Y(table2array(tabComp(tabComp.STI1=="IA" | tabComp.STI1=="IS","no"))))./N;
%     Lambdas = betas*sum(Y(table2array(tabComp(tabComp.syph1=="I1" | tabComp.syph1=="I2" | tabComp.syph1=="I3" ,"no"))))./N;

    Lambdah = betaIh*sum(Y([2     5     9    12    16    19    23    26    30    33    37    40    44    47    51    54    58    61    65    68    72 ...
        75    79    82    86    89    93    96   100   103   107   110   114   117   121   124   128   131   135   138]))./N +...
        betaCh*sum(Y([3     6    10    13    17    20    24    27    31    34    38    41    45    48    52    55    59    62    66    69    73 ...
        76    80    83    87    90    94    97   101   104   108   111   115   118   122   125   129   132   136   139]))./N;
    LambdaX = betaX*sum(Y([71    72    73    74    75    76    77    78    79    80    81    82    83    84    85    86    87    88    89    90    91 ...
        92    93    94    95    96    97    98    99   100   101   102   103   104   105   106   107   108   109   110   111   112 ...
        113   114   115   116   117   118   119   120   121   122   123   124   125   126   127   128   129   130   131   132   133 ...
        134   135   136   137   138   139   140]))./N;
    Lambdas = betas*sum(Y([15    16    17    18    19    20    21    22    23    24    25    26    27    28    29    30    31    32    33    34    35 ...
        50    51    52    53    54    55    56    57    58    59    60    61    62    63    64    65    66    67    68    69    70 ...
        85    86    87    88    89    90    91    92    93    94    95    96    97    98    99   100   101   102   103   104   105 ...
        120   121   122   123   124   125   126   127   128   129   130   131   132   133   134   135   136   137   138   139   140]))./N;

     dY(1) = Y(8)*rho + Y(15)*rho + Y(22)*rho + Y(36)*rho + Y(43)*rho + Y(50)*rho + Y(57)*rho + Y(64)*rho + Y(78)*rho + Y(85)*rho + Y(92)*rho + Y(99)*rho - Y(1)*(Lambdah + Lambdas + LambdaX + mu) - b*(ph - 1) + Y(106)*(gammaX + nuX) + Y(29)*(gamma3s + rho) + Y(71)*(nuX + rho);
     dY(2) = Lambdah*Y(1) - Y(2)*(Lambdas + LambdaX + mu + rho + sigmah) + Y(30)*gamma3s + Y(72)*nuX + Y(107)*(gammaX + nuX);
     dY(3) = Y(31)*gamma3s + Y(73)*nuX + Y(2)*sigmah + Y(108)*(gammaX + nuX) - Y(3)*(Lambdas + LambdaX + gammah + mu + rho);
     dY(4) = Y(11)*rho + Y(18)*rho + Y(25)*rho + Y(39)*rho + Y(46)*rho + Y(53)*rho + Y(60)*rho + Y(67)*rho + Y(81)*rho + Y(88)*rho + Y(95)*rho + Y(102)*rho + b*ph - Y(4)*(Lambdas + LambdaX + mu - Lambdah*(zetah - 1)) + Y(109)*(gammaX + nuX) + Y(32)*(gamma3s + rho) + Y(74)*(nuX + rho);
     dY(5) = Y(33)*gamma3s + Y(75)*nuX - Y(5)*(Lambdas + LambdaX + etah + mu + rho + sigmah) + Y(110)*(gammaX + nuX) - Lambdah*Y(4)*(zetah - 1);
     dY(6) = Y(34)*gamma3s + Y(76)*nuX + Y(5)*sigmah - Y(6)*(Lambdas + LambdaX + etah + gammah + mu + rho) + Y(111)*(gammaX + nuX);
     dY(7) = Y(2)*rho - Y(7)*(Lambdas + LambdaX + mu) + Y(9)*rho + Y(10)*rho + Y(12)*rho + Y(13)*rho + Y(14)*rho + Y(16)*rho + Y(17)*rho + Y(19)*rho + Y(20)*rho + Y(21)*rho + Y(23)*rho + Y(24)*rho + Y(26)*rho + Y(27)*rho + Y(28)*rho + Y(30)*rho + Y(31)*rho + Y(33)*rho + Y(34)*rho + Y(37)*rho + Y(38)*rho + Y(40)*rho + Y(41)*rho + Y(42)*rho + Y(44)*rho + Y(45)*rho + Y(47)*rho + Y(48)*rho + Y(49)*rho + Y(51)*rho + Y(52)*rho + Y(54)*rho + Y(55)*rho + Y(56)*rho + Y(58)*rho + Y(59)*rho + Y(61)*rho + Y(62)*rho + Y(63)*rho + Y(65)*rho + Y(66)*rho + Y(68)*rho + Y(69)*rho + Y(70)*rho + Y(72)*rho + Y(73)*rho + Y(75)*rho + Y(76)*rho + Y(79)*rho + Y(80)*rho + Y(82)*rho + Y(83)*rho + Y(84)*rho + Y(86)*rho + Y(87)*rho + Y(89)*rho + Y(90)*rho + Y(91)*rho + Y(93)*rho + Y(94)*rho + Y(96)*rho + Y(97)*rho + Y(98)*rho + Y(100)*rho + Y(101)*rho + Y(103)*rho + Y(104)*rho + Y(105)*rho + Y(6)*(etah + gammah + rho) + Y(112)*(gammaX + nuX) + Y(5)*(etah + rho) + Y(3)*(gammah + rho) + Y(35)*(gamma3s + rho) + Y(77)*(nuX + rho);
     dY(8) = Lambdas*Y(1) - Y(8)*(Lambdah + LambdaX + mu + rho + sigmas) + Y(78)*nuX + Y(113)*(gammaX + nuX);
     dY(9) = Lambdah*Y(8) + Lambdas*Y(2) + Y(79)*nuX - Y(9)*(LambdaX + mu + rho + sigmah + sigmas) + Y(114)*(gammaX + nuX);
     dY(10) = Lambdas*Y(3) - Y(10)*(LambdaX + gammah + mu + rho + sigmas) + Y(80)*nuX + Y(9)*sigmah + Y(115)*(gammaX + nuX);
     dY(11) = Lambdas*Y(4) + Y(81)*nuX - Y(11)*(LambdaX + mu + rho + sigmas - Lambdah*(zetah - 1)) + Y(116)*(gammaX + nuX);
     dY(12) = Lambdas*Y(5) + Y(82)*nuX - Y(12)*(LambdaX + etah + mu + rho + sigmah + sigmas) + Y(117)*(gammaX + nuX) - Lambdah*Y(11)*(zetah - 1);
     dY(13) = Lambdas*Y(6) + Y(83)*nuX + Y(12)*sigmah - Y(13)*(LambdaX + etah + gammah + mu + rho + sigmas) + Y(118)*(gammaX + nuX);
     dY(14) = Lambdas*Y(7) + Y(12)*etah + Y(10)*gammah + Y(84)*nuX - Y(14)*(LambdaX + mu + rho + sigmas) + Y(13)*(etah + gammah) + Y(119)*(gammaX + nuX);
     dY(15) = Y(85)*nuX - Y(15)*(Lambdah + LambdaX + mu + rho + taus) + Y(8)*sigmas + Y(120)*(gammaX + nuX);
     dY(16) = Lambdah*Y(15) + Y(86)*nuX + Y(9)*sigmas - Y(16)*(LambdaX + mu + rho + sigmah + taus) + Y(121)*(gammaX + nuX);
     dY(17) = Y(87)*nuX - Y(17)*(LambdaX + gammah + mu + rho + taus) + Y(16)*sigmah + Y(10)*sigmas + Y(122)*(gammaX + nuX);
     dY(18) = Y(88)*nuX + Y(11)*sigmas - Y(18)*(LambdaX + mu + rho + taus - Lambdah*(zetah - 1)) + Y(123)*(gammaX + nuX);
     dY(19) = Y(89)*nuX + Y(12)*sigmas - Y(19)*(LambdaX + etah + mu + rho + sigmah + taus) + Y(124)*(gammaX + nuX) - Lambdah*Y(18)*(zetah - 1);
     dY(20) = Y(90)*nuX + Y(19)*sigmah + Y(13)*sigmas - Y(20)*(LambdaX + etah + gammah + mu + rho + taus) + Y(125)*(gammaX + nuX);
     dY(21) = Y(19)*etah + Y(17)*gammah + Y(91)*nuX + Y(14)*sigmas - Y(21)*(LambdaX + mu + rho + taus) + Y(20)*(etah + gammah) + Y(126)*(gammaX + nuX);
     dY(22) = Y(92)*nuX - Y(22)*(Lambdah + LambdaX + mu + rho + thetas) + Y(15)*taus + Y(127)*(gammaX + nuX);
     dY(23) = Lambdah*Y(22) + Y(93)*nuX - Y(23)*(LambdaX + mu + rho + sigmah + thetas) + Y(16)*taus + Y(128)*(gammaX + nuX);
     dY(24) = Y(94)*nuX - Y(24)*(LambdaX + gammah + mu + rho + thetas) + Y(23)*sigmah + Y(17)*taus + Y(129)*(gammaX + nuX);
     dY(25) = Y(95)*nuX + Y(18)*taus - Y(25)*(LambdaX + mu + rho + thetas - Lambdah*(zetah - 1)) + Y(130)*(gammaX + nuX);
     dY(26) = Y(96)*nuX + Y(19)*taus - Y(26)*(LambdaX + etah + mu + rho + sigmah + thetas) + Y(131)*(gammaX + nuX) - Lambdah*Y(25)*(zetah - 1);
     dY(27) = Y(97)*nuX + Y(26)*sigmah + Y(20)*taus - Y(27)*(LambdaX + etah + gammah + mu + rho + thetas) + Y(132)*(gammaX + nuX);
     dY(28) = Y(26)*etah + Y(24)*gammah + Y(98)*nuX + Y(21)*taus - Y(28)*(LambdaX + mu + rho + thetas) + Y(27)*(etah + gammah) + Y(133)*(gammaX + nuX);
     dY(29) = Y(99)*nuX + Y(22)*thetas + Y(134)*(gammaX + nuX) - Y(29)*(Lambdah + LambdaX + gamma3s + mu + rho);
     dY(30) = Lambdah*Y(29) - Y(30)*(LambdaX + gamma3s + mu + rho + sigmah) + Y(100)*nuX + Y(23)*thetas + Y(135)*(gammaX + nuX);
     dY(31) = Y(101)*nuX - Y(31)*(LambdaX + gammah + gamma3s + mu + rho) + Y(30)*sigmah + Y(24)*thetas + Y(136)*(gammaX + nuX);
     dY(32) = Y(102)*nuX + Y(25)*thetas - Y(32)*(LambdaX + gamma3s + mu + rho - Lambdah*(zetah - 1)) + Y(137)*(gammaX + nuX);
     dY(33) = Y(103)*nuX + Y(26)*thetas - Y(33)*(LambdaX + etah + gamma3s + mu + rho + sigmah) + Y(138)*(gammaX + nuX) - Lambdah*Y(32)*(zetah - 1);
     dY(34) = Y(104)*nuX + Y(33)*sigmah + Y(27)*thetas - Y(34)*(LambdaX + etah + gammah + gamma3s + mu + rho) + Y(139)*(gammaX + nuX);
     dY(35) = Y(33)*etah + Y(31)*gammah + Y(105)*nuX + Y(28)*thetas - Y(35)*(LambdaX + gamma3s + mu + rho) + Y(34)*(etah + gammah) + Y(140)*(gammaX + nuX);
     dY(36) = LambdaX*Y(1) - Y(36)*(Lambdah + Lambdas + mu + rho + sigmaX) + Y(64)*gamma3s;
     dY(37) = Lambdah*Y(36) + LambdaX*Y(2) + Y(65)*gamma3s - Y(37)*(Lambdas + mu + rho + sigmah + sigmaX);
     dY(38) = LambdaX*Y(3) - Y(38)*(Lambdas + gammah + mu + rho + sigmaX) + Y(66)*gamma3s + Y(37)*sigmah;
     dY(39) = LambdaX*Y(4) + Y(67)*gamma3s - Y(39)*(Lambdas + mu + rho + sigmaX - Lambdah*(zetah - 1));
     dY(40) = LambdaX*Y(5) + Y(68)*gamma3s - Y(40)*(Lambdas + etah + mu + rho + sigmah + sigmaX) - Lambdah*Y(39)*(zetah - 1);
     dY(41) = LambdaX*Y(6) + Y(69)*gamma3s + Y(40)*sigmah - Y(41)*(Lambdas + etah + gammah + mu + rho + sigmaX);
     dY(42) = LambdaX*Y(7) + Y(40)*etah + Y(38)*gammah + Y(70)*gamma3s - Y(42)*(Lambdas + mu + rho + sigmaX) + Y(41)*(etah + gammah);
     dY(43) = Lambdas*Y(36) + LambdaX*Y(8) - Y(43)*(Lambdah + mu + rho + sigmas + sigmaX);
     dY(44) = Lambdah*Y(43) + Lambdas*Y(37) + LambdaX*Y(9) - Y(44)*(mu + rho + sigmah + sigmas + sigmaX);
     dY(45) = Lambdas*Y(38) + LambdaX*Y(10) + Y(44)*sigmah - Y(45)*(gammah + mu + rho + sigmas + sigmaX);
     dY(46) = Lambdas*Y(39) + LambdaX*Y(11) - Y(46)*(mu + rho + sigmas + sigmaX - Lambdah*(zetah - 1));
     dY(47) = Lambdas*Y(40) + LambdaX*Y(12) - Y(47)*(etah + mu + rho + sigmah + sigmas + sigmaX) - Lambdah*Y(46)*(zetah - 1);
     dY(48) = Lambdas*Y(41) + LambdaX*Y(13) + Y(47)*sigmah - Y(48)*(etah + gammah + mu + rho + sigmas + sigmaX);
     dY(49) = Lambdas*Y(42) + LambdaX*Y(14) + Y(47)*etah + Y(45)*gammah - Y(49)*(mu + rho + sigmas + sigmaX) + Y(48)*(etah + gammah);
     dY(50) = LambdaX*Y(15) + Y(43)*sigmas - Y(50)*(Lambdah + mu + rho + sigmaX + taus);
     dY(51) = Lambdah*Y(50) + LambdaX*Y(16) + Y(44)*sigmas - Y(51)*(mu + rho + sigmah + sigmaX + taus);
     dY(52) = LambdaX*Y(17) + Y(51)*sigmah + Y(45)*sigmas - Y(52)*(gammah + mu + rho + sigmaX + taus);
     dY(53) = LambdaX*Y(18) + Y(46)*sigmas - Y(53)*(mu + rho + sigmaX + taus - Lambdah*(zetah - 1));
     dY(54) = LambdaX*Y(19) + Y(47)*sigmas - Y(54)*(etah + mu + rho + sigmah + sigmaX + taus) - Lambdah*Y(53)*(zetah - 1);
     dY(55) = LambdaX*Y(20) + Y(54)*sigmah + Y(48)*sigmas - Y(55)*(etah + gammah + mu + rho + sigmaX + taus);
     dY(56) = LambdaX*Y(21) + Y(54)*etah + Y(52)*gammah + Y(49)*sigmas - Y(56)*(mu + rho + sigmaX + taus) + Y(55)*(etah + gammah);
     dY(57) = LambdaX*Y(22) - Y(57)*(Lambdah + mu + rho + sigmaX + thetas) + Y(50)*taus;
     dY(58) = Lambdah*Y(57) + LambdaX*Y(23) + Y(51)*taus - Y(58)*(mu + rho + sigmah + sigmaX + thetas);
     dY(59) = LambdaX*Y(24) + Y(58)*sigmah + Y(52)*taus - Y(59)*(gammah + mu + rho + sigmaX + thetas);
     dY(60) = LambdaX*Y(25) + Y(53)*taus - Y(60)*(mu + rho + sigmaX + thetas - Lambdah*(zetah - 1));
     dY(61) = LambdaX*Y(26) + Y(54)*taus - Y(61)*(etah + mu + rho + sigmah + sigmaX + thetas) - Lambdah*Y(60)*(zetah - 1);
     dY(62) = LambdaX*Y(27) + Y(61)*sigmah + Y(55)*taus - Y(62)*(etah + gammah + mu + rho + sigmaX + thetas);
     dY(63) = LambdaX*Y(28) + Y(61)*etah + Y(59)*gammah + Y(56)*taus - Y(63)*(mu + rho + sigmaX + thetas) + Y(62)*(etah + gammah);
     dY(64) = LambdaX*Y(29) - Y(64)*(Lambdah + gamma3s + mu + rho + sigmaX) + Y(57)*thetas;
     dY(65) = Lambdah*Y(64) + LambdaX*Y(30) + Y(58)*thetas - Y(65)*(gamma3s + mu + rho + sigmah + sigmaX);
     dY(66) = LambdaX*Y(31) + Y(65)*sigmah + Y(59)*thetas - Y(66)*(gammah + gamma3s + mu + rho + sigmaX);
     dY(67) = LambdaX*Y(32) + Y(60)*thetas - Y(67)*(gamma3s + mu + rho + sigmaX - Lambdah*(zetah - 1));
     dY(68) = LambdaX*Y(33) + Y(61)*thetas - Y(68)*(etah + gamma3s + mu + rho + sigmah + sigmaX) - Lambdah*Y(67)*(zetah - 1);
     dY(69) = LambdaX*Y(34) + Y(68)*sigmah + Y(62)*thetas - Y(69)*(etah + gammah + gamma3s + mu + rho + sigmaX);
     dY(70) = LambdaX*Y(35) + Y(68)*etah + Y(66)*gammah + Y(63)*thetas - Y(70)*(gamma3s + mu + rho + sigmaX) + Y(69)*(etah + gammah);
     dY(71) = Y(99)*gamma3s - Y(71)*(Lambdah + Lambdas + mu + nuX + rho) - Y(36)*sigmaX*(epsX - 1);
     dY(72) = Lambdah*Y(71) + Y(100)*gamma3s - Y(72)*(Lambdas + mu + nuX + rho + sigmah) - Y(37)*sigmaX*(epsX - 1);
     dY(73) = Y(101)*gamma3s - Y(73)*(Lambdas + gammah + mu + nuX + rho) + Y(72)*sigmah - Y(38)*sigmaX*(epsX - 1);
     dY(74) = Y(102)*gamma3s - Y(74)*(Lambdas + mu + nuX + rho - Lambdah*(zetah - 1)) - Y(39)*sigmaX*(epsX - 1);
     dY(75) = Y(103)*gamma3s - Y(75)*(Lambdas + etah + mu + nuX + rho + sigmah) - Lambdah*Y(74)*(zetah - 1) - Y(40)*sigmaX*(epsX - 1);
     dY(76) = Y(104)*gamma3s + Y(75)*sigmah - Y(76)*(Lambdas + etah + gammah + mu + nuX + rho) - Y(41)*sigmaX*(epsX - 1);
     dY(77) = Y(75)*etah + Y(73)*gammah + Y(105)*gamma3s - Y(77)*(Lambdas + mu + nuX + rho) + Y(76)*(etah + gammah) - Y(42)*sigmaX*(epsX - 1);
     dY(78) = Lambdas*Y(71) - Y(78)*(Lambdah + mu + nuX + rho + sigmas) - Y(43)*sigmaX*(epsX - 1);
     dY(79) = Lambdah*Y(78) + Lambdas*Y(72) - Y(79)*(mu + nuX + rho + sigmah + sigmas) - Y(44)*sigmaX*(epsX - 1);
     dY(80) = Lambdas*Y(73) + Y(79)*sigmah - Y(80)*(gammah + mu + nuX + rho + sigmas) - Y(45)*sigmaX*(epsX - 1);
     dY(81) = Lambdas*Y(74) - Y(81)*(mu + nuX + rho + sigmas - Lambdah*(zetah - 1)) - Y(46)*sigmaX*(epsX - 1);
     dY(82) = Lambdas*Y(75) - Y(82)*(etah + mu + nuX + rho + sigmah + sigmas) - Lambdah*Y(81)*(zetah - 1) - Y(47)*sigmaX*(epsX - 1);
     dY(83) = Lambdas*Y(76) + Y(82)*sigmah - Y(83)*(etah + gammah + mu + nuX + rho + sigmas) - Y(48)*sigmaX*(epsX - 1);
     dY(84) = Lambdas*Y(77) + Y(82)*etah + Y(80)*gammah - Y(84)*(mu + nuX + rho + sigmas) + Y(83)*(etah + gammah) - Y(49)*sigmaX*(epsX - 1);
     dY(85) = Y(78)*sigmas - Y(85)*(Lambdah + mu + nuX + rho + taus) - Y(50)*sigmaX*(epsX - 1);
     dY(86) = Lambdah*Y(85) + Y(79)*sigmas - Y(86)*(mu + nuX + rho + sigmah + taus) - Y(51)*sigmaX*(epsX - 1);
     dY(87) = Y(86)*sigmah + Y(80)*sigmas - Y(87)*(gammah + mu + nuX + rho + taus) - Y(52)*sigmaX*(epsX - 1);
     dY(88) = Y(81)*sigmas - Y(88)*(mu + nuX + rho + taus - Lambdah*(zetah - 1)) - Y(53)*sigmaX*(epsX - 1);
     dY(89) = Y(82)*sigmas - Y(89)*(etah + mu + nuX + rho + sigmah + taus) - Lambdah*Y(88)*(zetah - 1) - Y(54)*sigmaX*(epsX - 1);
     dY(90) = Y(89)*sigmah + Y(83)*sigmas - Y(90)*(etah + gammah + mu + nuX + rho + taus) - Y(55)*sigmaX*(epsX - 1);
     dY(91) = Y(89)*etah + Y(87)*gammah + Y(84)*sigmas - Y(91)*(mu + nuX + rho + taus) + Y(90)*(etah + gammah) - Y(56)*sigmaX*(epsX - 1);
     dY(92) = Y(85)*taus - Y(92)*(Lambdah + mu + nuX + rho + thetas) - Y(57)*sigmaX*(epsX - 1);
     dY(93) = Lambdah*Y(92) + Y(86)*taus - Y(93)*(mu + nuX + rho + sigmah + thetas) - Y(58)*sigmaX*(epsX - 1);
     dY(94) = Y(93)*sigmah + Y(87)*taus - Y(94)*(gammah + mu + nuX + rho + thetas) - Y(59)*sigmaX*(epsX - 1);
     dY(95) = Y(88)*taus - Y(95)*(mu + nuX + rho + thetas - Lambdah*(zetah - 1)) - Y(60)*sigmaX*(epsX - 1);
     dY(96) = Y(89)*taus - Y(96)*(etah + mu + nuX + rho + sigmah + thetas) - Lambdah*Y(95)*(zetah - 1) - Y(61)*sigmaX*(epsX - 1);
     dY(97) = Y(96)*sigmah + Y(90)*taus - Y(97)*(etah + gammah + mu + nuX + rho + thetas) - Y(62)*sigmaX*(epsX - 1);
     dY(98) = Y(96)*etah + Y(94)*gammah + Y(91)*taus - Y(98)*(mu + nuX + rho + thetas) + Y(97)*(etah + gammah) - Y(63)*sigmaX*(epsX - 1);
     dY(99) = Y(92)*thetas - Y(99)*(Lambdah + gamma3s + mu + nuX + rho) - Y(64)*sigmaX*(epsX - 1);
     dY(100) = Lambdah*Y(99) + Y(93)*thetas - Y(100)*(gamma3s + mu + nuX + rho + sigmah) - Y(65)*sigmaX*(epsX - 1);
     dY(101) = Y(100)*sigmah + Y(94)*thetas - Y(101)*(gammah + gamma3s + mu + nuX + rho) - Y(66)*sigmaX*(epsX - 1);
     dY(102) = Y(95)*thetas - Y(102)*(gamma3s + mu + nuX + rho - Lambdah*(zetah - 1)) - Y(67)*sigmaX*(epsX - 1);
     dY(103) = Y(96)*thetas - Y(103)*(etah + gamma3s + mu + nuX + rho + sigmah) - Lambdah*Y(102)*(zetah - 1) - Y(68)*sigmaX*(epsX - 1);
     dY(104) = Y(103)*sigmah + Y(97)*thetas - Y(104)*(etah + gammah + gamma3s + mu + nuX + rho) - Y(69)*sigmaX*(epsX - 1);
     dY(105) = Y(103)*etah + Y(101)*gammah + Y(98)*thetas - Y(105)*(gamma3s + mu + nuX + rho) + Y(104)*(etah + gammah) - Y(70)*sigmaX*(epsX - 1);
     dY(106) = Y(134)*gamma3s - Y(106)*(Lambdah + Lambdas + gammaX + mu + nuX) + Y(36)*epsX*sigmaX;
     dY(107) = Lambdah*Y(106) - Y(107)*(Lambdas + gammaX + mu + nuX + sigmah) + Y(135)*gamma3s + Y(37)*epsX*sigmaX;
     dY(108) = Y(136)*gamma3s - Y(108)*(Lambdas + gammah + gammaX + mu + nuX) + Y(107)*sigmah + Y(38)*epsX*sigmaX;
     dY(109) = Y(137)*gamma3s - Y(109)*(Lambdas + gammaX + mu + nuX - Lambdah*(zetah - 1)) + Y(39)*epsX*sigmaX;
     dY(110) = Y(138)*gamma3s - Y(110)*(Lambdas + etah + gammaX + mu + nuX + sigmah) + Y(40)*epsX*sigmaX - Lambdah*Y(109)*(zetah - 1);
     dY(111) = Y(139)*gamma3s + Y(110)*sigmah - Y(111)*(Lambdas + etah + gammah + gammaX + mu + nuX) + Y(41)*epsX*sigmaX;
     dY(112) = Y(110)*etah + Y(108)*gammah + Y(140)*gamma3s - Y(112)*(Lambdas + gammaX + mu + nuX) + Y(111)*(etah + gammah) + Y(42)*epsX*sigmaX;
     dY(113) = Lambdas*Y(106) - Y(113)*(Lambdah + gammaX + mu + nuX + sigmas) + Y(43)*epsX*sigmaX;
     dY(114) = Lambdah*Y(113) + Lambdas*Y(107) - Y(114)*(gammaX + mu + nuX + sigmah + sigmas) + Y(44)*epsX*sigmaX;
     dY(115) = Lambdas*Y(108) + Y(114)*sigmah - Y(115)*(gammah + gammaX + mu + nuX + sigmas) + Y(45)*epsX*sigmaX;
     dY(116) = Lambdas*Y(109) - Y(116)*(gammaX + mu + nuX + sigmas - Lambdah*(zetah - 1)) + Y(46)*epsX*sigmaX;
     dY(117) = Lambdas*Y(110) - Y(117)*(etah + gammaX + mu + nuX + sigmah + sigmas) + Y(47)*epsX*sigmaX - Lambdah*Y(116)*(zetah - 1);
     dY(118) = Lambdas*Y(111) + Y(117)*sigmah - Y(118)*(etah + gammah + gammaX + mu + nuX + sigmas) + Y(48)*epsX*sigmaX;
     dY(119) = Lambdas*Y(112) + Y(117)*etah + Y(115)*gammah - Y(119)*(gammaX + mu + nuX + sigmas) + Y(118)*(etah + gammah) + Y(49)*epsX*sigmaX;
     dY(120) = Y(113)*sigmas - Y(120)*(Lambdah + gammaX + mu + nuX + taus) + Y(50)*epsX*sigmaX;
     dY(121) = Lambdah*Y(120) + Y(114)*sigmas - Y(121)*(gammaX + mu + nuX + sigmah + taus) + Y(51)*epsX*sigmaX;
     dY(122) = Y(121)*sigmah + Y(115)*sigmas - Y(122)*(gammah + gammaX + mu + nuX + taus) + Y(52)*epsX*sigmaX;
     dY(123) = Y(116)*sigmas - Y(123)*(gammaX + mu + nuX + taus - Lambdah*(zetah - 1)) + Y(53)*epsX*sigmaX;
     dY(124) = Y(117)*sigmas - Y(124)*(etah + gammaX + mu + nuX + sigmah + taus) + Y(54)*epsX*sigmaX - Lambdah*Y(123)*(zetah - 1);
     dY(125) = Y(124)*sigmah + Y(118)*sigmas - Y(125)*(etah + gammah + gammaX + mu + nuX + taus) + Y(55)*epsX*sigmaX;
     dY(126) = Y(124)*etah + Y(122)*gammah + Y(119)*sigmas - Y(126)*(gammaX + mu + nuX + taus) + Y(125)*(etah + gammah) + Y(56)*epsX*sigmaX;
     dY(127) = Y(120)*taus - Y(127)*(Lambdah + gammaX + mu + nuX + thetas) + Y(57)*epsX*sigmaX;
     dY(128) = Lambdah*Y(127) + Y(121)*taus - Y(128)*(gammaX + mu + nuX + sigmah + thetas) + Y(58)*epsX*sigmaX;
     dY(129) = Y(128)*sigmah + Y(122)*taus - Y(129)*(gammah + gammaX + mu + nuX + thetas) + Y(59)*epsX*sigmaX;
     dY(130) = Y(123)*taus - Y(130)*(gammaX + mu + nuX + thetas - Lambdah*(zetah - 1)) + Y(60)*epsX*sigmaX;
     dY(131) = Y(124)*taus - Y(131)*(etah + gammaX + mu + nuX + sigmah + thetas) + Y(61)*epsX*sigmaX - Lambdah*Y(130)*(zetah - 1);
     dY(132) = Y(131)*sigmah + Y(125)*taus - Y(132)*(etah + gammah + gammaX + mu + nuX + thetas) + Y(62)*epsX*sigmaX;
     dY(133) = Y(131)*etah + Y(129)*gammah + Y(126)*taus - Y(133)*(gammaX + mu + nuX + thetas) + Y(132)*(etah + gammah) + Y(63)*epsX*sigmaX;
     dY(134) = Y(127)*thetas - Y(134)*(Lambdah + gamma3s + gammaX + mu + nuX) + Y(64)*epsX*sigmaX;
     dY(135) = Lambdah*Y(134) + Y(128)*thetas - Y(135)*(gamma3s + gammaX + mu + nuX + sigmah) + Y(65)*epsX*sigmaX;
     dY(136) = Y(135)*sigmah - Y(136)*(gammah + gamma3s + gammaX + mu + nuX) + Y(129)*thetas + Y(66)*epsX*sigmaX;
     dY(137) = Y(130)*thetas - Y(137)*(gamma3s + gammaX + mu + nuX - Lambdah*(zetah - 1)) + Y(67)*epsX*sigmaX;
     dY(138) = Y(131)*thetas - Y(138)*(etah + gamma3s + gammaX + mu + nuX + sigmah) + Y(68)*epsX*sigmaX - Lambdah*Y(137)*(zetah - 1);
     dY(139) = Y(138)*sigmah + Y(132)*thetas - Y(139)*(etah + gammah + gamma3s + gammaX + mu + nuX) + Y(69)*epsX*sigmaX;
     dY(140) = Y(138)*etah + Y(136)*gammah + Y(133)*thetas - Y(140)*(gamma3s + gammaX + mu + nuX) + Y(139)*(etah + gammah) + Y(70)*epsX*sigmaX;

    %          HIV1    syph1    STI1    no         X    
    %     ____    _____    ____    ___    _________
    % 
    %     "S"     "S"      "S"       1    [1×1 sym]
    %     "I"     "S"      "S"       2    [1×1 sym]
    %     "C"     "S"      "S"       3    [1×1 sym]
    %     "P"     "S"      "S"       4    [1×1 sym]
    %     "Ip"    "S"      "S"       5    [1×1 sym]
    %     "Cp"    "S"      "S"       6    [1×1 sym]
    %     "T"     "S"      "S"       7    [1×1 sym]
    %     "S"     "E"      "S"       8    [1×1 sym]
    %     "I"     "E"      "S"       9    [1×1 sym]
    %     "C"     "E"      "S"      10    [1×1 sym]
    %     "P"     "E"      "S"      11    [1×1 sym]
    %     "Ip"    "E"      "S"      12    [1×1 sym]
    %     "Cp"    "E"      "S"      13    [1×1 sym]
    %     "T"     "E"      "S"      14    [1×1 sym]
    %     "S"     "I1"     "S"      15    [1×1 sym]
    %     "I"     "I1"     "S"      16    [1×1 sym]
    %     "C"     "I1"     "S"      17    [1×1 sym]
    %     "P"     "I1"     "S"      18    [1×1 sym]
    %     "Ip"    "I1"     "S"      19    [1×1 sym]
    %     "Cp"    "I1"     "S"      20    [1×1 sym]
    %     "T"     "I1"     "S"      21    [1×1 sym]
    %     "S"     "I2"     "S"      22    [1×1 sym]
    %     "I"     "I2"     "S"      23    [1×1 sym]
    %     "C"     "I2"     "S"      24    [1×1 sym]
    %     "P"     "I2"     "S"      25    [1×1 sym]
    %     "Ip"    "I2"     "S"      26    [1×1 sym]
    %     "Cp"    "I2"     "S"      27    [1×1 sym]
    %     "T"     "I2"     "S"      28    [1×1 sym]
    %     "S"     "I3"     "S"      29    [1×1 sym]
    %     "I"     "I3"     "S"      30    [1×1 sym]
    %     "C"     "I3"     "S"      31    [1×1 sym]
    %     "P"     "I3"     "S"      32    [1×1 sym]
    %     "Ip"    "I3"     "S"      33    [1×1 sym]
    %     "Cp"    "I3"     "S"      34    [1×1 sym]
    %     "T"     "I3"     "S"      35    [1×1 sym]
    %     "S"     "S"      "E"      36    [1×1 sym]
    %     "I"     "S"      "E"      37    [1×1 sym]
    %     "C"     "S"      "E"      38    [1×1 sym]
    %     "P"     "S"      "E"      39    [1×1 sym]
    %     "Ip"    "S"      "E"      40    [1×1 sym]
    %     "Cp"    "S"      "E"      41    [1×1 sym]
    %     "T"     "S"      "E"      42    [1×1 sym]
    %     "S"     "E"      "E"      43    [1×1 sym]
    %     "I"     "E"      "E"      44    [1×1 sym]
    %     "C"     "E"      "E"      45    [1×1 sym]
    %     "P"     "E"      "E"      46    [1×1 sym]
    %     "Ip"    "E"      "E"      47    [1×1 sym]
    %     "Cp"    "E"      "E"      48    [1×1 sym]
    %     "T"     "E"      "E"      49    [1×1 sym]
    %     "S"     "I1"     "E"      50    [1×1 sym]
    %     "I"     "I1"     "E"      51    [1×1 sym]
    %     "C"     "I1"     "E"      52    [1×1 sym]
    %     "P"     "I1"     "E"      53    [1×1 sym]
    %     "Ip"    "I1"     "E"      54    [1×1 sym]
    %     "Cp"    "I1"     "E"      55    [1×1 sym]
    %     "T"     "I1"     "E"      56    [1×1 sym]
    %     "S"     "I2"     "E"      57    [1×1 sym]
    %     "I"     "I2"     "E"      58    [1×1 sym]
    %     "C"     "I2"     "E"      59    [1×1 sym]
    %     "P"     "I2"     "E"      60    [1×1 sym]
    %     "Ip"    "I2"     "E"      61    [1×1 sym]
    %     "Cp"    "I2"     "E"      62    [1×1 sym]
    %     "T"     "I2"     "E"      63    [1×1 sym]
    %     "S"     "I3"     "E"      64    [1×1 sym]
    %     "I"     "I3"     "E"      65    [1×1 sym]
    %     "C"     "I3"     "E"      66    [1×1 sym]
    %     "P"     "I3"     "E"      67    [1×1 sym]
    %     "Ip"    "I3"     "E"      68    [1×1 sym]
    %     "Cp"    "I3"     "E"      69    [1×1 sym]
    %     "T"     "I3"     "E"      70    [1×1 sym]
    %     "S"     "S"      "IA"     71    [1×1 sym]
    %     "I"     "S"      "IA"     72    [1×1 sym]
    %     "C"     "S"      "IA"     73    [1×1 sym]
    %     "P"     "S"      "IA"     74    [1×1 sym]
    %     "Ip"    "S"      "IA"     75    [1×1 sym]
    %     "Cp"    "S"      "IA"     76    [1×1 sym]
    %     "T"     "S"      "IA"     77    [1×1 sym]
    %     "S"     "E"      "IA"     78    [1×1 sym]
    %     "I"     "E"      "IA"     79    [1×1 sym]
    %     "C"     "E"      "IA"     80    [1×1 sym]
    %     "P"     "E"      "IA"     81    [1×1 sym]
    %     "Ip"    "E"      "IA"     82    [1×1 sym]
    %     "Cp"    "E"      "IA"     83    [1×1 sym]
    %     "T"     "E"      "IA"     84    [1×1 sym]
    %     "S"     "I1"     "IA"     85    [1×1 sym]
    %     "I"     "I1"     "IA"     86    [1×1 sym]
    %     "C"     "I1"     "IA"     87    [1×1 sym]
    %     "P"     "I1"     "IA"     88    [1×1 sym]
    %     "Ip"    "I1"     "IA"     89    [1×1 sym]
    %     "Cp"    "I1"     "IA"     90    [1×1 sym]
    %     "T"     "I1"     "IA"     91    [1×1 sym]
    %     "S"     "I2"     "IA"     92    [1×1 sym]
    %     "I"     "I2"     "IA"     93    [1×1 sym]
    %     "C"     "I2"     "IA"     94    [1×1 sym]
    %     "P"     "I2"     "IA"     95    [1×1 sym]
    %     "Ip"    "I2"     "IA"     96    [1×1 sym]
    %     "Cp"    "I2"     "IA"     97    [1×1 sym]
    %     "T"     "I2"     "IA"     98    [1×1 sym]
    %     "S"     "I3"     "IA"     99    [1×1 sym]
    %     "I"     "I3"     "IA"    100    [1×1 sym]
    %     "C"     "I3"     "IA"    101    [1×1 sym]
    %     "P"     "I3"     "IA"    102    [1×1 sym]
    %     "Ip"    "I3"     "IA"    103    [1×1 sym]
    %     "Cp"    "I3"     "IA"    104    [1×1 sym]
    %     "T"     "I3"     "IA"    105    [1×1 sym]
    %     "S"     "S"      "IS"    106    [1×1 sym]
    %     "I"     "S"      "IS"    107    [1×1 sym]
    %     "C"     "S"      "IS"    108    [1×1 sym]
    %     "P"     "S"      "IS"    109    [1×1 sym]
    %     "Ip"    "S"      "IS"    110    [1×1 sym]
    %     "Cp"    "S"      "IS"    111    [1×1 sym]
    %     "T"     "S"      "IS"    112    [1×1 sym]
    %     "S"     "E"      "IS"    113    [1×1 sym]
    %     "I"     "E"      "IS"    114    [1×1 sym]
    %     "C"     "E"      "IS"    115    [1×1 sym]
    %     "P"     "E"      "IS"    116    [1×1 sym]
    %     "Ip"    "E"      "IS"    117    [1×1 sym]
    %     "Cp"    "E"      "IS"    118    [1×1 sym]
    %     "T"     "E"      "IS"    119    [1×1 sym]
    %     "S"     "I1"     "IS"    120    [1×1 sym]
    %     "I"     "I1"     "IS"    121    [1×1 sym]
    %     "C"     "I1"     "IS"    122    [1×1 sym]
    %     "P"     "I1"     "IS"    123    [1×1 sym]
    %     "Ip"    "I1"     "IS"    124    [1×1 sym]
    %     "Cp"    "I1"     "IS"    125    [1×1 sym]
    %     "T"     "I1"     "IS"    126    [1×1 sym]
    %     "S"     "I2"     "IS"    127    [1×1 sym]
    %     "I"     "I2"     "IS"    128    [1×1 sym]
    %     "C"     "I2"     "IS"    129    [1×1 sym]
    %     "P"     "I2"     "IS"    130    [1×1 sym]
    %     "Ip"    "I2"     "IS"    131    [1×1 sym]
    %     "Cp"    "I2"     "IS"    132    [1×1 sym]
    %     "T"     "I2"     "IS"    133    [1×1 sym]
    %     "S"     "I3"     "IS"    134    [1×1 sym]
    %     "I"     "I3"     "IS"    135    [1×1 sym]
    %     "C"     "I3"     "IS"    136    [1×1 sym]
    %     "P"     "I3"     "IS"    137    [1×1 sym]
    %     "Ip"    "I3"     "IS"    138    [1×1 sym]
    %     "Cp"    "I3"     "IS"    139    [1×1 sym]
    %     "T"     "I3"     "IS"    140    [1×1 sym]
end

