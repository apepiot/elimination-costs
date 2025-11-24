function [U12,dU,P] = U12_SICTPSEIIS_v7(param1,param2,mu,b,paramRho,vecRho,c,f,solveMethod,optSolver)
    optSolver.tolP0=1e-5;
    [P,~,~,~,~] = P12_SICTPSEIIS_v7(param1,param2,mu,b,paramRho,vecRho,f,solveMethod,optSolver);

    U12 = vecRho*(P-c)*(vecRho>=0);
    dU = 0;


    
%     HIV1    STI1    no
%     ____    ____    __
% 
%     "S"     "S"      1
%     "I"     "S"      2
%     "C"     "S"      3
%     "P"     "S"      4
%     "Ip"    "S"      5
%     "Cp"    "S"      6
%     "T"     "S"      7
%     "S"     "E"      8
%     "I"     "E"      9
%     "C"     "E"     10
%     "P"     "E"     11
%     "Ip"    "E"     12
%     "Cp"    "E"     13
%     "T"     "E"     14
%     "S"     "IA"    15
%     "I"     "IA"    16
%     "C"     "IA"    17
%     "P"     "IA"    18
%     "Ip"    "IA"    19
%     "Cp"    "IA"    20
%     "T"     "IA"    21
%     "S"     "IS"    22
%     "I"     "IS"    23
%     "C"     "IS"    24
%     "P"     "IS"    25
%     "Ip"    "IS"    26
%     "Cp"    "IS"    27
%     "T"     "IS"    28
end