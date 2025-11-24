function [dY] = ODE_SICRSEIIISSEIIS(t,Y,M,betaI,betaC,betaSyph,betaSEIIS,indexIHIV,indexCHIV,indexInfSyph,...
                                    indexInfSEIIS,lambdaHIV,lambdaSyph,Lambda1,b)
    N = sum(Y);
    indexInf1 = indexInfSEIIS(1,:);    
    totInf1 = sum(Y(indexInf1));
    lambda1Num   = betaSEIIS(1)*totInf1/N;
    totInfIHIV = sum(Y(indexIHIV));
    totInfCHIV = sum(Y(indexCHIV));
    lambdaHIVNum = (betaI*totInfIHIV + betaC*totInfCHIV)/N;
    totInfSyph = sum(Y(indexInfSyph));
    lambdaSyphNum = betaSyph*totInfSyph/N;
    
    M = subs(M,lambdaHIV,lambdaHIVNum);
    M = subs(M,Lambda1,lambda1Num);
    M = subs(M,lambdaSyph,lambdaSyphNum);
    
    dY = double(M)*Y; dY(1) = dY(1)+b;
    
    t
    sum(Y)
end

