function [R1, R2, alpha1, alpha2, deltaMaxth1,deltaMaxth2,maxU,deltaUmax,max1,delta1Num,max2,delta2Num,max12,delta12Num] = findMaxU(beta1,beta2,gamma1,gamma2,vecDelta,s1,s2,b,mu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    R1 = beta1/(gamma1+mu);
    R2 = beta2/(gamma2+mu);
    alpha1 = beta1/s1*(1-1/R1);
    alpha2 = beta2/s2*(1-1/R2);
    
    deltaMaxth1 = beta1/s1*(1/sqrt(R1)-1/R1);
    deltaMaxth2 = beta2/(2*s2)*(1-1/R2);
    
    %Numerically
    [U,vecPrev, vecR10,vecR20,P1,P2,P12] = utilityFunctionSIRxSIS(beta1,beta2,gamma1,gamma2,vecDelta,vecDelta,s1,s2,b,mu);
    
    [maxU, iUmax]= max(U);
    deltaUmax = vecDelta(iUmax);
    
    [max1, i1max]= max(vecDelta.*P1);
    delta1Num = vecDelta(i1max);
    [max2, i2max] = max(vecDelta.*P2);
    delta2Num = vecDelta(i2max);
    [max12, i12max] = max(vecDelta.*P12);
    delta12Num = vecDelta(i12max);
    
    
    
end

