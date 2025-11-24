function [z,p] = myStudentTest(ech1,ech2)
    %t-test
    n1 = length(ech1);
    n2 = length(ech2);
    m1 = mean(ech1);
    m2 = mean(ech2);
    v1 = sum((ech1 - m1).^2)./(n1-1); %s1²
    v2 = sum((ech2 - m2).^2)./(n2-1); %s2²
    
    z = (m1-m2)./(sqrt(v1/(n1-1))+sqrt(v2/(n2-1)));
    p = normcdf(-abs(z),0,1)*2;
    %CI_low = 
end

