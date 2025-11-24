function [res] = findRhohat_v4_min(N,param2,mu,b,c,f)
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
    tab = findRhohat_v4(N,param2,mu,b,c,f);
    res = tab.rhohat.rhohat;
end

