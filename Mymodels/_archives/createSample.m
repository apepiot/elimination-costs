function [output] = createSample(minValue,meanValue,maxValue,n,type)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    if strcmp(type,"uniformly")
        output = linspace(minValue,maxValue,n);
    end
end