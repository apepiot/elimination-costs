function [newTab] = concatenateTabRhohat(oldTab, addedTab)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
newTab = oldTab;
oldC = oldTab.c;
newC = addedTab.c;
combinedC = [oldC,newC];
[newTab.c, sortedC] = sort(combinedC);

N = length(oldTab.one);

%concatenate 1 disease model
for k=1:N
    tempOne = [oldTab.one(k).rhohat;addedTab.one(k).rhohat];
    newTab.one(k).rhohat = tempOne(sortedC,:);
end

%concatenate 2 disease model
if N>=2
    for k=1:nchoosek(N,2) %N=1->0, N=2->1, N=3->3, N=4->6
        tempTwo = [oldTab.two(k).rhohat, addedTab.two(k).rhohat];
        newTab.two(k).rhohat = tempTwo(sortedC);
        tempNM = [oldTab.nm(k).rhohat, addedTab.nm(k).rhohat];
        newTab.nm(k).rhohat = tempNM(sortedC);
    end
end

%concatenate 3 disease model
if N>=3
    for k=1:nchoosek(N,3)
        tempThree = [oldTab.three(k).rhohat, addedTab.three(k).rhohat];
        newTab.three(k).rhohat = tempThree(sortedC);
        tempNMK = [oldTab.nmk(k).rhohat, addedTab.nmk(k).rhohat];
        newTab.nmk(k).rhohat = tempNMK(sortedC);
    end
end

%concatenate 4 disease model
if N==4
    tempFour = [oldTab.four.rhohat, addedTab.four.rhohat];
    newTab.four.rhohat = tempFour(sortedC);
    tempNMKL = [oldTab.nmkl; addedTab.nmkl];
    newTab.nmkl = tempNMKL(sortedC);
end

tempRHOHAT = [oldTab.rhohat.rhohat, addedTab.rhohat.rhohat];
newTab.rhohat.rhohat = tempRHOHAT(sortedC);
end

