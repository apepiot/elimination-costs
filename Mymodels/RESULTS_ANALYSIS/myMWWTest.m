function [Z,p] = myMWWTest(ech1,ech2)
%Mann-Whithey-Wilcon test
% distribution comparison
n1 = length(ech1);
n2 = length(ech2);
if n1<=n2
    echA=ech1;nA=n1;
    echB=ech2;nB=n2;
else
    echA=ech2;nA=n2;
    echB=ech1;nB=n1;
end

echtot = table('Size',[nA+nB,2],'VariableTypes',["double","string"],...
    'VariableNames',["data","sample"]);
echtot{1:nA,"data"} = echA;
echtot{1:nA,"sample"} = "A";
echtot{(nA+1):(nA+nB),"data"} = echB;
echtot{(nA+1):(nA+nB),"sample"} = "B";

echtot_sorted           = sortrows(echtot,"data");
echtot_sorted.data      = round(echtot_sorted.data,4);
[~,echtot_sorted.range] = sort(echtot_sorted.data);

%Rangs au sens du test : faire la moyenne des rangs si egalité
dataUnique = unique(echtot_sorted.data);
nUnique = length(dataUnique);
for i=1:nUnique
    thisValueRanges = echtot_sorted.range(echtot_sorted.data==dataUnique(i));
    echtot_sorted.rangeMWW(echtot_sorted.data==dataUnique(i)) = mean(thisValueRanges);
end

TA    = sum(echtot_sorted.rangeMWW(echtot_sorted.sample=="A"));
delta = TA - nA*(nA+nB+1)./2;
if delta>0
    TAp = TA-0.5;
else
    TAp = TA+0.5;
end

if nA<=10 && nB<=10
    Z = (TA-nA*(nA+nB+1)/2)./sqrt(nA*nB*(nA+nB+1)/12);
    % p a voir (cf table de wmm)
    p=NaN;
else
    Z = (TAp-nA*(nA+nB+1)/2)./sqrt(nA*nB*(nA+nB+1)/12); %follows a N(0,1)
    p = normcdf(-abs(Z),0,1)*2;
end

end

