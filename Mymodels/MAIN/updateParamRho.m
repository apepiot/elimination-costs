function [paramTab,paramRho] = updateParamRho(paramTab,paramRho,kit,newRho)

if abs(newRho)<1e-7
    newRho=0;
end

k=[]; 
for j=[3,4,1,2]
    if ismember(paramTab{j}.disease,kit)
        k=[k,paramTab{j}.mini_d];
    end
end

if isequal(kit,{'HIV'})
    paramTab{3}.rhob = newRho;
end
if isequal(kit,{'syphilis'})
    paramTab{4}.rhob = newRho;
end
if isequal(kit,{'Ct'})
    paramTab{1}.rhob = newRho;
end
if isequal(kit,{'Ng'})
    paramTab{2}.rhob = newRho;
end

paramRho.(['rho_',k]) = newRho;
end

