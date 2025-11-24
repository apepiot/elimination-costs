function param = getParams(k, paramTab)

for i=1:length(paramTab)
    if strcmp(paramTab{i}.mini_d,k) || strcmp(paramTab{i}.disease,k)
        param = paramTab{i};
    end 
end

end