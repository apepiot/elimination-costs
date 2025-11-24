function [kit] = kToKit(k)
kit={};
j=1;
if contains(k,'h')
    kit{j} ='HIV';
    j=j+1;
end
if contains(k,'s')
    kit{j} ='syphilis';
    j=j+1;
end
if contains(k,'c')
    kit{j} ='Ct';
    j=j+1;
end
if contains(k,'g')
    kit{j} ='Ng';
    j=j+1;
end

end

