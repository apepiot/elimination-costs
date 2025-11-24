function [new] = transformCharTex(chaine)
new=[];
for k=1:length(chaine)   
    if strcmp(chaine(k),'_')
        new = [new,'\',char(chaine(k))];
    else
        new = [new,char(chaine(k))];
    end
end
end

