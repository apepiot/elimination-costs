function tableToTex_damier(tab,pathW,caption,label)

replaceNaN = 1;
firstColumBold = 1;

fileID = fopen(pathW,'w');

[n,m] = size(tab);
varNames = tab.Properties.VariableNames;

% Changing data into char
for k=1:m
    varType = class(tab.(varNames{k}));    
    if strcmp(varType,'cell')
       tab.(varNames{k}) = char(tab.(varNames{k}));
    end   
    if strcmp(varType,'double')
        tab.(varNames{k}) = num2str(round((tab.(varNames{k})),3));
    end
end

%Initialization
begin = {'\begin{table}[h]', '\centering ','\footnotesize', ['\caption{',caption,'}'], ['\label{',label,'}'],...
    ['\begin{tabular}[h]{|l ',repmat('c',1,m-1),'|'],'} \hline'};


for k=1:length(begin)
    fprintf(fileID,'%12s\r\n',[begin{k}]);
end

% Variable Names
%varLine = ['\texttt{',strtrim(transformCharTex(varNames{1})),'}'];
varLine = ['\texttt{',strtrim(varNames{1}),'}'];
for k=2:size(tab,2)
    %var_k = transformCharTex(varNames{k});
    var_k = varNames{k};
    %varLine = [varLine,' & $\pmb{',var_k,'}$'];
    varLine = [varLine,' & $',var_k,'$'];
end
fprintf(fileID,'%12s\r\n ',varLine,'\\ ');

% Line by line
for j=1:n
    if firstColumBold
        LineJ = [' ${',strtrim(transformCharTex(tab(j,1).(varNames{1}))),'}$'];
    else
        LineJ = [' ',transformCharTex(tab(j,1).(varNames{1})),''];
    end
    for k=2:m
        if ismissing(tab(j,k).(varNames{k}))
            res = ' ';
            tab(j,k).(varNames{k}) = " ";
        end
            res = char(tab(j,k).(varNames{k}));
        
        if contains(res,'NaN') && replaceNaN
            LineJ = [LineJ,' & ',' '];
        else
            if mod(k,2) 
                if mod(j,2)
                    LineJ = [LineJ,' & ','\cellcolor{lightgray3}{',strtrim(tab(j,k).(varNames{k})),'}'];
                else
                    LineJ = [LineJ,' & ','\cellcolor{lightgray}{',strtrim(tab(j,k).(varNames{k})),'}'];
                end
            else
                if mod(j,2)
                    LineJ = [LineJ,' & ','\cellcolor{lightgray2}{',strtrim(tab(j,k).(varNames{k})),'}'];
                else
                    LineJ = [LineJ,' & ',strtrim(tab(j,k).(varNames{k}))];
                end
            end
        end
    end
    fprintf(fileID,'%12s\r\n ',join(LineJ),'\\ ');
end

fprintf(fileID,'%12s\r\n','\hline \end{tabular}');
fprintf(fileID,'%12s\r\n','\end{table}');

fclose(fileID);

end

