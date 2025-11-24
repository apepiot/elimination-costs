function tableToTex_col_row(tab,pathW,caption,label,longtable)

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
if longtable
    begin = {['\begin{longtable}[h]{','|p{0.8cm}|p{1.5cm}',repmat('|p{1cm}',1,m-3),'|p{1.5cm}|}'],['\caption{',caption,'}'], ['\label{',label,'} \\'],['\hline']};
else
    begin = {'\begin{table}[h]', '\centering ','\footnotesize', ['\caption{',caption,'}'], ['\label{',label,'}'],'\rowcolors{1}{}{lightgray}',...
        ['\begin{tabular}[h]{|',repmat(' l',1,m),'|'],'} \hline'};
end

for k=1:length(begin)
    fprintf(fileID,'%12s\r\n',[begin{k}]);
end

% Variable Names
varLine = ['\textbf{ ',transformCharTex(varNames{1}),'}'];
for k=2:size(tab,2)
    var_k = transformCharTex(varNames{k});
    varLine = [varLine,' & \textbf{ ',var_k,'}'];
end
fprintf(fileID,'%12s\r\n ',varLine,'\\');

% Line by line
for j=1:n
    if firstColumBold
        LineJ = [' \textbf{',transformCharTex(tab(j,1).(varNames{1})),'}'];
    else
        LineJ = [' ',transformCharTex(tab(j,1).(varNames{1})),''];
    end
    for k=2:m
        res = char(tab(j,k).(varNames{k}));
        if contains(res,'NaN') && replaceNaN
            LineJ = [LineJ,' & ','-'];
        else
            LineJ = [LineJ,' & ',tab(j,k).(varNames{k})];
        end
    end
     %if j<n
         fprintf(fileID,'%12s\r\n ',LineJ,'\\');
%      else
%          fprintf(fileID,'%12s\r\n ',LineJ,' ');
%      end
end

if longtable
    fprintf(fileID,'%12s\r\n','\end{longtable}');
else
    fprintf(fileID,'%15s\r\n','\hline \end{tabular}');
    fprintf(fileID,'%12s\r\n','\end{table}');
end
fclose(fileID);

end

