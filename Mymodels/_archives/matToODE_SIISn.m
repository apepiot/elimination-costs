function [C,dC] = matToODE_SIISn(n,M)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    nbBoxes = 3;
    
     nbCompartments = nbBoxes^n;
%     x = 1:nbBoxes;%'SIJ';                 %// Set of possible letters                     %// Length of each permutation
%     C = cell(n, 1);             %// Preallocate a cell arracompartments
%     [C{:}] = ndgrid(x);         %// Create K grids of values
%     compartments = cellfun(@(x){x(:)}, C); %// Convert grids to column vectors
%     compartments = [compartments{:}];
%     Month = num2str(compartments);
%     %Month = Month(find(~isspace(Month)))
%     %Month(Month == ' ') = [];
% 
%     monthsArray = Month ;%strsplit(Month,',');

    C = sym('I',[1 nbCompartments]);
    dC = sym('dI',[1 nbCompartments]);
    syms b;

    for i = 1:nbCompartments
%         %Varnames{i} = matlab.lang.makeValidName(strcat('I',monthsArray(i,:)));
         %myStruct.(Varnames{i}) = M(i,:)*eval(Varnames);
         dC(i) = sum(M(i,:).*C);
     end
    %myStruct.(Varnames{1,1}) % should give you a value of a random number
    %myStruct.C01 % same result above    
    dC(1) = dC(1)+b;
end