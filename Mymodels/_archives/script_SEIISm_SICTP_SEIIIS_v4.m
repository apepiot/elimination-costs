%% This code generates the matrix and the ODE system associated to a SICR^n x SEIIIS^p x SEIIS^m model
% combined and targeted testing are mixed
% from v_3 : correcting the algot that assigns voluntary tetsing rates
clear all;

%% Parameter initialization
% Nothing to change below, this code generates a 4-disease model.
n = 1; %number of SICTPrEP (min=?, max=1)
p = 1; %number of SEIIIS (min=?,max=1)
m = 2; %number of SEIIS (min=?, max m=4)

[nbCompartments,M,B] = createODEsystem_v4(n,p,m)

%% Converting the matrix product to ODE system
[X,dX,eqn,F] = matToODE_v2(nbCompartments,M,B);
eqn.'
tabComp.X = X;

%% Write the ODE system in a text file
pathW = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Mymodels\';
fileID = fopen([pathW,'ODE_SICTPSEIIISSEIIS2_vNEW.txt'],'w');
for k=1:size(tabComp,1)
    fprintf(fileID,'%12s%1s\r\n ',eqn(k),';');    
end
fclose(fileID);




if(0)







%% Computation of the gradient
%syms dXkdXi [nbCompartments nbCompartments]

dXkdXi = []%zeros(nbCompartments,nbCompartments);

clear all;
pathW = 'C:\Users\Moi\Documents\IPLESP\These\Codes\Mymodels\';
fileID = fopen([pathW,'gradF_ODE_SICTPSEIIISSEIIS_v4_2.txt'],'w');

F_ODESICTPSEIIISSEIIS2_v4_detailed
%syms dXkdXi [1 1]
for k=1:560
    for i=1:560
        %dXkdXi(k,i)=diff(F(k),X(i));
        res=diff(F(k),X(i));
        fprintf(fileID,'%10s%100s%1s\r\n ',['dF(',num2str(k),',',num2str(i),') = '],res,';'); 
    end
    disp([k,i]);
end
fclose(fileID);


%% Test du fichier du gradient:
tic
for i=1:10
grad = gradF_ODE_SICTPSEIIISSEIIS2_v4_fun(NaN,rand(560,1),paramTab{3}.betaI,paramTab{3}.betaC,paramTab{3}.sigma,paramTab{3}.theta,paramTab{3}.zeta,paramTab{3}.eta,paramTab{3}.p,...
    paramTab{4}.beta,paramTab{4}.sigma,paramTab{4}.gamma3,paramTab{4}.tau,paramTab{4}.theta,...
    0,paramTab{1}.gamma,paramTab{1}.nu,paramTab{1}.eps,paramTab{1}.sigma,...
    0,paramTab{2}.gamma,paramTab{2}.nu,paramTab{2}.eps,paramTab{2}.sigma,...
    mu,b,...
    0,0,0,0,...
    0,0,0,0,0,0.1,0,0,0,0,0,0,0,0,0,0,0);
disp(i);
toc
end
toc


end


