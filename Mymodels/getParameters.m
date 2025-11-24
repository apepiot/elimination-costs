
% Lecture des paramètres
pathParam = ['.\ParameterAnalysis\paramSets_',num2str(paramNo),'\round_',num2str(roundNo),'\'];
paramH = readtable([pathParam, 'allParametersSets_HIV.txt']);
paramS = readtable([pathParam, 'allParametersSets_syphilis.txt']);
paramC = readtable([pathParam, 'allParametersSets_Ct.txt']);
paramG = readtable([pathParam, 'allParametersSets_Ng.txt']);
ID_ech = paramH(paramH.p==pHIV & paramH.nbEch==nbEch,:).IDech_id;
paramTab{1} = table2struct(paramC(paramC.IDech_id==ID_ech,:));
paramTab{2} = table2struct(paramG(paramG.IDech_id==ID_ech,:));
paramTab{3} = table2struct(paramH(paramH.IDech_id==ID_ech & paramH.p==pHIV,:));
paramTab{4} = table2struct(paramS(paramS.IDech_id==ID_ech,:));
mu=paramTab{3}.mu; b = paramTab{3}.pi;

verbose=0;
paramSolver.tolP0 = 0.5e-4;
paramSolver.maxBndAlpha=20;
paramSolver.nbRelanceMax=5;
paramSolver.timeLimit = 20; %seconds
paramSolver.iterMaxDicho = 20;
paramSolver.tolAlpha = 1e-4;
paramSolver.method_alpha = 'dicho';
paramSolver.timeSolver = 20;
%opt.TolP0=0.5e-4;

pathRes = ['.\ParameterAnalysis\results_',num2str(paramNo),'\_round_',num2str(roundNo),'\'];
paramRho = table2struct(readtable([pathRes ,'paramRho.txt']));