backUpFolder = 'Runtest7';
roundFolder = '_round_1';
path2 = ['C:/Users/Moi/Documents/IPLESP/These/Codes/Simulations_StrategiesPC/',backUpFolder,'/',roundFolder,'/'];

paramCt = readtable([path2,'allParametersSets_Ct.txt']);
paramNg = readtable([path2,'allParametersSets_Ng.txt']);
paramHIV = readtable([path2,'allParametersSets_HIV.txt']);
paramS = readtable([path2,'allParametersSets_syphilis.txt']);

lastIter = paramCt(end,:).IDech;

paramTab{1} = table2struct(paramCt(paramCt.IDech==lastIter,:));
paramTab{2} = table2struct(paramNg(paramNg.IDech==lastIter,:));
paramTab{3} = table2struct(paramHIV(paramHIV.IDech==lastIter,:));
paramTab{4} = table2struct(paramS(paramS.IDech==lastIter,:));

pHIV = paramTab{3}.p;
mu = paramTab{3}.mu;

