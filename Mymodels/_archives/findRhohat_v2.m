function [tabRhohat,tabC0,tabCnn,tabTimes] = findRhohat_v2(nSEIIS,nSICR,nSEIIIS,paramSEIIS,paramSICR,paramSEIIIS,mu,b,vecC,f)
%for any combination of SIS^nSIS x SIR^nSIR x SICAT^nSICAT model, this
%function finds the argmax of single diseases models, two-disease
%models,...,N diseases models
tic;middlevecC = median(vecC);
%checks
if nSEIIS~=size(paramSEIIS,1)
   error('Error. The number of SEIIS diseases does not match with the parameters')
end
if nSICR~=size(paramSICR,1)
   error('Error. The number of SICR diseases does not match with the parameters')
end
if nSEIIIS~=size(paramSEIIIS,1)
   error('Error. The number of SEIIIS diseases does not match with the parameters')
end

N = nSEIIS + nSICR + nSEIIIS; %total number of diseases in the model
modelTypes = repelem(["SEIIS","SICR", "SEIIIS"], [nSEIIS nSICR nSEIIIS ]);

d = 1;
if nSEIIS>=1
    listParam{d} = paramSEIIS(1,:);
    d = d+1;
    if nSEIIS==2
        listParam{d} = paramSEIIS(2,:);
        d = d+1;
    end
end
if nSICR==1
    listParam{d} = paramSICR;
    d = d+1;
end
if nSEIIIS==1
    listParam{d} = paramSEIIIS;
    %d=d+1;
end
% 
% if nSEIIS==2 && nSICR==1 && nSEIIIS==1
%     listParam{2} = paramSEIIS(2,:);
%     listParam{3} = paramSICR;
%     listParam{4} = paramSEIIIS;
% end
% if nSEIIS==1 && nSICR==1 && nSEIIIS==0
%     listParam{2} = paramSICR;
% end

options = optimset('Display','off'); %options for minsearch

%vecBeta=[];vecGamma=[];
vecC0n = []; %stock c0 values (when rhohatn=0)
vecCnn = []; %stock cn values (when rhohatn=alphan)

lC = length(vecC);

%N3 = nchoosek(N, 3); % nb of comb. of 3 disease models  (e.g. SIS^2 x SIR)
%vecRhohat3d = zeros(lC,N3);
vecAlpha=[]; toc0 = toc;
%rhohat of 1-disease models
vecRhohat1d = zeros(lC,N); %initialization of the matrix that contains rhohat's values
for n=1:N
    type = modelTypes(n);
    disp(type);

    [vecRhohat1d_n,c0n,cnn,alpha] = findRhohat1d_v2(inputArg1,inputArg2);
    
    vecRhohat1d(:,n) = vecRhohat1d_n;
    vecC0n = [vecC0n,c0n];
    vecCnn = [vecCnn,cnn];
    vecAlpha =[vecAlpha,alpha];
    clear c0n cnn rhohatn alpha i 
    toc
end
tabRhohat.one = vecRhohat1d;
tabC0.one  = vecC0n;
tabCnn.one = vecCnn;
tabRhohat.single = 1:N;
toc1 = toc;
tabTimes.t1 = toc1 - toc0;

if N>=2 %if there are at least 2 diseases
    %rhohat of 2 disease models
    N2 = nchoosek(N, 2);        % nb of comb. of 2 disease models (e.g. SIS x SIR)
    duos = nchoosek(1:N,2);     % duets of diseases (sorted)
    vecRhohat2d = zeros(lC,N2); % initialization of the matrix that contains rhohat's values for the 2 disease model
    vecRhohatnm = zeros(lC,N2); % initialization of the matrix that contains rhohatij's values for the 2 disease model
    %possibilities : SIS2, SIRxSIR, SIR2, SISxSICAT, SICAT2, SIRxSICAT
    for nm = 1:N2
        n = duos(nm,1); m = duos(nm,2); %assigned numeros on diseases    
        disp([modelTypes(n),'x',modelTypes(m)]);
%         if (n==1 && m==3);             
%             disp([modelTypes(n),'x',modelTypes(m)]);
%         end
        alphan = vecAlpha(n); alpham = vecAlpha(m);
        maxAlpha = max(alphan,alpham);
        %omegan = paramSICAT(
        %omegan=0;omegam=0;sigman=0;sigmam=0;
        i=0;  
        for c=vecC
            c
            i=i+1;
            [rhohat,rhohatnm] = findRhohat2d_v2(c,modelTypes(n),modelTypes(m),vecRhohat1d(i,n),vecRhohat1d(i,m),listParam{n},listParam{m},mu,b,alphan,alpham,f);
            vecRhohat2d(i,nm) = rhohat;  
            vecRhohatnm(i,nm) = rhohatnm; %argmaxUnm
         
            if c==middlevecC
                disp('one half of c values left')
            end
        end
       %to do : to find c0 and cnn
       toc
    end
tabRhohat.two = vecRhohat2d;
tabRhohat.nm  = vecRhohatnm;
tabRhohat.duos = duos;
end

toc2 = toc;
tabTimes.t2 = toc2 - toc1;

if N>=3
    N3      = nchoosek(N, 3);  % nb of comb. of 3 disease models (e.g. SIS^2 x SIR)
    trios   = nchoosek(1:N,3); % trios of diseases (sorted)
    vecRhohat3d  = zeros(lC,N3); % initialization of the matrix that contains rhohat's values for the 3 disease model
    vecRhohatnmk = zeros(lC,N3); % initialization of the matrix that contains rhohatij's values for the 2 disease model
    for nmk = 1:N3
        n = trios(nmk,1); m = trios(nmk,2);k=trios(nmk,3);
        %reading parameters
%         betan = vecBeta(n); betam = vecBeta(m);betak = vecBeta(k);
%         gamman = vecGamma(n); gammam = vecGamma(m);gammak = vecGamma(k);
%         sn = vecSe(n); sm = vecSe(m);sk = vecSe(k);
        alphan = vecAlpha(n); alpham = vecAlpha(m);alphak = vecAlpha(k);
        disp([modelTypes(n),'x',modelTypes(m),'x',modelTypes(k)]);

        i=0;
        for c=vecC
            i=i+1;
            [rhohat,rhohatnmk] = findRhohat3d_v2(c,modelTypes(n),modelTypes(m),modelTypes(k),...
                vecRhohat1d(i,:),vecRhohat2d(i,:),listParam{n},listParam{m},listParam{k},mu,b,...
                alphan,alpham,alphak,f); %vecRhohat2d(i,:) for the value of c
            vecRhohat3d(i,nmk) = rhohat;  
            vecRhohatnmk(i,nmk) = rhohatnmk; %argmaxUnmk
            %c
             if c==middlevecC
                disp('one half of c values left')
            end
        end
        tabRhohat.three = vecRhohat3d;
        tabRhohat.nmk  = vecRhohatnmk;
        tabRhohat.trios = trios;
        toc
    end
    toc3 = toc;
    tabTimes.t3 = toc3 - toc2;
end

if N==4
    vecRhohatnmkl = zeros(lC,1); % initialization of the matrix that contains rhohatij's values for the 2 disease model
    vecRhohat4d = zeros(lC,1);
    disp('4 diseases')
    i=0;
    n=1;m=2;k=3;l=4; nmkl=1;
    alpha1 = vecAlpha(n);alpha2 = vecAlpha(m);alpha3=vecAlpha(k);alpha4=vecAlpha(l);
%         toc2=toc;
    for c=vecC
        c
        i=i+1;
        [rhohat,rhohatnmkl] = findRhohat4d_v2(c,modelTypes(1),modelTypes(2),modelTypes(3),modelTypes(4),...
            vecRhohat1d(i,:),vecRhohat2d(i,:),vecRhohat3d(i,:),listParam{n},listParam{m},listParam{k},listParam{l},mu,b,...
            alpha1,alpha2,alpha3,alpha4,f);
        vecRhohat4d(i,nmkl) = rhohat;  
        vecRhohatnmkl(i,nmkl) = rhohatnmkl; %argmaxUnmk

        tabRhohat.four = vecRhohat4d;
        tabRhohat.nmkl  = vecRhohatnmkl;
        %tabRhohat.quad = trios;
%             toc1 = toc;
%             toc1 - toc2
%             toc2 = toc1;
        if c==middlevecC
            disp('one half of c values left')
        end
    end
    toc
end

toc4 = toc;
tabTimes.t4 = toc4 - toc3;


if N==1
    tabRhohat.rhohat = vecRhohat1d;
end
if N==2
    tabRhohat.rhohat = vecRhohat2d;
end
if N==3
    tabRhohat.rhohat = vecRhohat3d;
end
if N==4
    tabRhohat.rhohat = vecRhohat4d;
end

end

