function [tabRhohat,tabC0,tabCnn] = findRhohat(nSIS,nSIR,nSICAT,paramSIS,paramSIR,paramSICAT,mu,b,vecC)
%for any combination of SIS^nSIS x SIR^nSIR x SICAT^nSICAT model, this
%function finds the argmax of single diseases models, two-disease
%models,...,N diseases models

%checks
if nSIS~=size(paramSIS,1)
   error('Error. The number of SIS diseases does not match with the parameters')
end
if nSIR~=size(paramSIR,1)
   error('Error. The number of SIR diseases does not match with the parameters')
end
if nSICAT~=size(paramSICAT,1)
   error('Error. The number of SICAT diseases does not match with the parameters')
end

N = nSIS + nSIR + nSICAT; %total number of diseases in the model
modelTypes = repelem(["SIS", "SIR", "SICAT"], [nSIS nSIR nSICAT]);

vecBeta=[];vecGamma=[];vecSe=[];
vecC0n = []; %stock c0 values (when rhohatn=0)
vecCnn = []; %stock cn values (when rhohatn=alphan)
if ~isempty(paramSIS)
    vecBeta  =  [vecBeta ,paramSIS(:,1)'];
    vecGamma =  [vecGamma,paramSIS(:,2)'];
    vecSe    =  [vecSe,paramSIS(:,3)'];
end
if ~isempty(paramSIR)
    vecBeta  =  [vecBeta ,paramSIR(:,1)'];
    vecGamma =  [vecGamma,paramSIR(:,2)'];
    vecSe    =  [vecSe,paramSIR(:,3)'];
end
vecOmega=zeros(nSIS+nSIR,1);
vecSigma=zeros(nSIS+nSIR,1);
if ~isempty(paramSICAT)
    vecBeta  =  [vecBeta ,paramSICAT(:,1)'];
    vecGamma =  [vecGamma,paramSICAT(:,2)'];
    vecSe    =  [vecSe,paramSICAT(:,3)'];
    vecOmega =  [vecOmega,paramSICAT(:,4)'];
    vecSigma =  [vecSigma,paramSICAT(:,5)'];
end

vecR0SISR    = vecBeta./(vecGamma+mu); %R0 SIS and SIR
%to do : R0 and alpha SICAT
%vecR0SICAT   = ;
vecAlpha     = vecBeta./vecSe.*(1-1./vecR0SISR); %rho' SIS and SIR
%to do : same for SICAT

%initialization of matrices that contain rhohat's value
lC = length(vecC);

%N3 = nchoosek(N, 3); % nb of comb. of 3 disease models  (e.g. SIS^2 x SIR)
%vecRhohat3d = zeros(lC,N3);

%rhohat of 1-disease models
vecRhohat1d = zeros(lC,N); %initialization of the matrix that contains rhohat's values
for n=1:N
    type = modelTypes(n);
    if type=="SIS"
        betan = paramSIS(n,1); sn = paramSIS(n,3);
        alphan = vecAlpha(n);
        R0n    = vecR0SISR(n);
        cnn = 1/R0n-1; c0n = -cnn;
        rhohatn = (alphan/2-vecC*betan/(2*sn)).*(vecC<c0n & vecC>cnn) + alphan.*(vecC<=cnn);
        vecRhohat1d(:,n) = rhohatn;
    elseif type=="SIR"
        m = n-nSIS; %mth SIR disease
        betam  = paramSIR(m,1); sm = paramSIR(m,3);
        alpham = vecAlpha(n);
        R0m    = vecR0SISR(m);
        cnn = -mu/betam*(1-1/R0m);c0n = mu/betam*(R0m-1);
        rhohatn = (betam./R0m/sm.*(sqrt(R0m*mu./(mu+betam*vecC))-1)).*(vecC<c0n & vecC>cnn) + alpham.*(vecC<=cnn);
        vecRhohat1d(:,n) = rhohatn;
    elseif type=="SICAT"
        % to be written
    end
    vecC0n = [vecC0n,c0n];
    vecCnn = [vecCnn,cnn];
    clear R0n cnn rhohatn betan gamman alphan
end
tabRhohat.one = vecRhohat1d;
tabC0.one  = vecC0n;
tabCnn.one = vecCnn;
tabRhohat.single = 1:N;

if N>=2 %if there are at least 2 diseases
    %rhohat of 2 disease models
    N2 = nchoosek(N, 2);        % nb of comb. of 2 disease models (e.g. SIS x SIR)
    duos = nchoosek(1:N,2);     % duets of diseases (sorted)
    vecRhohat2d = zeros(lC,N2); % initialization of the matrix that contains rhohat's values for the 2 disease model
    vecRhohatnm = zeros(lC,N2); % initialization of the matrix that contains rhohatij's values for the 2 disease model
    %possibilities : SIS2, SIRxSIR, SIR2, SISxSICAT, SICAT2, SIRxSICAT
    for nm = 1:N2
        n = duos(nm,1); m = duos(nm,2); %assigned numeros on diseases        
        %reading parameters
        betan = vecBeta(n); betam = vecBeta(m);
        gamman = vecGamma(n); gammam = vecGamma(m);
        sn = vecSe(n); sm = vecSe(m);
        alphan = vecAlpha(n); alpham = vecAlpha(m);
        %omegan = paramSICAT(
        omegan=0;omegam=0;sigman=0;sigmam=0;

        i=0;
        for c=vecC
            i=i+1;
            [rhohat,rhohatnm] = findRhohat2d(c,modelTypes(n),modelTypes(m),vecRhohat1d(i,n),vecRhohat1d(i,m),betan,betam,gamman,gammam,sn,sm,mu,b,omegan,omegam,sigman,sigmam,alphan,alpham);
            vecRhohat2d(i,nm) = rhohat;  
            vecRhohatnm(i,nm) = rhohatnm; %argmaxUnm
        end
       %to do : to find c0 and cnn
    end
tabRhohat.two = vecRhohat2d;
tabRhohat.nm  = vecRhohatnm;
tabRhohat.duos = duos;
end

if N>=3
    N3      = nchoosek(N, 3);  % nb of comb. of 3 disease models (e.g. SIS^2 x SIR)
    trios   = nchoosek(1:N,3); % trios of diseases (sorted)
    vecRhohat3d  = zeros(lC,N3); % initialization of the matrix that contains rhohat's values for the 3 disease model
    vecRhohatnmk = zeros(lC,N3); % initialization of the matrix that contains rhohatij's values for the 2 disease model
    for nmk = 1:N3
        n = trios(nmk,1); m = trios(nmk,2);k=trios(nmk,3);
        %reading parameters
        betan = vecBeta(n); betam = vecBeta(m);betak = vecBeta(k);
        gamman = vecGamma(n); gammam = vecGamma(m);gammak = vecGamma(k);
        sn = vecSe(n); sm = vecSe(m);sk = vecSe(k);
        alphan = vecAlpha(n); alpham = vecAlpha(m);alphak = vecAlpha(k);

        %a modifier avec l'arrivee de sicat
        omegan=0;omegam=0;sigman=0;sigmam=0;
        i=0;
        for c=vecC
            i=i+1;
            [rhohat,rhohatnmk] = findRhohat3d(c,modelTypes(n),modelTypes(m),modelTypes(k),...
                vecRhohat1d(i,:),vecRhohat2d(i,:),betan,betam,betak,gamman,gammam,gammak,sn,sm,sk,mu,b,...
                omegan,omegam,sigman,sigmam,alphan,alpham,alphak);
            vecRhohat3d(i,nmk) = rhohat;  
            vecRhohatnmk(i,nmk) = rhohatnmk; %argmaxUnmk
        end
        tabRhohat.three = vecRhohat3d;
        tabRhohat.nmk  = vecRhohatnmk;
        tabRhohat.trios = trios;
    end
end

if N==1
    tabRhohat.rhohat = vecRhohat1d;
end
if N==2
    tabRhohat.rhohat = vecRhohat2d;
end
if N==3
    tabRhohat.rhohat = vecRhohat3d;
end
 %a modifier
%to do : for 3 diseases and more

end

