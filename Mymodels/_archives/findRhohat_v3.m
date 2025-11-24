function [tabRhohat,tabC0,tabCnn,tabTimes] = findRhohat_v3(N,paramTab,mu,b,vecC,f)
%for any combination of SIS^nSIS x SIR^nSIR x SICAT^nSICAT model, this
%function finds the argmax of single diseases models, two-disease
%models,...,N diseases models
%difference with findRhohat_v2 : parameters are stored in structures
tic; middlevecC = median(vecC);
vecC0n = []; %stock c0 values (when rhohatn=0)
vecCnn = []; %stock cn values (when rhohatn=alphan)

lC = length(vecC);
vecAlpha=[]; toc0 = toc;
vecRhohat1d = []; %initialization of the matrix that contains rhohat's values
for n=1:N
    type = paramTab{n}.modelType;
    %disp(type);
    [vecRhohat1d_n,c0n,cnn] = findRhohat1d_v3(type,paramTab{n},mu,vecC,f);   
    vecRhohat1d(n).rhohat = vecRhohat1d_n;
    vecRhohat1d(n).modelType = type;
    vecRhohat1d(n).disease = paramTab{n}.disease;
    
    vecC0n = [vecC0n,c0n];
    vecCnn = [vecCnn,cnn];
    %vecAlpha =[vecAlpha,alpha];
    %clear c0n cnn rhohatn alpha i 
    %toc;
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
    %vecRhohat2d = zeros(1,N2); % initialization of the matrix that contains rhohat's values for the 2 disease model
    %vecRhohatnm = zeros(1,N2); % initialization of the matrix that contains rhohatij's values for the 2 disease model
    %possibilities : SIS2, SIRxSIR, SIR2, SISxSICAT, SICAT2, SIRxSICAT
    for nm = 1:N2
        n = duos(nm,1); m = duos(nm,2); %assigned numeros on diseases 
        typen = paramTab{n}.modelType;
        typem = paramTab{m}.modelType;
        %disp([typen,'x',typem]);
        alphan = paramTab{n}.alpha; alpham = paramTab{m}.alpha;
        i=0;
        
        vecRhohat2d(nm).rhohat = zeros(1,lC);
        vecRhohatnm(nm).rhohat = zeros(1,lC);
        vecRhohat2d(nm).type1 = typen;
        vecRhohat2d(nm).type2 = typem;
        vecRhohatnm(nm).type1 = typen;
        vecRhohatnm(nm).type2 = typem;
        
        for c=vecC
            %disp(c);
            i=i+1;
            [rhohat,rhohatnm] = findRhohat2d_v2(c,typen,typem,vecRhohat1d(n).rhohat(i),...
                vecRhohat1d(m).rhohat(i),paramTab{n},paramTab{m},mu,b,alphan,alpham,f);
            vecRhohat2d(nm).rhohat(i) = rhohat;
            vecRhohatnm(nm).rhohat(i) = rhohatnm;
            %if c==middlevecC
            %    disp('one half of c values left')
            %end
        end
       %to do : to find c0 and cnn
       %toc;
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
    %vecRhohat3d  = zeros(1,N3); % initialization of the matrix that contains rhohat's values for the 3 disease model
    %vecRhohatnmk = zeros(1,N3); % initialization of the matrix that contains rhohatij's values for the 2 disease model
    for nmk = 1:N3
        n = trios(nmk,1); m = trios(nmk,2);k=trios(nmk,3);
        typen = paramTab{n}.modelType;
        typem = paramTab{m}.modelType;
        typek = paramTab{k}.modelType;
        alphan = paramTab{n}.alpha; alpham = paramTab{m}.alpha; alphak = paramTab{k}.alpha;
        %disp([typen,'x',typem,'x',typek]);
        
        vecRhohat3d(nmk).rhohat = zeros(1,lC);
        vecRhohatnmk(nmk).rhohat = zeros(1,lC);
        vecRhohat3d(nmk).type1 = typen;
        vecRhohat3d(nmk).type2 = typem;
        vecRhohat3d(nmk).type3 = typek;
        vecRhohatnmk(nmk).type1 = typen;
        vecRhohatnmk(nmk).type2 = typem;
        vecRhohatnmk(nmk).type3 = typek;
        
        rhohat1d_reshape = [vecRhohat1d([n,m,k]).rhohat]; %on recupere les rhohat a 1 maladie
        vecRhohat2d_nmk = vecRhohat2d((tabRhohat.duos(:,1)==n|tabRhohat.duos(:,1)==m|tabRhohat.duos(:,1)==k)&...
            (tabRhohat.duos(:,2)==n|tabRhohat.duos(:,2)==m|tabRhohat.duos(:,2)==k));
        rhohat2d_reshape = reshape([vecRhohat2d_nmk(:).rhohat],[],3);
        i=0;
        for c=vecC
            i=i+1;
            [rhohat,rhohatnmk] = findRhohat3d_v3(c,typen,typem,typek,...
                rhohat1d_reshape(i,:),rhohat2d_reshape(i,:),paramTab{n},paramTab{m},paramTab{k},mu,b,...
                alphan,alpham,alphak,f); 
            vecRhohat3d(nmk).rhohat(i) = rhohat;
            vecRhohatnmk(nmk).rhohat(i) = rhohatnmk; %argmaxUnmk
            %c
             %if c==middlevecC
             %   disp('one half of c values left')
             %end
        end
        tabRhohat.three = vecRhohat3d;
        tabRhohat.nmk  = vecRhohatnmk;
        tabRhohat.trios = trios;
        %toc;
    end
end

toc3 = toc;
tabTimes.t3 = toc3 - toc2;

if N==4
    vecRhohatnmkl = zeros(lC,1); % initialization of the matrix that contains rhohatij's values for the 2 disease model
    vecRhohat4d = zeros(lC,1);
    %disp('4 diseases')
    i=0;
    rhohat3d_reshape = reshape([vecRhohat3d(:).rhohat],[],N3);
    for c=vecC
        %c
        i=i+1;
        [rhohat,rhohatnmkl] = findRhohat4d_v2(c,paramTab{1}.modelType,paramTab{2}.modelType,paramTab{3}.modelType,paramTab{4}.modelType,...
            rhohat1d_reshape(i),rhohat2d_reshape(i),rhohat3d_reshape(i),paramTab{1},paramTab{2},paramTab{3},paramTab{4},mu,b,...
            paramTab{1}.alpha,paramTab{2}.alpha,paramTab{3}.alpha,paramTab{4}.alpha,f);
      
        vecRhohat4d(i,1) = rhohat;  
        vecRhohatnmkl(i,1) = rhohatnmkl; %argmaxUnmk

        tabRhohat.four = vecRhohat4d;
        tabRhohat.nmkl  = vecRhohatnmkl;
        %tabRhohat.quad = trios;
%             toc1 = toc;
%             toc1 - toc2
%             toc2 = toc1;
        %if c==middlevecC
        %    disp('one half of c values left')
        %end
    end
    %toc;
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

