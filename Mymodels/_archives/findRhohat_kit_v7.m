function [vec_rho_hat,ES,status] = findRhohat_kit_v7(paramTab_all,mu,b,f,paramRho,kit,vecC,alphas,verbose,mySeed,log_path)
%for any combination of SICTP x SEIIIS x SEIIS^m, this
%function finds the argmax of single diseases models, two-disease
%models,...,N diseases models

GlobalOpt.tolP0=1e-5;
GlobalOpt.verbose=0;
GlobalOpt.TolFun=1e-8;
GlobalOpt.inf_bnd_alpha=0;

includeHIV = ~alphas.elim_h;

Nk = length(kit); %nombre d'infections dans le kit
for i=1:4
    for j=1:Nk
        if strcmp(paramTab_all{i}.disease,kit{j})
            paramTab{j} = paramTab_all{i};
        end
    end
end

vecC = sort(vecC);
vecRhohat1d = []; %initialization of the matrix that contains rhohat's values

vecC0n = []; %store c0 values (when rhohatn=0)
vecCnn = []; %store cn values (when rhohatn=alphan)
lC = length(vecC);

%% Calcul de rhohat
%Single-disease sub-model
tic;
%On calcule rhohat pour les modèles à 1 infection,
%i.e. quand il ne reste plus que celle là
opt = GlobalOpt;
for n=1:Nk
    type = paramTab{n}.modelType;
    k = paramTab{n}.mini_d;
    mod = k;
    if includeHIV && ~contains(k,'h')
        mod = ['h',k];
    end
    i=1;
    opt.up_bnd_alpha  = alphas.([k,'_',mod]);
    opt.inf_bnd_alpha = 0;
    for c=vecC
        [vecRhohat1d_n,c0n,cnn] = fminU_knitro_v7_bis(kit,mod,paramTab,paramRho,mu,b,c,opt);
        %[vecRhohat1d_n,c0n,cnn] = findRhohat1d_v7(type,paramTab{n},mu,b,vecC,f);
        vecRhohat1d(n).rhohat(i) = vecRhohat1d_n;
        i=i+1;
    end
    vecRhohat1d(n).modelType = type;
    vecRhohat1d(n).disease = paramTab{n}.disease;
    vecC0n = [vecC0n,c0n];
    vecCnn = [vecCnn,cnn];
end
tabRhohat.one = vecRhohat1d;
tabC0.one  = vecC0n;
tabCnn.one = vecCnn;
tabRhohat.single = 1:Nk;
toc1 = toc;
tabTimes.t1 = toc1;

%% Two-disease submodel
if  Nk>=2 
    %rhohat of 2 disease models
    N2 = nchoosek(Nk, 2);        % nb of comb. of 2 disease models (e.g. SIS x SIR)
    duos = nchoosek(1:Nk,2);     % duets of diseases (sorted)
    %vecRhohat2d = zeros(1,N2); % initialization of the matrix that contains rhohat's values for the 2 disease model
    %vecRhohatnm = zeros(1,N2); % initialization of the matrix that contains rhohatij's values for the 2 disease model
    %possibilities : SICTPxSEIIIS, SICTPxSEIIS(Ct), SICTPxSEIIS(Ng),
    %SEIIISxSEIIS(Ct), SEIIISxSEIIS(Ng), SEIIS2
    for nm = 1:N2
        n = duos(nm,1); m = duos(nm,2); %assigned numeros on diseases 
        typen = paramTab{n}.modelType;
        typem = paramTab{m}.modelType;
        disp([typen,'x',typem]);
        %alphan = paramTab{n}.alpha; alpham = paramTab{m}.alpha;
        vecRhohat2d(nm).rhohat = zeros(1,lC);
        vecRhohatnm(nm).rhohat = zeros(1,lC);
        vecRhohat2d(nm).type1 = typen;
        vecRhohat2d(nm).type2 = typem;
        vecRhohatnm(nm).type1 = typen;
        vecRhohatnm(nm).type2 = typem;
        vecRhohat2d(nm).inf1 = paramTab{n}.disease;
        vecRhohat2d(nm).inf2 = paramTab{m}.disease;
         
        %[alphaSorted,iSorted] = sort([paramTab{n}.alpha, paramTab{m}.alpha]);
        %alphaMin = alphaSorted(1);
        numDis = [n,m];
        vecRhohat = zeros(1,length(vecC));
        i=0;
        for c=flip(vecC)
            %disp(c); 
            %[rhohat,rhohatnm] = findRhohat2d_v4(c,typen,typem,vecRhohat1d(n).rhohat(lC-i),...
            %    vecRhohat1d(m).rhohat(lC-i),paramTab{n},paramTab{m},mu,b,f);
            opt.up_bnd_alpha = alphas.([paramTab{n}.mini_d,paramTab{m}.mini_d,'_',paramTab{n}.mini_d,paramTab{m}.mini_d]);
            [rhohat] = findRhohat2d_v7(c,typen,typem,paramTab{n},paramTab{m},paramRho,mu,b,f,kit,opt);
            vecRhohat2d(nm).rhohat(lC-i) = rhohat;
            %vecRhohatnm(nm).rhohat(lC-i) = rhohatnm;
            %if c==middlevecC
            %    disp('one half of c values left')
            %end
            i=i+1;

        end
        % tout le reste du vecteur, on le met à NA (ce sont des valeurs
        % supposees plus grandes de rhohat, mais comme rhohat>=min(alpha),
        % alors on passe directement au modèle à 1 maladie.
        vecRhohat2d(nm).rhohat(1:(lC-i)) = vecRhohat1d(numDis(iSorted(2))).rhohat(1:(lC-i));%rhohat m ou n; 
        vecRhohatnm(nm).rhohat(1:(lC-i)) = NaN(1,lC-i);%rhohatnm;
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
        
        %alphan = paramTab{n}.alpha; alpham = paramTab{m}.alpha; alphak = paramTab{k}.alpha;
        %disp([typen,'x',typem,'x',typek]);
        
        vecRhohat3d(nmk).rhohat = zeros(1,lC);
        vecRhohatnmk(nmk).rhohat = zeros(1,lC);
        vecRhohat3d(nmk).type1 = typen;
        vecRhohat3d(nmk).type2 = typem;
        vecRhohat3d(nmk).type3 = typek;
        vecRhohatnmk(nmk).type1 = typen;
        vecRhohatnmk(nmk).type2 = typem;
        vecRhohatnmk(nmk).type3 = typek;
        vecRhohat3d(nmk).inf1 = paramTab{n}.disease;
        vecRhohat3d(nmk).inf2 = paramTab{m}.disease;
        vecRhohat3d(nmk).inf3 = paramTab{k}.disease;
        
        rhohat1d_reshape = [vecRhohat1d([n,m,k]).rhohat]; %on recupere les rhohat a 1 maladie
        vecRhohat2d_nmk = vecRhohat2d((tabRhohat.duos(:,1)==n|tabRhohat.duos(:,1)==m|tabRhohat.duos(:,1)==k)&...
            (tabRhohat.duos(:,2)==n|tabRhohat.duos(:,2)==m|tabRhohat.duos(:,2)==k));
        rhohat2d_reshape = reshape([vecRhohat2d_nmk(:).rhohat],[],3);
        
        [alphaSorted,iSorted] = sort([paramTab{n}.alpha, paramTab{m}.alpha,paramTab{k}.alpha]);
        alphaMin = alphaSorted(1);
        numDis = [n,m,k];
        
        rhohat = 0; i=0;
        %for c=vecC
        while (rhohat<alphaMin & i<lC)
            c = vecC(lC-i);
            [rhohat,rhohatnmk] = findRhohat3d_v4(c,typen,typem,typek,...
                rhohat1d_reshape(lC-i,:),rhohat2d_reshape(lC-i,:),...
                paramTab{n},paramTab{m},paramTab{k},mu,b,f); 
            vecRhohat3d(nmk).rhohat(lC-i) = rhohat;
            vecRhohatnmk(nmk).rhohat(lC-i) = rhohatnmk; %argmaxUnmk
            %c
             %if c==middlevecC
             %   disp('one half of c values left')
             %end
             i=i+1;
        end
        % tout le reste du vecteur, on le met à NA (ce sont des valeurs
        % supposees plus grandes de rhohat, mais comme rhohat>=min(alpha),
        % alors on passe directement au modèle à 2 maladies 
        % (avec les 2 maladies qui restent).
        %numDis(iSorted(1:2)) %2 diseases that remains
        rhohatRemain2d = vecRhohat2d((tabRhohat.duos(:,1)==numDis(iSorted(2))&tabRhohat.duos(:,2)==numDis(iSorted(3)))|...
            (tabRhohat.duos(:,1)==numDis(iSorted(3))&tabRhohat.duos(:,2)==numDis(iSorted(2)))).rhohat;
        vecRhohat3d(nmk).rhohat(1:(lC-i)) = rhohatRemain2d(1:(lC-i));
        vecRhohatnmk(nmk).rhohat(1:(lC-i)) = NaN(1,lC-i);%rhohatnm; we don't need this information
        
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
    %vecRhohat4d = zeros(lC,1);
    %disp('4 diseases')
    vecRhohat4d.type1='SEIIS';
    vecRhohat4d.type2='SEIIS';
    vecRhohat4d.type3='SICR';
    vecRhohat4d.type4='SEIIIS';
    
    rhohat3d_reshape = reshape([vecRhohat3d(:).rhohat],[],N3);
    
    [alphaSorted,iSorted] = sort([paramTab{1}.alpha, paramTab{2}.alpha,...
                                  paramTab{3}.alpha,paramTab{4}.alpha]);
    alphaMin = alphaSorted(1);
    numDis = [1,2,3,4];

    rhohat = 0; i=0;
    while (rhohat<alphaMin & i<lC)
    %for c=vecC      
        c = vecC(lC-i);
        
        [rhohat,rhohatnmkl] = findRhohat4d_v4(c,...
            paramTab{1}.modelType,paramTab{2}.modelType,...
            paramTab{3}.modelType,paramTab{4}.modelType,...
            rhohat1d_reshape(lC-i),rhohat2d_reshape(lC-i),rhohat3d_reshape(lC-i),...
            paramTab{1},paramTab{2},paramTab{3},paramTab{4},mu,b,f);
      
        vecRhohat4d.rhohat(lC-i) = rhohat;
        %vecRhohat4d(i,1) = rhohat;  
        vecRhohatnmkl(lC-i,1) = rhohatnmkl; %argmaxUnmk   
        %if c==middlevecC
        %    disp('one half of c values left')
        %end
        i=i+1;
    end
    NumDisRem = (1:4)~=iSorted(1);
    rhohatRemain3d = vecRhohat3d(mean(tabRhohat.trios(:,:)==numDis(NumDisRem),2)==1).rhohat;
    vecRhohat4d.rhohat(1:(lC-i)) = rhohatRemain3d(1:(lC-i));
    vecRhohatnmkl(1:(lC-i)) = NaN(1,lC-i);%rhohatnm; we don't need this information
    %toc;
    tabRhohat.four = vecRhohat4d;
    tabRhohat.nmkl  = vecRhohatnmkl;
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

tabRhohat.c=vecC;

end

