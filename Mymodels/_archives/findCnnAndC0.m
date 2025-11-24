function [cMax0,cMin1,t0,t1] = findCnnAndC0(cMin,cMax,Tmax,errMax,N,paramTab,mu,b,f)
    % GENERALISATION OF findCnnAndC0_2d
    % -------------
    % of the N diseases, and c0 is the cost where rhohat=0
    %%
    % PROBLEME : dans ce code on essaie d'encadrer les valeurs de c 
    % pour la modele à 4 maladies, mais pas pour les sous modèles
    %% 
    % find the two thresholds cnn and c0 where cnn is the cost for elimination
    %%
    vecAlpha=zeros(1,N);
    for i=1:N
        vecAlpha(i) = paramTab{i}.alpha;
    end
    alphaMax = max(vecAlpha);
    %------------------
    %Finding c0 
    disp('recherche de c0')
    t0=0;
    step0 = errMax*2;
    cMin0=0; cMax0=cMax; %on pose cMin=0, car en general c0>0
    
    %on verifie d'abord que les bornes encadrent bien c00 :
    %sinon, on refinie les bornes : 
    %en cMax : on doit rhohat=0 et en cMin : rhohat>0
    %disp('pass1')
    borneSupOK = 0;k=1;
    while(~borneSupOK)
        [tabRhohat,~,~,~] = findRhohat_v4(N,paramTab,mu,b,cMax0,f);
        rhohat = tabRhohat.rhohat.rhohat;
        if (rhohat==0) %if rhohat=0 then c0<cMax0
            borneSupOK = 1;
        else
            cMax0 = abs(cMax0+errMax); %*2? %if rhohat>0 then c0>cMax0, we looking for bigger values of cMax0
        end
        k=k+1;
    end
    borneInfOK = 0; k=1;
    %disp('pass2')
    while(~borneInfOK)
        [tabRhohat,~,~,~] = findRhohat_v4(N,paramTab,mu,b,cMin0,f);
        rhohat = tabRhohat.rhohat.rhohat;
        k=k+1;
        if (rhohat>0)
            borneInfOK = 1;
        else %rhohat==0 means that cMin0>c0, we need to look at values of cMin0 such that cMin0<c0
            cMin0 = -2*abs(cMin0+errMax);
        end
    end
    %disp('pass3') 
    k=0;
    %les bornes encadrent bien la valeur de c00, donc on procede par
    %dichotomie
    while t0<Tmax && abs(step0)>errMax
        k=k+1;
        %disp(['cmin=', num2str(cMin0), ' cmax=', num2str(cMax0)])
        c0 = (cMin0 + cMax0)/2;
        [tabRhohat,~,~,~] = findRhohat_v4(N,paramTab,mu,b,c0,f);    
        %[rhohat,~] = findRhohat2d_v2(c0,type1,type2,rhohat1,rhohat2,param1,param2,mu,b,alpha1,alpha2,f);   
        rhohat = tabRhohat.rhohat.rhohat;
        if rhohat>0
            cMin0 = c0;
        elseif rhohat<=0
            cMax0 = c0;
        end  
       t0=t0+1;
       step0 = abs(cMax0-cMin0);
    end
    %c0 = (cMax0+cMin0)/2;
    
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    %Finding cnn
    %cnn = (cMax1+cMin1)/2;
    disp('recherche de cnn')
  
    %on calcule cnn de la maladie eliminee en dernier
    [~,imax] = max(vecAlpha);
    type = paramTab{imax}.modelType;

    if isequal(type,'SICR')
        %[~,c0i] = U1_SICR_v4(paramTab{imax},mu,0,0,f); %c2 = c0
        [~,cii] = U1_SICR_v4(paramTab{imax},mu,b,paramTab{imax}.alpha,0,f); %c1=c1
    elseif isequal(type,'SEIIIS')
        %[~,c0i] = U1_SEIIIS_v4(paramTab{imax},mu,0,0,f);
        [~,cii] = U1_SEIIIS_v4(paramTab{imax},mu,b,paramTab{imax}.alpha,0,f);
    elseif isequal(type,'SEIIS')
        %[~,c0i] = U_SEIIS_v4(paramTab{imax},mu,0,0,f); %c2 = c0
        [~,cii] = U1_SEIISv4(paramTab{imax},mu,b,paramTab{imax}.alpha,0,f); %c1=c1
    end
    %disp('pass4')
    %on affine autour de cii
    t1=0;
    step1 = errMax*2;
    cMin1 = -1.2*abs(cii); cMax1 = 0;
    %on verifie d'abord que les bornes encadrent bien cnn :
    %sinon, on refinie les bornes : 
    %en cMax, on doit avoir : rhohat<alphaMax et en cMin : rhohat=alphaMax
    borneInfOK = 0; k=0;
    while(~borneInfOK)
        [tabRhohat,~,~,~] = findRhohat_v4(N,paramTab,mu,b,cMin1,f);
        k=k+1;
        rhohat = tabRhohat.rhohat.rhohat;
        if (rhohat>=alphaMax)
            borneInfOK = 1;
        else
            cMin1 = -2*abs(cMin1-errMax); %le terme errMax vient palier le pb si cMin=0
        end
    end
    %disp('pass5')
    k=0;
    borneSupOK = 0;
    while(~borneSupOK)
        k=k+1;
        [tabRhohat,~,~,~] = findRhohat_v4(N,paramTab,mu,b,cMax1,f);
        rhohat = tabRhohat.rhohat.rhohat;
        if (rhohat<alphaMax)
            borneSupOK = 1;
        else
            cMax1 = 2*abs(cMax1+errMax);
        end
    end
    %disp('pass6')
    while t1<Tmax && abs(step1)>errMax
        %disp(['cmin=', num2str(cMin1), ' cmax=', num2str(cMax1)])
        cnn = (cMin1 + cMax1)/2;
        [tabRhohat,~,~,~] = findRhohat_v4(N,paramTab,mu,b,cnn,f);  
        rhohat = tabRhohat.rhohat.rhohat;
        if rhohat>=alphaMax
            cMin1 = cnn;
        elseif rhohat<alphaMax
            cMax1 = cnn;
        end  
        t1=t1+1;
        step1 = abs(cMax1-cMin1);
    end
    cMin1;
    cMax0;
    %disp(['cnn<-(',num2str(cnn),'), c0<-(',num2str(c0),')'])
    %disp(['t0=',num2str(t0), 't1=', num2str(t1)])
end

