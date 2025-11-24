function [c0,cnn,t0,t1] = findCnnAndC0_2d(cMin,cMax,Tmax,errMax,param1,param2,mu,b,f)
    % USE findCnnAndC0(cMin,cMax,Tmax,errMax,param1,param2,mu,b,f) INSTEAD
    %find the two thresholds cnn and c0 where cnn is the cost for elimination
    %of the 2 diseases, and c0 is the cost where rhohat=0
    paramTab{1}=param1;paramTab{2}=param2;
    N=2;
    alpha1 = param1.alpha;alpha2=param2.alpha;
    alphaMax = max(alpha1,alpha2);
    %------------------
    %Finding c0 
    t0=0;
    step0 = errMax*2;
    cMin0=0; cMax0=cMax; %on pose cMin=0, car en general c0>0
    
    %on verifie d'abord que les bornes encadrent bien c00 :
    %sinon, on refinie les bornes : 
    %en cMax : on doit rhohat=0 et en cMin : rhohat>0
    borneSupOK = 0;
    while(~borneSupOK)
        [tabRhohat,~,~,~] = findRhohat_v3(N,paramTab,mu,b,cMax0,f);
        rhohat = tabRhohat.two;
        if (rhohat==0)
            borneSupOK = 1;
        else
            cMax0 = 2*abs(cMax+errMax);
        end
    end
    borneInfOK = 0;
    while(~borneInfOK)
        [tabRhohat,~,~,~] = findRhohat_v3(N,paramTab,mu,b,cMin0,f);
        rhohat = tabRhohat.two;
        if (rhohat>0)
            borneInfOK = 1;
        else
            cMin0 = -2*abs(cMin+errMax);
        end
    end
    
    %les bornes encadrent bien la valeur de c00, donc on procede par
    %dichotomie
    while t0<Tmax && abs(step0)>errMax
        c0 = (cMin0 + cMax0)/2;
        [tabRhohat,~,~,~] = findRhohat_v3(N,paramTab,mu,b,c0,f);    
        %[rhohat,~] = findRhohat2d_v2(c0,type1,type2,rhohat1,rhohat2,param1,param2,mu,b,alpha1,alpha2,f);   
        rhohat = tabRhohat.two;
        if rhohat>0
            cMin0 = c0;
        elseif rhohat<=0
            cMax0 = c0;
        end  
       t0=t0+1;
       step0 = abs(cMax0-cMin0);
    end
    c0 = (cMax0+cMin0)/2;
    
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    %Finding cnn
    t1=0;
    step1 = errMax*2;
    cMin1 = cMin; cMax1 = 0;
    %on verifie d'abord que les bornes encadrent bien cnn :
    %sinon, on refinie les bornes : 
    %en cMax : on doit rhohat<alphaMax et en cMin : rhohat=alphaMax
    borneInfOK = 0;
    while(~borneInfOK)
        [tabRhohat,~,~,~] = findRhohat_v3(N,paramTab,mu,b,cMin1,f);
        rhohat = tabRhohat.two;
        if (rhohat>=alphaMax)
            borneInfOK = 1;
        else
            cMin1 = -2*abs(cMin1+errMax); %le terme errMax vient palier le pb si cMin=0
        end
    end
    
    borneSupOK = 0;
    while(~borneSupOK)
        [tabRhohat,~,~,~] = findRhohat_v3(N,paramTab,mu,b,cMax1,f);
        rhohat = tabRhohat.two;
        if (rhohat<alphaMax)
            borneSupOK = 1;
        else
            cMax1 = 2*abs(cMax+errMax);
        end
    end
    
    while t1<Tmax && abs(step1)>errMax
        cnn = (cMin1 + cMax1)/2;
        [tabRhohat,~,~,~] = findRhohat_v3(N,paramTab,mu,b,cnn,f);  
        rhohat = tabRhohat.two;
        if rhohat>=alphaMax
            cMin1 = cnn;
        elseif rhohat<alphaMax
            cMax1 = cnn;
        end  
       t1=t1+1;
       step1 = abs(cMax1-cMin1);
    end
    cnn = (cMax1+cMin1)/2;
end

