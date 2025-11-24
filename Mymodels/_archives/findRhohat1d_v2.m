function [outputArg1,outputArg2] = findRhohat1d_v2(inputArg1,inputArg2)
%n'a  aps ete terminé 

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    vecRhohat1d = zeros(lC);
    if type=="SEIIS" %beta1,gamma1,nu1,sigma1,eps1
        betan = paramSEIIS(n,1); gamman = paramSEIIS(n,2);
        nun = paramSEIIS(n,3); sigman = paramSEIIS(n,4); 
        epsn = paramSEIIS(n,5); 
        alpha = paramSEIIS(n,6);
        [~,c0n] = U_SEIISv2(epsn, betan, sigman, gamman, mu, nun, 0, 0,alpha, f); %c2 = c0
        [~,cnn] = U_SEIISv2(epsn, betan, sigman, gamman, mu, nun, alpha, 0,alpha, f); %c1=c1
        i=1;
        for c=vecC    
            fun = @(rho) -U_SEIISv2(epsn, betan, sigman, gamman, mu, nun, rho, c,alpha ,f);
            vecRhohat1d(i) = min(max(fminsearch(fun,0,options),0),alpha);
            i=i+1; %c
        end
    elseif type=="SICR"
        m = n-nSEIIS; %mth SIR disease %[betaISICR,betaCSICR,gammaSICR,sigmaSICR,thetaSICR,alphaSICR,RSICR];
        betaIm = paramSICR(m,1);betaCm  = paramSICR(m,2);
        gammam = paramSICR(m,3);sigmam  = paramSICR(m,4);
        thetam = paramSICR(m,5);
        alpha = paramSICR(m,6);
        [~,c0n] = U_SICR(betaIm, betaCm, sigmam, gammam, thetam, mu, 0, 0, alpha, f); %c2 = c0
        [~,cnn] = U_SICR(betaIm, betaCm, sigmam, gammam, thetam, mu, alpha, 0, alpha, f); %c1=c1
        i=1;
        for c=vecC 
            fun = @(rho) -U_SICR(betaIm, betaCm, sigmam, gammam, thetam, mu, rho, c,alpha, f);
            vecRhohat1d(i) = min(max(fminsearch(fun,0,options),0),alpha);
            i=i+1; %c
        end
    elseif type=="SEIIIS" %[betaS,sigmaS,tauS,gamma1S,thetaS,gamma3S,nuS,alphaS,SICR];
        k = n - nSEIIS - nSICR;
        betak = paramSEIIIS(k,1);sigmak = paramSEIIIS(k,2);
        tauk = paramSEIIIS(k,3);gamma1k = paramSEIIIS(k,4);
        thetak = paramSEIIIS(k,5);gamma3k = paramSEIIIS(k,6);
        nuk = paramSEIIIS(k,7);
        alpha = paramSEIIIS(k,8);
        [~,c0n] = U_SEIIIS(b, betak, sigmak, gamma1k, gamma3k, tauk, thetak, nuk, mu, 0, 0, alpha) ;
        [~,cnn] = U_SEIIIS(b, betak, sigmak, gamma1k, gamma3k, tauk, thetak, nuk, mu, alpha, 0,alpha) ;
        i=1;
        for c=vecC
            fun = @(rho) -U_SEIIIS(b, betak, sigmak, gamma1k, gamma3k, tauk, thetak, nuk, mu, rho, c,alpha);
            vecRhohat1d(i) = min(max(fminsearch(fun,0,options),0),alpha);
            i=i+1; %c
        end
    end
end

