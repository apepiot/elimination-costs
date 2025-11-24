function [tabcs] = findThresholds_v5(nSEIIS,nSICR,nSEIIIS,tab,vecAlpha,vecC)
    N = nSEIIS+nSICR+nSEIIIS; 
    if N==1
        param = tab.param;
        if nSEIIS==1
            [~,tabcs.c0] = U1_SEIISv4(param,param.mu,param.b,0,0,1);
            [~,tabcs.cs1dis] = U1_SEIISv4(param,param.mu,param.b,param.alpha,0,1); 
        elseif nSEIIIS==1
            [~,tabcs.c0] = U1_SEIIIS_v4(param,param.mu,param.b,0,0,1);
            [~,tabcs.cs1dis] = U1_SEIIIS_v4(param,param.mu,param.b,param.alpha,0,1); 
        elseif nSICR==1
            [~,tabcs.c0] = U1_SICTP(param,param.mu,param.b,0,0,1);
            [~,tabcs.cs1dis] = U1_SICTP(param,param.mu,param.b,param.alpha,0,1); 
        end
                
    end
    
    if(N>1)
        if(N==2)
            %U=Uij if rho<rhoj', U=Ui if rho>=rhoj'
            [vecAlpha_sorted,tabcs.order] = sort(vecAlpha);
            alphai = vecAlpha_sorted(2); %rhoi'>rhoj'
            alphaj = vecAlpha_sorted(1);

            %cii i.e. c (max) such that rhohati=rhoi'
            %already given by tabcn

            %cjij i.e. c such that rhohatij = rhoj' (%diff from cjixj i.e. c
            %such that rhohatixj = rhoj')
%             if nSIS==2 && nSIR==0 && nSICAT==0
%                 [U120,P120,cjij] = U12_SIS2(alphaj,beta1,beta2,gamma1,gamma2,s1,s2,b,mu,0);
%             elseif nSIS==1 && nSIR==1 && nSICAT==0
%                 [U120,P120,cjij] = U12_SIS2(alphaj,beta2,beta1,gamma2,gamma1,s2,s1,b,mu,0);
%             end

            %switch from one zone to another (Uij to Ui)
            %csijtoi = findCswitch(tab.rhohat,tab.one(:,ii),tab.nm,vecC,pres);

            %costs of disease elimination
            tabcs.cs2dis = max(vecC(tab.two.rhohat>=alphai|abs(tab.two.rhohat-alphai)<100*eps));%+step/2; %when 2 diseases are eliminated
            tabcs.cs1dis = max(vecC(tab.two.rhohat>=alphaj|abs(tab.two.rhohat-alphaj)<100*eps));%+step/2; %when only one disease is eliminated
    
        elseif (N==3)
            %alphai>=alphaj>=alphak
            [vecAlpha_sorted,tabcs.order] = sort(vecAlpha); 
            alphak = vecAlpha_sorted(1); %k : first disease eliminated
            alphaj = vecAlpha_sorted(2);
            alphai = vecAlpha_sorted(3);
            tabcs.cs1dis = max(vecC(tab.rhohat.rhohat>=alphak));
            tabcs.cs2dis = max(vecC(tab.rhohat.rhohat>=alphaj));
            tabcs.cs3dis = max(vecC(tab.rhohat.rhohat>=alphai));

    %         %cii
    %         vecCii = [c11,c22,c33]; %listCii = {' $c_1^1=$',' $c_2^2=$',' $c_3^3=$'}; 
    %         cii = vecCii(ii);
    %         %cjij
    %         vecCjij(:,:,1) = [0,c112,c113;c112,0,0;c113,0,0];
    %         vecCjij(:,:,2) = [0,c212,0;c212,0,c223;0,c223,0];
    %         vecCjij(:,:,3) = [0,0,c313;0,0,c323;c313,c323,0];
    %         cjij = vecCjij(ii,jj,jj);
    %         %listCjij...
    %         %ck123
    %         vecCk123 = [c1123,c2123,c3123];
    %         ck123 = vecCk123(kk);
    % 
    %         %cswitch
    %         %Il doit y avoir 2 switches dans ce modèle. Un entre rhoihat et rhoijhat et
    %         %un entre rhoijhat et rhoijkhat.
    %         tabRhohats = [vecRhohat1;vecRhohat2;vecRhohat3;vecRhohat12;vecRhohat13;vecRhohat23];
    %         vecRhohati = tabRhohats(ii,:);
    %         vecRhohatij = tabRhohats(ii+jj+1,:); %voir les notes du 23/02/22 (fin de page du verso)
    % 
    %         cs123toij = max(vecC(vecRhohat==vecRhohatij & vecRhohatij~=0));%ok (affiner le c)
    %         csijtoi = max(vecC(vecRhohat==vecRhohati & vecRhohati~=0));%ok (affiner le c)
    %         cs123toi = max(vecC(vecRhohat==vecRhohatij & vecRhohatij~=0));%ok (affiner le c)
    %         tabcs.cs3dis = max(vecC(vecRhohat==maxalpha)); %when 3 diseases are eliminated
    %         tabcs.cs2dis = max(vecC(vecRhohat<maxalpha & vecRhohat>=alphak & vecRhohat>=alphaj))+pasC/2; %when 2 diseases are eliminated
    %         tabcs.cs1dis = max(vecC(vecRhohat>=minalpha))+pasC/2; %when only one disease is eliminated
        elseif (N==4)
            [vecAlpha_sorted,tabcs.order] = sort(vecAlpha); 
            alphal = vecAlpha_sorted(1); %k : first disease eliminated
            alphak = vecAlpha_sorted(2);
            alphaj = vecAlpha_sorted(3);
            alphai = vecAlpha_sorted(4);
            tabcs.cs1dis = max(vecC(tab.rhohat.rhohat>=alphal));
            tabcs.cs2dis = max(vecC(tab.rhohat.rhohat>=alphak));
            tabcs.cs3dis = max(vecC(tab.rhohat.rhohat>=alphaj));
            tabcs.cs4dis = max(vecC(tab.rhohat.rhohat>=alphai));
        else
            error('the case where N>4 has not been coded yet')
        end

        tabcs.c0     = min(vecC(tab.rhohat.rhohat==0));
    end
end

