function [cnn,tabP,msg] = find_Cnn_kit_2(paramTab,mu,b,paramRho,f,kit,paramsC,vecAlpha,afficherOutput,log_path,paramSolver,aprioriCnn,ampl_models_dir)
%This function finds the cost of elimination of each disease of the kit
%when the kit is used
addTime=5;
N = length(kit);
tolP0 = paramSolver.tolP0;
tolC  = paramsC.tolC;
k = indexKit(kit);
msg = [];
supBound = paramsC.sup;
infBound = -abs(paramsC.inf)*5;
iterMax = paramsC.iterMax;

includeHIV = ~vecAlpha.elim_h;
if includeHIV && ~contains(k,'h')
    mod = [k,'_','h',k];
    alpha_max_mod = getAlphaMax(vecAlpha,mod,N);
else
    mod = [k,'_',k];
    alpha_max_mod = getAlphaMax(vecAlpha,mod,N);
end


%% Preliminary verification
%Chercher rhohat_k a-t-il un sens ? Non, si pour rho_k=0, les infections
%sont déjà éliminées

% [paramTab,paramRho] = updateParamRho(paramTab,paramRho,kit,0); %replace paramRho.(['rho_',k]) by 0
% [ES,~,~,~,~,~,msg_err] = ...
%     P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,f,'knitroampl',20);
% msg = [msg,msg_err];
%disp('#####################################')
%disp('ajouter contraintes dans FindArgmax')
%disp('####################################')

fprintf('  -----------------------------------------------------------------------------------------------------------------\n')
fprintf(['      c_',k,'  |   * P_',k,' *  |     P_h    |     P_s    |     P_c    |     P_g    |   rhohat  |   status   |  time\n']);

%----------------------------------%
bnd_sup_c = supBound;
bnd_inf_c = paramsC.inf;
tabP.c       = [];
tabP.rhohat  = [];
tabP.U       = [];

tabP.P_k     = [];
tabP.P_h     = [];
tabP.P_s     = [];
tabP.P_c     = [];
tabP.P_g     = [];
tabP.status  = [];

err = inf;

j = 0;
ampl_c = 0;

% if length(k)==1
%     newmod=0;
% else
newmod=0;
% end

inf_bnd_not_found=1;
sup_bnd_not_found=1;

%% -1. Computing around the costs from apriori
if ~isempty(aprioriCnn)
    fprintf('  ----------------------------------- A priori -------------------------------------------------------------------\n')
    aprioriCnn_simplied = unique(round(aprioriCnn,4));
    aprioriCnn_up_and_down = sort([aprioriCnn_simplied+tolC*0.9, aprioriCnn_simplied-tolC*0.95]);
    i=1; j=0;
    while i<=length(aprioriCnn_up_and_down) && j<=iterMax
        tic;
        c = aprioriCnn_up_and_down(i);
        [rhohat_all,Cval,ESval,msg_solver,ampl_c] = findRhohat_kit_v9(paramTab,mu,b,f,paramRho,kit,c,vecAlpha,afficherOutput,j,log_path,paramSolver,ampl_c,newmod,ampl_models_dir);
        rhohat = rhohat_all.tot;
        ES = ESval.tot;
        P_k = P_kit_v2(ES,kit);           tabP.P_k = [tabP.P_k;P_k];
        P_h = P_kit_v2(ES,{'HIV'});       tabP.P_h = [tabP.P_h;P_h];
        P_s = P_kit_v2(ES,{'syphilis'});  tabP.P_s = [tabP.P_s;P_s];
        P_c = P_kit_v2(ES,{'Ct'});        tabP.P_c = [tabP.P_c;P_c];
        P_g = P_kit_v2(ES,{'Ng'});        tabP.P_g = [tabP.P_g;P_g];

        tabP.rhohat = [tabP.rhohat;rhohat];
        tabP.U      = [tabP.U;-Cval.tot];
        tabP.c      = [tabP.c;c];
        tabP.status = [tabP.status;str2double(msg_solver)];
        tabP = sortTable(tabP);

        fprintf(['    %+2.4f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.4f  |  %3.0f  |  %4.0f  \n'],c,P_k,P_h,P_s,P_c,P_g,rhohat_all.tot,str2double(msg_solver),toc);
        if str2double(msg_solver)==0
            i=i+1;
        else
            paramSolver.timeSolver = min(paramSolver.timeSolver+addTime,120);
        end
        j=j+1;
    end
    c_Pk_0            = tabP.c((tabP.P_k-tolP0)<tolP0 & tabP.status==0);
    inf_bnd_not_found = isempty(c_Pk_0);
    if ~inf_bnd_not_found
        bnd_inf_c = max(c_Pk_0);
    end 
%     for j=1:length(tabP.c)
%         sup_bnd_not_found = 0;
%         for dis=kit
%             inf = indexKit(dis);
%             P_i = tabP.(['P_',inf]);
%             if (P_i(j)-tolP0)<tolP0         %if one infection is eliminated
%                 sup_bnd_not_found = 1;
%
%             end
%         end
%     end
%     if ~sup_bnd_not_found
%         bnd_sup_c = max(c_Pk_0);
%     end
if str2double(msg_solver)==0
    msg = [msg, 'a0'];
else
    msg = [msg, 'a1'];
end
end


msg = [msg,'-'];

%% 0. Finding the interval where the solutions are supposed to be
%0.a) inf bound
fprintf('  ----------------------------------- Borne inf -------------------------------------------------------------------\n')
while inf_bnd_not_found && j<iterMax && bnd_inf_c>=infBound
    tic;
    [rhohat,Cval,ESval,msg_solver,ampl_c] = findRhohat_kit_v9(paramTab,mu,b,f,paramRho,kit,bnd_inf_c,vecAlpha,afficherOutput,j,log_path,paramSolver,ampl_c,newmod,ampl_models_dir);
    rho_inf_bnd = rhohat.tot;
    ES = ESval.tot;
    %[paramTab,paramRho] = updateParamRho(paramTab,paramRho,kit,rho_inf_bnd);
    %paramRhobis = paramRho;
    %paramRhobis.rho_hs = rhohat.hs_hs;
    %[P,ES,status] = P_mod_v7(paramTab,paramRhobis,b,mu,f,mod,verbose,paramSolver,log_path,ampl_models_dir);
    
    P_k = P_kit_v2(ES,kit);            tabP.P_k = [tabP.P_k;P_k];
    P_h = P_kit_v2(ES,{'HIV'});        tabP.P_h = [tabP.P_h;P_h];
    P_s = P_kit_v2(ES,{'syphilis'});   tabP.P_s = [tabP.P_s;P_s];
    P_c = P_kit_v2(ES,{'Ct'});         tabP.P_c = [tabP.P_c;P_c];
    P_g = P_kit_v2(ES,{'Ng'});         tabP.P_g = [tabP.P_g;P_g];
    fprintf(['    %+2.4f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.4f  |  %3.0f  |  %4.0f  \n'],bnd_inf_c,P_k,P_h,P_s,P_c,P_g,rhohat.tot,str2double(msg_solver),toc);
    
    %Update vec_c and vec_rho_hat
    tabP.c = [tabP.c;bnd_inf_c];
    tabP.rhohat = [tabP.rhohat;rho_inf_bnd];
    tabP.U = [tabP.U;-Cval.tot];
    tabP.status  = [tabP.status;str2double(msg_solver)];
    tabP = sortTable(tabP);
    
    %If the low boundary has not been found, we look for another one
    %c_Pk_0            = tabP.c((tabP.P_k-tolP0)<tolP0 & tabP.status==0);
    %inf_bnd_not_found = isempty(c_Pk_0);
    inf_bnd_not_found = ~(((P_k-tolP0)<tolP0 || abs(rhohat.tot-alpha_max_mod)<tolP0) & str2double(msg_solver)==0);
    
    if inf_bnd_not_found
        bnd_inf_c = bnd_inf_c-1;
    end
    
    if ~isequal(msg_solver,'0')
        paramSolver.timeSolver = min(paramSolver.timeSolver+addTime,120);
    end
    
    j=j+1;
end

if inf_bnd_not_found
    warning(['la borne inf de c_',k,' n a pas ete trouvee'])
    msg = [msg, 'i1'];
else
    msg = [msg, 'i0'];
end
msg = [msg,'-'];

fprintf('  --------------------------------------------- Borne sup ---------------------------------------------------------\n')

%0.b) sup bound
j=0;
while sup_bnd_not_found && j<iterMax %&& (bnd_sup_c<supBound || j==0)
    tic;
    [rhohat,Cval,ESval,msg_solver,ampl_c] = findRhohat_kit_v9(paramTab,mu,b,f,paramRho,kit,bnd_sup_c,vecAlpha,afficherOutput,j,log_path,paramSolver,ampl_c,newmod,ampl_models_dir);
    rho_sup_bnd = rhohat.tot;
    ES = ESval.tot;
    
    %[paramTab,paramRho] = updateParamRho(paramTab,paramRho,kit,rho_sup_bnd);
    %[ES,~,~,~,~,~] = ...
    %    P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,f,'knitroampl',20);
    P_k = P_kit_v2(ES,kit);            tabP.P_k = [tabP.P_k;P_k];
    P_h = P_kit_v2(ES,{'HIV'});        tabP.P_h = [tabP.P_h;P_h];
    P_s = P_kit_v2(ES,{'syphilis'});   tabP.P_s = [tabP.P_s;P_s];
    P_c = P_kit_v2(ES,{'Ct'});         tabP.P_c = [tabP.P_c;P_c];
    P_g = P_kit_v2(ES,{'Ng'});         tabP.P_g = [tabP.P_g;P_g];
    
    %Update vec_c and vec_rho_hat
    tabP.c = [tabP.c;bnd_sup_c];
    tabP.rhohat = [tabP.rhohat;rho_sup_bnd];
    tabP.U = [tabP.U;-Cval.tot];
    tabP.status  = [tabP.status;str2double(msg_solver)];
    tabP = sortTable(tabP);
    
    % the sup boundary needs to satisfy the condition, where all infections
    % are endemic (i.e. P_i > 0)
    %If the high boundary has not been found, we look for another one
    if isequal(msg_solver,'0')
        sup_bnd_not_found = 0;
        for dis=kit
            P_i = P_kit_v2(ES,dis);
            if (P_i-tolP0)<tolP0     %if one infection is eliminated
                sup_bnd_not_found = 1;
            end
        end
    else
        sup_bnd_not_found = 1;
        paramSolver.timeSolver = min(paramSolver.timeSolver+addTime,120);
        bnd_sup_c = min(bnd_sup_c+1/10,1);
    end
    
    fprintf(['    %+2.4f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.4f  |  %3.0f  |  %4.0f  \n'],bnd_sup_c,P_k,P_h,P_s,P_c,P_g,rhohat.tot,str2double(msg_solver),toc);
    j = j+1;
end
if sup_bnd_not_found
    warning(['la borne sup de c_',k,' n a pas ete trouvee'])
    msg = [msg, 's1'];
else
    msg = [msg, 's0'];
end

msg = [msg,'-'];

fprintf(['  ----------------------------------------- Recherche de c_',k,' ----------------------------------------------------\n'])

if ~sup_bnd_not_found && ~inf_bnd_not_found
    c_Pk_0 = tabP.c((abs(tabP.P_k-tolP0)<tolP0 | abs(tabP.rhohat-alpha_max_mod)<tolP0) & tabP.status==0);
    c_Pk_not_0 = tabP.c((abs(tabP.P_k-tolP0)>=tolP0 & abs(tabP.rhohat-alpha_max_mod)>tolP0 ) & tabP.status==0);
    bnd_inf_c  = c_Pk_0(end);
    bnd_sup_c  = c_Pk_not_0(1);
    c = (bnd_sup_c+bnd_inf_c)/2;
    j = 0;
    %% 1. First objective: finding c such that P_kit(rho_hat_kit)=0
    while (err>tolC && j<iterMax)
        tic;
        %Updating vec_rho_hat and vec_c
        [rhohat_all,Cval,ESval,msg_solver,ampl_c] = findRhohat_kit_v9(paramTab,mu,b,f,paramRho,kit,c,vecAlpha,afficherOutput,j,log_path,paramSolver,ampl_c,newmod,ampl_models_dir);
        rhohat = rhohat_all.tot;
        ES = ESval.tot;
        
        %Condition to continue looking for new rho_hat
        %if the prevalence_kit(rho_kit) different from zero, we continue
        %%[paramTab,paramRho] = updateParamRho(paramTab,paramRho,kit,rhohat);
        %%[ES,~,~,~,~,~] = ...
        %%    P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,f,'knitroampl',20);
        P_k = P_kit_v2(ES,kit);           tabP.P_k = [tabP.P_k;P_k];
        P_h = P_kit_v2(ES,{'HIV'});       tabP.P_h = [tabP.P_h;P_h];
        P_s = P_kit_v2(ES,{'syphilis'});  tabP.P_s = [tabP.P_s;P_s];
        P_c = P_kit_v2(ES,{'Ct'});        tabP.P_c = [tabP.P_c;P_c];
        P_g = P_kit_v2(ES,{'Ng'});        tabP.P_g = [tabP.P_g;P_g];
        
        tabP.rhohat  = [tabP.rhohat;rhohat];
        tabP.U       = [tabP.U;-Cval.tot];
        tabP.c       = [tabP.c;c];
        tabP.status  = [tabP.status;str2double(msg_solver)];
        tabP = sortTable(tabP);
        
        fprintf(['    %+2.4f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.4f  |  %3.0f  |  %4.0f  \n'],c,P_k,P_h,P_s,P_c,P_g,rhohat_all.tot,str2double(msg_solver),toc);
        
        err = (bnd_sup_c-bnd_inf_c);
        
        if ~isnan(P_k) && isequal(msg_solver,'0')
            if abs(P_k-tolP0)<=tolP0 || abs(rhohat_all.tot-alpha_max_mod)<tolP0 % abs(rhohat_all.tot-vecAlpha.(mod))<=tolP0
                bnd_inf_c = c;
                c = (c+bnd_sup_c)/2;
            else
                bnd_sup_c = c;
                c = (bnd_inf_c+c)/2;
            end
        elseif ~isequal(msg_solver,'0')
            %c=0.95*c;
            paramSolver.timeSolver = min(paramSolver.timeSolver+addTime,120);
        else
            %continuer avec le même c et relancer
            %c = c*(1+0.001);
        end
        
        j   = j+1;
        tabP = sortTable(tabP);
    end
    cnn.kit = c; %cost of elimination of all the infections in the kit, i.e. cost of elimination of the last infection
    fprintf(['  c_%s=%2.4f \n'],k,c)
     if abs(c-infBound)<tolC
        msg = [msg, 'ik1'];      %la borne inf trouvee pour c est le minimum possible.
    else
        msg = [msg, 'ik0'];
    end   
else %la borne superieure en 1, a donné l'elimination des infections
    cnn.kit = Inf;
    fprintf(['  c_%s=%2.4f \n'],k,cnn.kit)
end


msg = [msg,'-'];

%% 2. Second objective: finding for every infection i, c such that P_i(rho_hat_kit)=0
for i=1:N
    dis = indexKit(kit{i});
    fprintf('  -------------------------------------------------------------------------------------------------------------------\n')
    fprintf(['      c_',k,'  |     P_',k,'    | *  P_%c  * ||     P_h    |     P_s    |     P_c    |     P_g    |   rhohat  |   status  |  time \n'],dis);
    
    Pi = tabP.(['P_',dis]);
    %limit between Pi>0 and Pi=0;
    j = 0;
    c_Pi_0    = tabP.c(abs(Pi-tolP0)<=tolP0 & tabP.status==0);
    if isempty(c_Pi_0)
        bnd_inf_c = infBound;
    else
        bnd_inf_c = c_Pi_0(end);    %Pi=0
    end
    c_Pi_not0 = tabP.c(abs(Pi-tolP0)>tolP0 & bnd_inf_c<=tabP.c & tabP.status==0);
    if isempty(c_Pi_not0)
        bnd_sup_c = supBound;
    else
        bnd_sup_c = c_Pi_not0(1);   %Pi>0
    end
    err = (bnd_sup_c-bnd_inf_c);
    %On affine pour trouver le cout d'elimination de l'infection i
    c = (bnd_inf_c+bnd_sup_c)/2;

    if err<0 %problème de monotonie de \Pi_i(\hat\rho(c))
        % a voir
    end
    
    while err>tolC && j<iterMax
        tic;
        [rhohat_all,Cval,ESval,msg_solver,ampl_c] = findRhohat_kit_v9(paramTab,mu,b,f,paramRho,kit,c,vecAlpha,afficherOutput,j,log_path,paramSolver,ampl_c,newmod,ampl_models_dir);
        rhohat = rhohat_all.tot;
        ES = ESval.tot;
        
        tabP.rhohat = [tabP.rhohat;rhohat];
        tabP.c      = [tabP.c;c];
        tabP.U      = [tabP.U;-Cval.tot];
        tabP.status  = [tabP.status;str2double(msg_solver)];
        
        %[paramTab,paramRho] = updateParamRho(paramTab,paramRho,kit,rhohat);
        %[ES,~,~,~,~,~] = ...
        %    P1234_SICTPSEIIISSEIIS2_v5(paramTab{3},paramTab{4},paramTab{1},paramTab{2},mu,b,paramRho,f,'knitroampl',20);
        
        P_k = P_kit_v2(ES,kit);            tabP.P_k = [tabP.P_k;P_k];
        P_h = P_kit_v2(ES,{'HIV'});        tabP.P_h = [tabP.P_h;P_h];
        P_s = P_kit_v2(ES,{'syphilis'});   tabP.P_s = [tabP.P_s;P_s];
        P_c = P_kit_v2(ES,{'Ct'});         tabP.P_c = [tabP.P_c;P_c];
        P_g = P_kit_v2(ES,{'Ng'});         tabP.P_g = [tabP.P_g;P_g];
        Pi  = tabP.(['P_',dis]);
        
        %fprintf(['   %+2.4f  |  %1.6f  |  %1.6f  ||  %1.6f  |  %1.6f  |  %1.6f  |  %1.6f  |  %3.0f   \n'],c,P_k,Pi(end),P_h,P_s,P_c,P_g,str2double(msg));
        fprintf(['    %+2.4f  |  %1.6f  |  %1.6f  ||  %1.6f  |  %1.6f  |  %1.6f  |  %1.6f  |  %1.4f  |  %3.0f  |  %4.0f  \n'],...
            c,P_k,Pi(end),P_h,P_s,P_c,P_g,rhohat_all.tot,str2double(msg_solver),toc);
        
        tabP = sortTable(tabP);
        Pi = tabP.(['P_',dis]);
        
        c_Pi_0    = tabP.c((Pi-tolP0)<tolP0);
        c_Pi_not0 = tabP.c((Pi-tolP0)>=tolP0);
        
        if isempty(c_Pi_0)
            bnd_inf_c = infBound;
        else
            bnd_inf_c = c_Pi_0(end);    %Pi=0
        end
        if isempty(c_Pi_not0)
            bnd_sup_c = supBound;
        else
            bnd_sup_c = c_Pi_not0(1);   %Pi>0
        end
        
        err = (bnd_sup_c-bnd_inf_c);
        
        if err<0
            c_Pi_0    = tabP.c((Pi-tolP0*10)<tolP0*10);
            c_Pi_not0 = tabP.c((Pi-tolP0*10)>=tolP0*10);
            bnd_inf_c = c_Pi_0(end);
            bnd_sup_c = c_Pi_not0(1);   %Pi>0
            err = (bnd_sup_c-bnd_inf_c);
        end
            
        
        if isequal(msg_solver,'0')
            c = (bnd_inf_c+bnd_sup_c)/2;
        else
            %c=c%0.95*c; %seed is different
            paramSolver.timeSolver = min(paramSolver.timeSolver+addTime,120);
        end
        j = j+1;
    end
    
    if abs(err)<=tolC
        msg = [msg, 'i',dis,'0-'];      %la borne inf trouvee pour c est le minimum possible.
    else
        msg = [msg, 'i',dis,'1-'];
    end
    
    if sup_bnd_not_found
        cnn.(kit{i}) = Inf;
    else
        cnn.(kit{i}) = c;
    end
    fprintf('  c_%s=%2.4f \n',dis,c)
end
ampl_c.close();
fprintf('  -------------------------------------------------------------------------------------------------------------------\n')
end

