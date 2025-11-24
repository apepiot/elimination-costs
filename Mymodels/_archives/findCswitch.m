function [cs] = findCswitch(rhohat, rhohatj,rhohatij,vecC,pres)
    %finding c where rhohat switch from one curve to another
    
    % Finding cs1 and cs2 such that cs in [cs1,cs2]
    % rhohatj rhohatj rhohatj rhohatij rhohatij rhohatij
    % some c  some c  cs1     cs2      some c   some c
    
    cs2 = min(vecC(rhohat ~= rhohatj));
    cs1 = max(vecC(rhohat~=rhohatij)) ;
    
    err = abs(cs1-cs2);
    
    if(err==0) %then cs1=cs2=cs
        cs2;
    else 
        %if(err>pres)
            %refaire le processus en precisant cs
        %end
        cs = (cs2+cs1)/2;
    end
    
    % cs = (cs2+cs1)/2;
    % while abs((cs1-cs2)/(vecCmin-vecCmax))>pres then affiner cs

    
end

