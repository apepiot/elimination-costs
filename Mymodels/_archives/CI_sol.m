function [F] = CI_sol(Y, S_h,I_h,C_h,P_h,Ip_h,Cp_h,T_h,...
                         S_s,E_s,I1_s,I2_s,I3_s,...
                         S_c,E_c,IA_c,IS_c,...
                         S_g,E_g,IA_g,IS_g,...
                         N)
     
    %Contraintes HIV                
    F(1) = sum(Y(1:7:554)) - S_h;
    F(2) = sum(Y(2:7:555)) - I_h;
    F(3) = sum(Y(3:7:556)) - C_h;
    F(4) = sum(Y(4:7:557)) - P_h;
    F(5) = sum(Y(5:7:558)) - Ip_h;
    F(6) = sum(Y(6:7:559)) - Cp_h;
    F(7) = sum(Y(7:7:560)) - T_h;
    
    %Contraintes syphilis
    F(8)  = sum(Y(reshape((repmat(1:7,560/(7*5),1)+(0:35:526)')', [1,560/5]))) - S_s;
    F(9)  = sum(Y(reshape((repmat(1:7,560/(7*5),1)+(7:35:533)')', [1,560/5]))) - E_s;
    F(10) = sum(Y(reshape((repmat(1:7,560/(7*5),1)+(14:35:540)')',[1,560/5]))) - I1_s;
    F(11) = sum(Y(reshape((repmat(1:7,560/(7*5),1)+(21:35:547)')',[1,560/5]))) - I2_s;
    F(12) = sum(Y(reshape((repmat(1:7,560/(7*5),1)+(28:35:554)')',[1,560/5]))) - I3_s;
    
    %Contraintes Ct
    F(13) = sum(Y(reshape((repmat(1:35,560/(7*5*4),1)+(0:(7*5*4):421)')',  [1,560/4]))) - S_c;
    F(14) = sum(Y(reshape((repmat(1:35,560/(7*5*4),1)+(35:(7*5*4):456)')', [1,560/4]))) - E_c;
    F(15) = sum(Y(reshape((repmat(1:35,560/(7*5*4),1)+(70:(7*5*4):491)')', [1,560/4]))) - IA_c;
    F(16) = sum(Y(reshape((repmat(1:35,560/(7*5*4),1)+(105:(7*5*4):526)')',[1,560/4]))) - IS_c;
    
    %Contraintes Ng
    F(17) = sum(Y(1:140))   - S_g;
    F(18) = sum(Y(141:280)) - E_g;
    F(19) = sum(Y(281:420)) - IA_g;
    F(20) = sum(Y(421:560)) - IS_g;
    
    F(21) = sum(Y)-N;
                                            
end

