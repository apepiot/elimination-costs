function [beta1,beta2,gamma1,gamma2,s1,s2,b,mu,rho] = random_parameters(R1_sup1, R2_sup1)
% random_parameters(R1_supf1, R2_sup1)
% Cette fonction tire un jeu de paramètres (beta1,beta2,gamma1,gamma2,s1,s2,b,mu), 
% sous les contraintes R1>1 ou R2>1 si true.
    
    beta1 = myrand(1,0.1,10,2);
    beta2 = myrand(1,0.1,10,2);

    mu  = myrand(1,0.01,0.9*min(beta1,beta2),3);
    b   = myrand(1,0.01,0.5,2);
    s1  = myrand(1,0.3,1,2);
    s2  = myrand(1,0.3,1,2);

    if(R1_sup1)
        gamma1 = myrand(1,0.001,0.99*(beta1-mu),4);
    else
        gamma1 = myrand(1,0.0001,5,3);
    end

    if(R2_sup1)
        gamma2 = myrand(1,0.00001,0.99*(beta2-mu),4);
    else
        gamma2 = myrand(1,0.000001,5,2);
    end
    
    rho = myrand(1,0.001, 0.9*min((beta1-gamma1-mu), (beta2-gamma2-mu)), 3);


end

