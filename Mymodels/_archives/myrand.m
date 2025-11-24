function [res] = myrand(n, mymin,mymax, k)
    % myrand(n, mymin,mymax, k)
    % Cette fonction tire n nombre(s) aleatoire(s) entre mymin et mymax avec une precision de
    % precision de k chiffres après la virgule.
    
    res = round(rand(1,n)*(mymax-mymin) + mymin,k);
end

