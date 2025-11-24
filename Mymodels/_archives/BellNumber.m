function [m] = BellNumber(n)
%% this function gives the number of combination of n disease
%%e.g. {A,B,C} gives 5 combinations : A,B,C ; A+B,C ; A+C,B ; A,B+C ; A+B+C
    if n==0
        m = 1;
    else
        m = 0;
        for k=0:(n-1)    
            m = m + nchoosek(n-1,k)*BellNumber(k);
        end
    end
end