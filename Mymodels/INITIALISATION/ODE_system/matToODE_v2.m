function [X,dX,eqn,dXright] = matToODE_v2(n,M,B)
    syms X [n 1]  
    syms dXright [n 1] 
    
    dXright = M*X+B;
 
    syms dX [n 1]%('dX',[1 n]);
    eqn = dX==dXright;%dCbis==dC;
end