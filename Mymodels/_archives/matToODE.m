function [C,dC,eqn] = matToODE(nbCompartments,M)

    C = sym('X',[1 nbCompartments]);
    syms b;
    
    dC = sym('dX',[1 nbCompartments]);

    for i = 1:nbCompartments
         dC(i) = sum(M(i,:).*C);
    end  
    dC(1) = dC(1)+b;
 
    dCbis = sym('dX',[1 nbCompartments]);
    eqn = dCbis==dC;
    
end