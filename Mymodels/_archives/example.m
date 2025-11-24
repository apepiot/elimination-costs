n=10
greyred = subdivisedColormap([[0.,0.,0.];[1,1,1];[0.85,0.325,0.098]],n, 'quad'); %2^n+1

N=2^n+1;
M = 0.1*[1 2 3 ; -1 -0.5 ,8 ; 0 -1 -2]

x=[1 2 3; 1 2 3;1 2 3];
surf1 = surf(x,x',M, 'FaceColor','interp');

colormap(greyred)
view(2)
minM = min(min(M))
maxM = max(max(M))

%greyred2 = greyred(2^n*0.5*(1-abs(minM)/maxM):end,:)
truncated_len = (minM - (-maxM)) / (2 * maxM);
truncation_ix = round(N * truncated_len);
greyred2 = greyred(truncation_ix:end, :);
colormap(greyred2)
colorbar
