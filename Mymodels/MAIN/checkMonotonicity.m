function [res] = checkMonotonicity(arr,direction,tol)
% arr is an array
% we check the direction at each column
tol=abs(tol);
diffArray   = arr(2:end,:) - arr(1:(end-1),:);
punctualInc = sum(diffArray>=-tol,1);
punctualDec = sum(diffArray<=tol,1);

if strcmp(direction,'increasing')
    res = (punctualInc==size(arr,1)-1); %if 1:increasing
elseif strcmp(direction, 'decreasing')
    res = (punctualDec==size(arr,1)-1); %if 1:decreasing
else
    warning(["direction input has not been correctly addressed, choose 'increasing' or 'decreasing'"])
    res=NaN;
end


end