function [newarray] = subdivisedColormap(array, n, type) %m x 3
m = size(array,1);
if contains('quad', type)
    if m==1
        %disp('m=1')
        newarray = array;
        %disp(newarray)
    elseif m==2 && n==1
        newarray = [array(1,:);array(1,:)+(array(2,:)-array(1,:))./2;array(2,:)];
        %disp(newarray)
    elseif m==2 && n>1
        col1 = array(1,:); col2 = array(2,:);
        mid = (col2-col1)/2+col1;
        col_left  = subdivisedColormap([col1;mid], n-1, type);
        col_right = subdivisedColormap([mid;col2], n-1, type);
        newarray  = [col_left;col_right(2:end,:)];
        %disp(newarray)
    else
        newarray = array(1,:);
        for i=1:(m-1)
            col_left  = array(i,:);
            col_right = array(i+1,:);
            add_colors = subdivisedColormap([col_left;col_right], n-1, type);
            newarray = [newarray;add_colors(2:end,:)];
        end
        %disp(n)
        %disp(newarray)
    end        
    
elseif strcmp('linear',type)
    newarray = array(1,:);
    for i=1:(m-1)
        col_left  = array(i,:);
        col_right = array(i+1,:);
        new_r = linspace(col_left(1),col_right(1),n);
        new_g = linspace(col_left(2),col_right(2),n);
        new_b = linspace(col_left(3),col_right(3),n);
        new_cols = [new_r;new_g;new_b]';
        newarray = [newarray;new_cols(2:end,:)];
    end
elseif strcmp('log',type)
    %highlights more colors
    if n==1 || m==1
        newarray = array;
    elseif n>=2 && m==2
        array(array==0) = 0.0001;
        newarray = zeros(n,3);
        newarray(:,1) = exp(linspace(log(array(1,1)),log(array(2,1)),n));
        newarray(:,2) = exp(linspace(log(array(1,2)),log(array(2,2)),n));
        newarray(:,3) = exp(linspace(log(array(1,3)),log(array(2,3)),n));
    else
        array(array==0) = 0.0001;
        array_log  = log(array);
        array_log_lin = subdivisedColormap(array_log,n,'linear');
        newarray   = exp(array_log_lin);
    end  
elseif strcmp('log10',type)
    %highlights more colors
    if n==1 || m==1
        newarray = array;
    elseif n>=2 && size(array,1)==2
        array(array==0) = 0.0001;
        newarray = zeros(n,3);
        newarray(:,1) = 8.^(linspace(log10(array(1,1)),log10(array(2,1)),n));
        newarray(:,2) = 8.^(linspace(log10(array(1,2)),log10(array(2,2)),n));
        newarray(:,3) = 8.^(linspace(log10(array(1,3)),log10(array(2,3)),n));
    else
        warning('a faire')
    end
elseif strcmp('straigthen')
    
end
end

