function V2 = calcSmoothnessPrior(ds, beta)
    ds = double(ds);
    [m,n]=size(ds);
    V2 = 0;
    ds = double(ds);
    % row sweeping
    for i = 1:m-1   
        V2 = V2+sum(abs(ds(i,:)-ds(i+1,:))>1e-3);
    end
    % col sweeping
    for i = 1:n-1   
        V2 = V2+sum(abs(ds(:,i)-ds(:,i+1))>1e-3);
    end
    V2 = V2.*beta;
end

