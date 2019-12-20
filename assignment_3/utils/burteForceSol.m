function V2 = burteForceSol(ds,beta)
%The most straightforward way to compute the prior cost, only for validatation.
    [m,n]=size(ds);
    V2 = 0;

    ds = double(ds);
    for i = 1:m
        for j = 1:n
%             if i - 1 > 0
%                 V2 = V2 + double(abs(ds(i,j) - ds(i-1,j))>1e-3);
%             end
            if i + 1 <= m
                V2 = V2 + double(abs(ds(i,j) - ds(i+1,j))>1e-3);
            end
%             if j - 1 > 0 
%                 V2 = V2 + double(abs(ds(i,j) - ds(i,j-1))>1e-3);
%             end
            if j + 1 <= n 
                V2 = V2 + double(abs(ds(i,j) - ds(i,j+1))>1e-3);
            end
        end
    end
    V2 = V2.*beta;
end