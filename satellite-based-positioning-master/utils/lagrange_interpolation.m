function yint = lagrange_interpolation(xs,ys,xint)
%% lagrange_interpolation performs langrange interpolation in 1D.
    K = size(xs,1);
    %% 1st, compute basis function
    yint = zeros(numel(xint),1);
    
    for i = 1:numel(xint)
        for j = 1:K
            xj = xs(j);yj = ys(j);
            xo = xs;
            xo(j) = [];
            f = @(x) (yj.*cumprod((repmat(x,size(xo,1),1)-xo)./(repmat(xj,size(xo,1),1)-xo)));
            %% interpolation
            y1 = f(xint(i));
            yint(i) = yint(i) + y1(end);
        end
    end
end

