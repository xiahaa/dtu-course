function ker = boxker(hsize)
% box filter
    M = 2*hsize+1;
    ker = ones(M,M);
    ker = ker ./ sum(ker(:));
end