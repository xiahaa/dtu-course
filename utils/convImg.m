function If = convImg(I, kernel)
% just to try convolution.
    [m,n] = size(kernel);
    [M,N] = size(I);
    % 
    m2 = floor(m*0.5);
    n2 = floor(n*0.5);
    Iv = zeros(M+m2*2,N+n2*2);
    Iv(m2+1:m2+M,n2+1:n2+N) = I;
    Iv(1:m2,:) = [repmat(I(1,1),m2, n2) ones(m2, 1)*I(1,1:N) repmat(I(1,end),m2, n2)];
    Iv(m2+M+1:m2+M+m2,:) = [repmat(I(end,1),m2, n2) ones(m2, 1)*I(end,1:N) repmat(I(end,end),m2, n2)];
    Iv(1:M, 1:n2) = repmat(I(:,1), 1, n2);
    Iv(1:M, n2+N+1:1:n2+N+n2) = repmat(I(:,end), 1, n2);
    % staightforward way
    [mx, my] = meshgrid(-m2:1:m2, -n2:1:n2);
    piv = Iv(:);
    fkernel = fliplr(kernel);
    pkernel = fkernel(:);
    If = zeros(size(I));
    for i = m2+1:m2+M        
        for j = n2+1:n2+N
            sub = [mx+i, my+j];
            ind = sub2ind(size(Iv), sub(:,1),sub(:,2));
            eval = piv(ind)'*pkernel;
            If(i-m2,j-n2) = eval;
        end
    end
end