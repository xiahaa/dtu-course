function x = lssolver(A,b)
%% solve Ax = b

    solverType = 'svd';%% svd is faster

    if norm(b) < 0.1 || strcmp(solverType,'svd') == 1
        %% near zero, solve Abar xbar = 0
        Abar = [A -b.*ones(size(A,1),1)];
        [U,S,V] = svd(Abar);
        xbar = V(:,end);
        xbar = xbar ./ xbar(end);%% normalization
        x = xbar(1:end-1);
    elseif strcmp(solverType,'qr') == 1
        %% use QR decomposition
        [Q,R] = qr(A);
        x = back_substitution(Q,R,b);
    end
end

function x = back_substitution(Q,R,b)
    y = Q'*b;
    [m,n] = size(R);
    x = zeros(n,1);
    x(end) = y(n)./R(n,n);
    for i = n-1:-1:1
        x(i) = 1/R(i,i)*(y(i)-R(i,i+1:end)*x(i+1:end));
    end
end