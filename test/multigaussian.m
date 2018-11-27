function prob = multigaussian(x, u, Omega, w)
    n = size(u,1);
    m = size(u,2);
    prob = 0;
    for i = 1:n
        R = Omega(i,:,:);
        R = reshape(R,m,m,1);
        prob = prob + w(i).*1/(sqrt((2*pi)^(m)*det(R)))* ...
                      exp(-0.5*(x-u(i,:)')'*inv(R)*(x-u(i,:)'));
    end
end