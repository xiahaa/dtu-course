
function Bint = regularization(a, b, N)
    L1 = spdiags([-2.*ones(N,1) ones(N,1) ones(N,1)],[0,1,-1],N,N);
    L1(1,N) = 1;
    L1(N,1) = 1;
    L1 = L1.*a;
    L2 = spdiags([-1.*ones(N,1) 4.*ones(N,1) -6.*ones(N,1) 4.*ones(N,1) -1.*ones(N,1)], ...
                 [-2,-1,0,1,2],N,N);
    L2(1,N) = 4;L2(1,N-1)=-1;
    L2(2,N) = -1;
    L2(N,1) = 4;L2(N,2) = -1;
    L2(N-1,1) = -1;
    L2 = L2.*b;
    Bint = (eye(N) - (L1+L2));
end
