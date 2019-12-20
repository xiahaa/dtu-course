function ps = curve_smoothing(p, type)
% This file does the curve smoothing by applying the following methods:
% 1. Michael Kass et al. Snakes: Active contour models.                    
% International Journal of Computer Vision, 1(4):321?331, 1988
% 2. Chenyang Xu, Dzung L Pham, and Jerry L Prince. 
% Image segmentation using deformable models. Handbook of medical imaging, 2:129?174, 2000a
% Author: xiahaa@space.dtu.dk
    lambda = 0.5;
    N = size(p,1);
    if type == 1
        L = spdiags([-2.*ones(N,1) ones(N,1) ones(N,1)],[0,1,-1],N,N);
        L(1,N) = 1;
        L(N,1) = 1;
        ps = (p + lambda.*L*p);
    elseif type == 2
        L = spdiags([-2.*ones(N,1) ones(N,1) ones(N,1)],[0,1,-1],N,N);
        L(1,N) = 1;
        L(N,1) = 1;
        ps = (eye(N)-lambda.*L)\p;
    elseif type == 3
        a = 0.5;b = 0.5;
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
        ps = (eye(N) - (L1+L2))\p;
    end
end