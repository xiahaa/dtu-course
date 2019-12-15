function [L1,L2] = estab_l1l2(n)
    e = ones(1,n)';
    L2 = spdiags([e,-2*e,e],-1:1,n,n);
    L1 = spdiags([-e,e],[-1,0],n,n);
end