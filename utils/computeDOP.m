function QDOP = computeDOP(A,P,var_prior)
    if isempty(P) == 1
        QDOP = inv(A'*A);
    else
        QDOP = inv(A'*P*A);
    end
    
    QDOP = QDOP/var_prior;
    
end