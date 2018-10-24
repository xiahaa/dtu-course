function QDOP = computeDOP(A,P,var_prior)
    if isempty(P) == 1
        QDOP = A'*A;
    else
        QDOP = A'*P*A;
    end
    
    QDOP = QDOP/var_prior;
    
    
end