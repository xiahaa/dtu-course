function V1 = calcLikelihood(miuf, ds, fs, alpha)
%compute likelihood cost
    ds = double(ds);
    V1 = alpha*sum((miuf(fs(:))-ds(:)).^2);
end