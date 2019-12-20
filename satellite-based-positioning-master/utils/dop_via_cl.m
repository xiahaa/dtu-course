function dop = dop_via_cl(H)
    M = H'*H;
    h1 = trace(M);
    h2 = trace(M*M);
    h3 = trace(M*M*M);
    h4 = det(M);
    dop = sqrt((0.5*h1*h1*h1-1.5*h1*h2+h3)/(3*h4));
end