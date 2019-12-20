function PDOP = PDOP_calc(M)
    Minv = inv(M);
    PDOP = sqrt(trace(Minv(1:3,1:3)));
end