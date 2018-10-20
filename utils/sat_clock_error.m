function d_satclk = sat_clock_error(sp3t, sp3cr, PRN, ts)
    c = 299792458;% m / s;
    
    %% find the corresponding sp3 indices
    id = sp3(:,1) == PRN;
    tsp3 = sp3(id,3);
    clkrec = sp3(id,5);
    
    
    
    d_satclk = c * clk_rec*1e-6;
end