function [d_satclk,clkerr] = sat_clock_error(sp3SatInfo, id, ts)
    c = 299792458;% m / s;
    %% find the corresponding sp3 indices
    tsp3 = sp3SatInfo{id}(:,1);
    clkrec = sp3SatInfo{id}(:,5);
    %% interpolation the clkrec
    clkerr = interp1(tsp3,clkrec,ts,'linear');
    d_satclk = c * clkerr*1e-6;
end