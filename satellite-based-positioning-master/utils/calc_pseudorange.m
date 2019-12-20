function [visibility, pr, R, d_iono, d_trop, d_satclk, d_recclk, clkerr] = calc_pseudorange(lat, lon, h, sat_pos, sp3SatInfo, id, ts, TECU, recerr)    
    
    consParams = struct('a',6378137.0,'f',1/298.257223563); % some constants
    [xo,yo,zo] = llhtoCartesian(lat, lon, h, consParams);% to ECEF

    dx = sat_pos(1) - xo;
    dy = sat_pos(2) - yo;
    dz = sat_pos(3) - zo;
    %% conversion
    [e,n,u] = WGS842ENU(lat, lon, dx, dy, dz);
    %% compute azimuth and zenith
    [azimuth, zenith, elevation] = calcAzimuthZenithElevation(e,n,u);
    
    if elevation < 5
        visibility = 0;
        pr = inf;
        R = inf;
        d_iono = inf;
        d_trop = inf;
        d_satclk = inf;
        d_recclk = inf;
        clkerr = inf;
        return;
    else
        visibility = 1;
    end

    %% calculate real geometric distance: R
    R = sqrt(dx^2+dy^2+dz^2);
    
    %% compute elevation angle
    el = elevation;
    
    %% compute ionosphere delay
    d_iono = iono_delay_first_order_group(el, TECU);
    
    %% compute estimated troposhere delay 
    d_trop = tropo_delay_via_saastamoinen_model(lat, h, el, '1');
    
    %% compute satellite clokc delay
    [d_satclk,clkerr] = sat_clock_error(sp3SatInfo, id, ts);
    
    %% compute receiver clock delay
    d_recclk = rec_clock_error(recerr);
    
    %% compose pseudo-range
    pr = R - d_satclk + d_recclk + d_iono + d_trop;
end