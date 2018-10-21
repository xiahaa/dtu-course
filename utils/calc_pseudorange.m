function pr = calc_pseudorange(at_pos, sp3SatInfo, id, ts, TECU, recerr)
    
    %% 1 interpolation satllite position
    
    %% 2 calculate real geometric distance: R
    R = 0;
    
    %% 3 compute elevation angle
    el = 0;
    
    %% compute ionosphere delay
    d_iono = iono_delay_first_order_group(el, TECU);
    
    %% compute estimated troposhere delay 
    d_trop = tropo_deley_via_saastamoinen_model(at_pos(1), at_pos(3), el, '1');
    
    %% compute satellite clokc delay
    d_satclk = sat_clock_error(sp3SatInfo, id, ts);
    
    %% compute receiver clock delay
    d_recclk = rec_clock_error(recerr);
    
    pr = R + d_satclk - d_recclk + d_iono + d_trop;
    
end