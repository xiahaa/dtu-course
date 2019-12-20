function d_recclk = rec_clock_error(recerr)

    c = 299792458;% m / s;
    d_recclk = c * recerr*1e-3;

end