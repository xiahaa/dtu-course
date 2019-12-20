function Ps = calc_ps_via_barometric_formula(h)
%calc_ps_via_barometric_formula calculate pressure using barometric formula
%model and constants from
%https://en.wikipedia.org/wiki/Atmospheric_pressure
%Author: xiahaa@space.dtu.dk
    L = 0.0065;%Temperature lapse rate
    P0 = 101325;%sea level standard atmospheric pressure
    T0 = 288.15;%sea level standard temperature
    g = 9.80665;%earth-surface gravitational acceleration
    M = 0.0289644;%Molar mass of dry air
    R0 = 8.31447;%Universal gas constant
    Ps = P0*(1-L*h/T0)^(g*M/R0/L);% pressure in Pa
    Ps = Ps * 0.01;% convert to hPa
end