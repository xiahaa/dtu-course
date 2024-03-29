function d_trop = tropo_deley_via_saastamoinen_model(lati, h, el, type)    
%tropo_deley_via_saastamoinen_model estimate troposhere using the saastamoinen_model 
% inputs:
%   lati: latitude.
%   h: altitude.
%   el: elevation angle.
%   type: modeling option.
% outputs:
%   Ps: pressure in HPa.
%    T: temperature in kelvin.
%Author: xiahaa@space.dtu.dk
    % compute p and t
    [Ps, T] = calc_ps_T_via_earth_atmosphere_model(h,type);
    % compute e
    e = calc_partial_pressure(T);
    
    % the Saastamoinen model
    if type == '1'
        % type 1 from slides
        D = 1+0.0026*cos(deg2rad(lati)*2)+0.00028*h/1000;
        d_trop = 0.002277*D*(Ps+(1255/T+0.05)*e);
        % consider elevation angle
        d_trop = d_trop * 1/sin(deg2rad(el));
    elseif type == '2'
        % from Global Positioning System.
        d_dry = 0.002277*(1+0.0026*cos(deg2rad(lati)*2)+0.00028*h/1000)*Ps;
        d_wet = 0.002277*(1255/T+0.05)*e;
        d_trop = d_dry * 1/(sin(deg2rad(el))+(0.00143/(tan(deg2rad(el))+0.0445))) ...
            + d_wet * 1/(sin(deg2rad(el))+(0.00035/(tan(deg2rad(el))+0.017)));
    elseif type == '3'
        % from ESA Navipedia: Gelileo Troposheric Correction Model.
        d_trop = 0.002277*Ps/(1-0.0026*cos(deg2rad(lati)*2)-0.00028*h/1000) ...
            + 0.002277 * (1255/T+0.05)*e;
        d_trop = d_trop * 1/sin(deg2rad(el));
    end
end

function [Ps, T] = calc_ps_T_via_earth_atmosphere_model(h,type)
%calc_ps_T_via_earth_atmosphere_model calculate pressure using  earth 
%atmosphere model given by NASA.
%model and constants from:
%   https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
% inputs:      
%   h: altitude
% outputs:
%   Ps: pressure in HPa.
%    T: temperature in kelvin.
%Author: xiahaa@space.dtu.dk
    if type == '2' || type == '3'
        %% use nasa parameters
        T = 15.04 - 0.00649*h;% temperature in Celsius degrees
        T = T + 273.1;
        Ps = 101.29*((T)/288.08)^(5.256);% pressure in KPa
        Ps = Ps * 10;% KPa to HPa
    elseif type == '1'
        %% use parameters given in the slides.
        T = 18;
        T = T + 273.1;
        Ps = 1013;%hpa to kpa
    end
end

function e = calc_partial_pressure(T)
%calc_partial_pressure calculate partial pressure due to water vapor.
%model from: Global Positioning System.
% inputs:
%   T: temperature in kelvin.
% ouputs:
%   e: partial pressure.
%Author: xiahaa@space.dtu.dk
    RH = 0.5;%assume relative humidity is 0.5
    e = 6.108*RH*exp((17.15*T - 4684)/(T-38.45));%mbars=HPa
end
