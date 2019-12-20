function gast = jd2gast(jd)
%% julian date to approximate Greenwich Apparent Sidereal Time
% reference: https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html
% Author: xiahaa@space.dtu.dk
    jd0 = floor(jd) + 0.5;
    d0 = jd - jd0;
    H0 = d0 .* 24;
    d1 = jd - 2451545.0;%Julian centuries from 2000 Jan. 1 12h UT1
    T = d1 ./ 36525;
    D0 = jd0 - 2451545.0;
    gmst = 6.697374558 ...
         + 0.06570982441908 * D0 ...
         + 1.00273790935 * H0 ...
         + 0.000026 * T^2
    %% Greenwich Apparent Sidereal Time (GAST) is Greenwich Mean Sidereal 
    %   Time (GMST) corrected for the shift in the position of the vernal 
    %   equinox due to nutation. 
    % use the approximation
    gast = mod(gmst,24).*15;
end