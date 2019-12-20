function [azimuth, zenith, elevation] = calcAzimuthZenithElevation(e,n,u)
%   [azimuth, zenith, elevation] = calcAzimuthZenithElevation(e,n,u) does
%   the computation of the azimuth, zenith, elevation angle for a given
%   point.
%   Author: xiahaa@space.dtu.dk
    azimuth = rad2deg(atan2(e,n));
    zenith = 90 - rad2deg(asin(u / sqrt(e^2+n^2+u^2)));
    elevation = rad2deg(asin(u / sqrt(e^2+n^2+u^2)));
end