function [x,y,z] = llhtoCartesian(lat, lon, height, consParams)
%   [x,y,z] = llhtoCartesian(lat, lon, height, consParams) does the 
%   conversion from latitude, longitude, altitude to the Cartesian 
%   coordiantes in the ECEF coordiante system.
%   lat, lon, height: coordinates in the ellipsoidal coordinate system.
%   consParams: struct of relative constant parameters
%       a: length of semi-major axis in meters.
%       f: flattening of the ellipsoid.
%   Author: xiahaa@space.dtu.dk

    a = consParams.a;
    f = consParams.f;
    % fistly, convert f to e
    e = sqrt(2*f-f*f);
    rlat = deg2rad(lat);
    rlon = deg2rad(lon);
    e2 = e*e;
    % then compute N
    N = a / sqrt(1-e2*sin(rlat)*sin(rlat));
    % XYZ, zc = height*sin(rlat);
    hsum = N+height;
    x = hsum*cos(rlat)*cos(rlon);
    y = hsum*cos(rlat)*sin(rlon);
    z = (N*(1-e2)+height)*sin(rlat);
end