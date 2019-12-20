function [x,y,z,latr,lonr,heightr] = coordinate_transformation(lat,lon,height)
%   [x,y,z,latr,lonr,heightr] = coordinate_transformation(lat,lon,height)
%   does the forward and backward conversion from Cartesian coordiantes 
%   to latitude, longitude, altitude.
%   lat,lon,height: coordinates in the ellipsoidal coordinate system.
%   x,y,z: coordinates in the ECEF coordinate system.
%   latr,lonr,heightr: reconstructed coordinates in the 
%   ellipsoidal coordinate system.
%   Author: xiahaa@space.dtu.dk

    consParams = struct('a',6378137.0,'f',1/298.257223563);
    addpath('../utils/')
    [x,y,z] = llhtoCartesian(lat, lon, height, consParams);
     %% reverse
    [latr, lonr, heightr] = Cartesian2llh(x,y,z,consParams);
end

