function [lat, lon, height] = Cartesian2llh(x,y,z,consParams)
%   [lat, lon, height] = Cartesian2llh(x,y,z,consParams) does the 
%   conversion from Cartesian coordiantes to latitude, longitude,
%   altitude.
%   x, y, z: coordinates in the ECEF coordinate system.
%   consParams: struct of relative constant parameters
%       a: length of semi-major axis in meters.
%       f: flattening of the ellipsoid.
%   Author: xiahaa@space.dtu.dk

    a = consParams.a;
    f = consParams.f;
    % fistly, convert f to e
    e = sqrt(2*f-f*f);
    e2 = e*e;
    x2y2sqr = sqrt(x*x+y*y);
    
    iter = 1;
    vold = [0,0,0];
    vnew = zeros(1,3);
    tolerance = 1e-6;
    maxIter = 1e6;
    
    while iter < maxIter
%         disp(strcat('iter:',num2str(iter)));
        N = a / sqrt(1-e2*sin(vold(1))*sin(vold(1)));
        % update
        vnew(1) = atan2(z, x2y2sqr*(1-e2*N/(N+vold(3)+1e-6)));
        vnew(2) = atan2(y,x);
        vnew(3) = x2y2sqr/cos(vnew(1)) - N;
        % check tolerance
        verr = vnew - vold;
        if norm(verr) < tolerance
            break;
        else
            vold = vnew;
            iter = iter + 1;
        end
    end
    % to deg
    lat = rad2deg(vnew(1));
    lon = rad2deg(vnew(2));
    height = vnew(3);
end