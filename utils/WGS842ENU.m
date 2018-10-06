function [e,n,u] = WGS842ENU(lati, longi, dx, dy, dz)
%WGS842ENU: convert coordinates from WGS84 to ENU.
%   [e,n,u] = WGS842ENU(lati, longi, dx, dy, dz) does the conversion
%   from WGS84 coordinates to ENU coordinates.
%   Author: xiahaa@space.dtu.dk
    R1 = consR(deg2rad(90-lati),1);
    R3 = consR(deg2rad(90+longi),3);
    R = R1*R3;
    v1 = R*[dx;dy;dz];
    e = v1(1);n = v1(2);u = v1(3);
end