function [xx,yy] = forward(xx,yy,k1,k2,k3,p1,p2)
% function [xx,yy] = forward(xx,yy,k1,k2,k3,p1,p2)
%
% Apply distorsiton forward.
%   Inputs:
%       xx: x coordinates.
%       yy: y coordinates.
% k1,k2,k3: radial distorsion coefficients.
%    p1,p2: tangent distorsion coefficients.
%   Outputs:
%    xx,yy: distorted coordiantes..
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    r2 = xx.^2 + yy.^2;
    r4 = r2.*r2;
    r6 = r4.*r2;
    xy = xx.*yy;
    
    dradial = 1+r2.*k1+r4.*k2+r6.*k3;
    dtangentx = 2*p1.*xy + p2.*(r2 + 2.*xx.^2);
    dtangenty = p1.*(r2 + 2.*yy.^2) + 2*p2.*xy;

    xx = xx.*dradial + dtangentx;
    yy = yy.*dradial + dtangenty;
end