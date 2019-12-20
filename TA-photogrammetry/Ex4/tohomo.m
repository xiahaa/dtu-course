function xh = tohomo(x)
% function xh = tohomo(x)
%
%From inhomogeneous coordinate to homogeneous coordiante.
%   Inputs:
%       x: inhomogeneous coordinate.
%   Outputs:
%       xh: homogeneous coordiante.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.

    xh = [x;ones(1,size(x,2))];
end