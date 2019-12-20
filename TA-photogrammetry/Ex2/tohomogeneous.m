function ph = tohomogeneous(p)
% function ph = tohomogeneous(p)
%
% From inhomogeneous coordinate to homogeneous coordinate.
%   Inputs:
%       p: inhomogeneous coordinates.
%   Outputs:
%       ph: homogeneous coordiantes.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    ph = [p;ones(1,size(p,2))];
end
