function p = fromhomogeneous(ph)
% function p = fromhomogeneous(ph)
%
% From homogeneous coordinate to inhomogeneous coordinate.
%   Inputs:
%       ph: homogeneous coordinates.
%   Outputs:
%       ph: inhomogeneous coordiantes.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    p = ph(1:end-1,:) ./ ph(end,:);
end