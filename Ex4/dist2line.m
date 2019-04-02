
function d = dist2line(ph, params)
% function d = dist2line(ph, params)
%
% Compute point to line distance.
%   Inputs:
%       ph: the homogeneous coordinates of points.
%   params: line parameters.
%   Outputs:
%       d: point to line distance.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    params = params ./ norm(params(1:2));
    d = abs(ph' * params);
end