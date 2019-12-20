
function params = line_estimation(ph)
% function params = line_estimation(ph)
%
% Generate simulated data for 2D line fitting.
%   Inputs:
%       ph: the homogeneous coordinates of points.
%   Outputs:
%   params: line parameters.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    if size(ph,2) == 2
        % minimum set, use cross product
        params = cross(ph(:,1),ph(:,2));
    else
        A = ph';
        [~,~,V] = svd(A);
        params = V(:,end);
    end
end