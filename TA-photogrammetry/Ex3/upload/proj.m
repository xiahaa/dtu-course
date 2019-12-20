
function [uv1, in] = proj(R, t, p, K)
% function [uv1, in] = proj(R, t, p, K)
%
% Projection of 3D points.
%   Inputs:
%       R,t: extrinsics.
%       K: camera intrinsics.
%       p: 3D points.
%   Outputs:
%       uv1: projected image pixels.
%       in??? array of indicators. Each indicator indicate whether this point is inside the image: 1 inside, 0 otherwise.
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    P1 = K*([R t]);
    phomo = tohomo(p);
    uv1 = P1*phomo;
    uv1 = uv1./uv1(3,:);
    in = uv1(1,:) > 0 & uv1(1,:) < K(1,3)*2 & uv1(2,:) > 0 & uv1(2,:) < K(2,3)*2;
end
