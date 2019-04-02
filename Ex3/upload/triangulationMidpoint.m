
function P = triangulationMidpoint(x1,P1,x2,P2)
% function P = triangulationMidpoint(x1,P1,x2,P2)
%
%Implementation of mid-point triangulation method.
%   Inputs:
%       x1: point in image 1, homogeneous coordinates.
%       P1: camera projection matrix of image 1.
%       x2: point in image 2, homogeneous coordinates.
%       P2: camera projection matrix of image 2.
%   Outputs:
%       P: Triangulated 3D point.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    M1 = P1(1:3,1:3);
    c1 = -M1\P1(1:3,4);

    M2 = P2(1:3,1:3);
    c2 = -M2\P2(1:3,4);

    % ray
    r1 = M1\x1;
    r2 = M2\x2;

    A = [r1 -r2];
    b = c2 - c1;

    x = A\b;

    P = (c1 + x(1).*r1 + c2 + x(2).*r2)*0.5;
end
