function [R,t] = lssol(p,q)
% function [R,t] = lssol(p,q)
%
% Implementation of the simple linear solution for absolute orientation
% proposed in
% Haralick B M, Lee C N, Ottenberg K, et al. 
% Review and analysis of solutions of the three point perspective pose estimation problem[J]. 
% International journal of computer vision, 1994, 13(3): 331-356.
%
%   Inputs:
%       p: 3D points in World coordinate system - 3xN.
%       q: 3D points in Camera coordinate system, - 3xN. 
%   Outputs:
%       R: rotation matrix.
%       t: translation vector.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.

    % transformation
    v1 = p(:,2) - p(:,1);v1 = v1./norm(v1);
    v2 = p(:,3) - p(:,1);v2 = v2./norm(v2);
    v3 = cross(v1,v2);v3 = v3./norm(v3);
    v4 = [0;0;1];
    v5 = cross(v3,v4);v5 = v5./norm(v5);
    theta = acos(dot(v4,v3));
    K = [0 -v5(3) v5(2);v5(3) 0 -v5(1);-v5(2) v5(1) 0];
    R1 = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
    
    pnew = R1*p;
    d = mean(pnew(3,:));
    A = [pnew(1,1) pnew(2,1) 0 0 0 0 1 0 0; ...
         0 0 pnew(1,1) pnew(2,1) 0 0 0 1 0; ...
         0 0 0 0 pnew(1,1) pnew(2,1) 0 0 1; ...
         pnew(1,2) pnew(2,2) 0 0 0 0 1 0 0; ...
         0 0 pnew(1,2) pnew(2,2) 0 0 0 1 0; ...
         0 0 0 0 pnew(1,2) pnew(2,2) 0 0 1; ...
         pnew(1,3) pnew(2,3) 0 0 0 0 1 0 0; ...
         0 0 pnew(1,3) pnew(2,3) 0 0 0 1 0; ...
         0 0 0 0 pnew(1,3) pnew(2,3) 0 0 1];
    B = q(:);
    X = A\B;
    r11 = X(1);r12 = X(2);r21 = X(3);r22 = X(4);r31 = X(5);r32 = X(6);
    r13 = r21*r32-r22*r31;
    r23 = r12*r31-r11*r32;
    r33 = r11*r22-r12*r21;
    R2 = [r11 r12 r13;r21 r22 r23;r31 r32 r33];
    R = R2*R1;
    t = -R2*[0;0;d]+[X(7);X(8);X(9)];
end