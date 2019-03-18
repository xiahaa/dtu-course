function rot = Rxyz(roll, pitch, yaw)
% function rot = Rxyz(roll, pitch, yaw)
%
% Create a rotation matrix corresponding to the given Euler angles. The
% rotation order is assumed to be fistly along z-axis, then y-axis, and
% finally x-axis.
%   Inputs:
%       roll: rotation angle along x-axis, rad.
%       pitch: rotation angle along y-axis, rad.
%       yaw: rotation angle along z-axis, rad.
%   Outputs:
%       rot: Rotation Matrix
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    rot1 = [cos(yaw) sin(yaw) 0;-sin(yaw) cos(yaw) 0; 0 0 1];
    rot2 = [cos(pitch) 0 -sin(pitch);0 1 0; sin(pitch) 0 cos(pitch)];
    rot3 = [1 0 0;0 cos(roll) sin(roll);0 -sin(roll) cos(roll)];

    rot = rot3*rot2*rot1;
    rot = rot';
end