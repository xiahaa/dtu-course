function R = consR(rad, axis)
%consR: construct a rotation matrix.
%   R = consR(rad, axis) return the corresponding rotation matrix.
%   rad: rotation angle.
%   axis: 1-x axis, 2-y axis, 3-z axis.
%   R: rotation matrix.
%   Author:xiahaa@space.dtu.dk
    if axis == 1
        R = [1 0 0; ...
             0 cos(rad) sin(rad); ...
             0 -sin(rad) cos(rad)];
    elseif axis == 3
        R = [cos(rad) sin(rad) 0; ...
            -sin(rad) cos(rad) 0; ...
            0 0 1];
    elseif axis == 2
        R = [cos(rad) 0 -sin(rad); ...
             0 1 0; ...
             sin(rad) 0 cos(rad)];
    end
end