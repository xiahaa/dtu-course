
function xyz = stereoTriangulate(disp, f, cx, cy, b,id)
% function xyz = stereoTriangulate(disp, f, cx, cy, b,id)
%
% Triangulation for stereo vision using disparity map.
%   z = f*b/d;
%   x = （x-cx)*z/f
%   y = （y-cy)*z/f
%
%   Inputs:
%       disp: disparity map.
%       f: focus length.
%       cx,cy: principal points.
%       b: baseline.
%       id: valid pixels in the disparity map.
%   Outputs:
%       xyz: Triangulated 3D points.
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    z = f * b ./ disp(id);
    [yy,xx] = find(id);
    x = (xx - cx).*z / f;
    y = (yy - cy).*z / f;
    xyz = [x y z];
end
