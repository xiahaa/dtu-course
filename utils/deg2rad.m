function rad = deg2rad(deg)
%deg2rad    conversion from degree to radian.
%   rad = deg2rad(deg) produces a corresponding radian value to a 
%   given degree value.
%   Author: xiahaa@space.dtu.dk
    rad = deg .* pi ./ 180.0;
end