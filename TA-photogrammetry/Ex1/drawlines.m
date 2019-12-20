function I = drawlines(I,q,indices)
% function I = drawlines(I,q,order)
%
% Draw lines on a image I.
%   Inputs:
%       I: image - MxNx3.
%       q: Point set, contains point coordinates - 2xN. 
%       indices, index to choose the start and end point - Nx2.
%   Outputs:
%       I: image with pixels on lines colored - MxNx3.
% This function has been vectorized for boosting the performance.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    ccmap = jet(size(indices,1));
    for i = 1:size(indices,1)
        I = drawline(I,[q(1,indices(i,1)) q(1,indices(i,2))], [q(2,indices(i,1)) q(2,indices(i,2))], ccmap(i,:));
    end
end

function I = drawline(I,x,y,color)
    M = size(I,1);N = size(I,2);
    nPoints = max(abs(diff(x)), abs(diff(y)))+1;  % Number of points in line
    rIndex = round(linspace(y(1), y(2), nPoints));  % Row indices
    cIndex = round(linspace(x(1), x(2), nPoints));  % Column indices
    index = sub2ind([M,N], rIndex, cIndex);     % Linear indices
    nn = [-M +M -1 +1];
    index = [index index+nn(1) index+nn(2) index+nn(3) index+nn(4)];
    index = index(index > 0 & index < M*N);
    I(index) = color(1);  % Set the line pixels to the max value of 255 for uint8 types
    I(index+M*N) = color(2);
    I(index+M*N*2) = color(3);
end