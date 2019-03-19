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

