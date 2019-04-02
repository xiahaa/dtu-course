function F = DLT8pt(x1,x2)
% function F = DLT8pt(x1,x2)
%
%Implementation of 8 point algorithm for fundamental matrix estimation.
%
%   Inputs:
%       x1: points in image 1, at least 8.
%       x2: corresponding points in image 2.
%   Outputs:
%       F: fundamental matrix.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    if size(x1,1) ~= 3
        x1h = tohomo(x1);
    end
    if size(x2,1) ~= 3
        x2h = tohomo(x2);
    end

    A = zeros(size(x1,2),9);
    for i = 1:size(x1,2)
        A(i,:) = kron(x1h(:,i)', x2h(:,i)');
    end

    [~,~,V] = svd(A);
    F = V(:,end);
    F = F([1 4 7;2 5 8;3 6 9]);
    F = F./F(3,3);
end