function ptres = undistortPoint(ptsrc, K, k1, k2, k3, p1, p2)
% function ptres = undistortPoint(ptsrc, K, k1, k2, k3, p1, p2)
%
% Undistorb points.
%
%   Inputs:
%       ptsrc: raw distorted 2D points in image - 2xN.
%       k1, k2, k3, p1, p2: radial and tangent parameters.
%		K: camera intrinsics.
%   Outputs:
%       ptres: undistorted 2D points.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    % normalized
    x0 = (ptsrc(1,:) - K(1,3))./K(1,1);
    y0 = (ptsrc(2,:) - K(2,3))./K(2,2);
    
    ptres = ptsrc;
    
    for i = 1:size(x0, 2)
        iter = 1;
        maxiter = 20;
        % iteration
        x = x0(i);y = y0(i);
        while iter <= maxiter
            r2 = x^2 + y^2;
            r4 = r2*r2;
            r6 = r4*r2;
            xy = x*y;
            
            % backward
            newxx = (x0(i) - (2*p1*xy + p2*(r2 + 2*x^2)))/(1 + r2*k1+r4*k2+r6*k3);
            newyy = (y0(i) - (p1*(r2 + 2*y^2) + 2*p2*xy))/(1 + r2*k1+r4*k2+r6*k3);
    
            % forward
            [x1,y1] = forward(newxx,newyy,k1,k2,k3,p1,p2);
            if abs(x1-x0(i)) < 1e-3 && abs(y1-y0(i)) < 1e-3
                break;
            end         
            
            x = newxx;
            y = newyy;
            
            iter = iter + 1;
        end
        ptres(:,i) = [x;y];
    end
    ptres(1,:) = ptres(1,:).*K(1,1) + K(1,3);
    ptres(2,:) = ptres(2,:).*K(2,2) + K(2,3);
end


     




