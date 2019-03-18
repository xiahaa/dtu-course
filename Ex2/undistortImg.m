function Irec = undistortImg(I, K, k1, k2, k3, p1, p2)
    [xx0,yy0] = meshgrid(1:size(I,2),1:size(I,1));
    xx0 = xx0(:)' - 1;
    yy0 = yy0(:)' - 1;
    
    xx1 = (xx0 - K(1,3))./K(1,1);
    yy1 = (yy0 - K(2,3))./K(2,2);
    
%     iter = 1;
%     maxiter = 1;
     
%     r2 = xx.^2 + yy.^2;
%     r4 = r2.*r2;
%     r6 = r4.*r2;
%     xy = xx.*yy;
    
    [xx,yy] = forward(xx1,yy1,k1,k2,k3,p1,p2);
    
    xx = round(xx.*K(1,1) + K(1,3));
    yy = round(yy.*K(2,2) + K(2,3));
    valid = xx > 0 & xx <= size(I,2) & yy > 0 & yy <= size(I,1);
    m = size(I,1);n = size(I,2);
    ind0 = sub2ind([m,n], round(yy(valid)), round(xx(valid)));
    ind1 = sub2ind([m,n], round(yy0(valid)+1), round(xx0(valid)+1));
    Irec = uint8(zeros(m,n));
    Irec(ind1) = I(ind0);
    
    
    % backward
%     newxx = (xx - (2*p1.*xy + p2.*(r2 + 2.*xx.^2)))./(1 + r2.*k1+r4.*k2+r6.*k3);
%     newyy = (yy - (p1.*(r2 + 2.*yy.^2) + 2*p2.*xy))./(1 + r2.*k1+r4.*k2+r6.*k3);
%     
%     while iter <= maxiter
%         % forward
%         [xx1,yy1] = forward(newxx,newyy,k1,k2,k3,p1,p2);
%         if max(abs(xx1-xx)) < 1e-3 && max(abs(yy1-yy)) < 1e-3
%             break;
%         end
%         
%         r2 = newxx.^2 + newyy.^2;
%         r4 = r2.*r2;
%         r6 = r4.*r2;
%         xy = newxx.*newyy;
%         
%         newxx = (xx - (2*p1.*xy + p2.*(r2 + 2.*newxx.^2)))./(1 + r2.*k1+r4.*k2+r6.*k3);
%         newyy = (yy - (p1.*(r2 + 2.*newyy.^2) + 2*p2.*xy))./(1 + r2.*k1+r4.*k2+r6.*k3);
% 
%         iter = iter + 1;
%     end
    
%     xx = round(newxx.*K(1,1) + K(1,3));
%     yy = round(newyy.*K(2,2) + K(2,3));
%     
%     valid = xx > 0 & xx <= size(I,2) & yy > 0 & yy <= size(I,1);
%     m = size(I,1);n = size(I,2);
%     ind0 = sub2ind([m,n], round(yy(valid)), round(xx(valid)));
%     Irec = uint8(zeros(m,n));
%     Irec(ind0) = I(ind0);
end

function [xx,yy] = forward(xx,yy,k1,k2,k3,p1,p2)
    r2 = xx.^2 + yy.^2;
    r4 = r2.*r2;
    r6 = r4.*r2;
    xy = xx.*yy;
    
    dradial = 1+r2.*k1+r4.*k2+r6.*k3;
    dtangentx = 2*p1.*xy + p2.*(r2 + 2.*xx.^2);
    dtangenty = p1.*(r2 + 2.*yy.^2) + 2*p2.*xy;

    xx = xx.*dradial + dtangentx;
    yy = yy.*dradial + dtangenty;
end
