function Irec = undistortImg(I, K, k1, k2, k3, p1, p2)
% function Irec = undistortImg(I, K, k1, k2, k3, p1, p2)
%
% Rectify image.
%   Inputs:
%       I: input image, gray or color.
%       K: camera intrinsics.
%       k1,k2,k3: radial distortion parameters.
%       p1,p2: tangent distortion parameters.
%   Outputs:
%       Irec: Rectified image.
%
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.

    [xx0,yy0] = meshgrid(1:size(I,2),1:size(I,1));
    xx0 = xx0(:)';
    yy0 = yy0(:)';
    % normalized
    xx1 = (xx0 - K(1,3))./K(1,1);
    yy1 = (yy0 - K(2,3))./K(2,2);
    % projection and distortion
    [xx,yy] = forward(xx1,yy1,k1,k2,k3,p1,p2);
    % interpolation
    xxf = xx.*K(1,1) + K(1,3);
    yyf = yy.*K(2,2) + K(2,3);
    xx = floor(xxf);
    yy = floor(yyf);
    valid = xx > 1 & xx <= size(I,2)-1 & yy > 1 & yy <= size(I,1)-1;
    xx = xx(valid);xxf = xxf(valid);
    yy = yy(valid);yyf = yyf(valid);

    m = size(I,1);n = size(I,2);
    ind0 = sub2ind([m,n], yy, xx);
    ind1 = sub2ind([m,n], yy, xx+1);
    ind2 = sub2ind([m,n], yy+1,xx);
    ind3 = sub2ind([m,n], yy+1,xx+1);

    % bilinear interpolation
    dx = xxf - xx;
    dy = yyf - yy;
    d0 = (1-dx).*(1-dy);
    d1 = dx.*(1-dy);
    d2 = (1-dx).*dy;
    d3 = dx.*dy;

    ind4 = sub2ind([m,n], round(yy0(valid)), round(xx0(valid)));

    if size(I,3) == 1
        Irec = uint8(zeros(m,n));
        Irec(ind4) = uint8(double(I(ind0)).*d0+double(I(ind1)).*d1+double(I(ind2)).*d2+double(I(ind3)).*d3);
    else
        Ir = I(:,:,1);Ig = I(:,:,2);Ib = I(:,:,3);
        Irecr = uint8(zeros(m,n));
        Irecg = uint8(zeros(m,n));
        Irecb = uint8(zeros(m,n));
        Irecr(ind4) = uint8(double(Ir(ind0)).*d0+double(Ir(ind1)).*d1+double(Ir(ind2)).*d2+double(Ir(ind3)).*d3);
        Irecg(ind4) = uint8(double(Ig(ind0)).*d0+double(Ig(ind1)).*d1+double(Ig(ind2)).*d2+double(Ig(ind3)).*d3);
        Irecb(ind4) = uint8(double(Ib(ind0)).*d0+double(Ib(ind1)).*d1+double(Ib(ind2)).*d2+double(Ib(ind3)).*d3);
        Irec = cat(3,Irecr,Irecg,Irecb);
    end
end


