function It = gradIntensity(x,y,flowu,flowv,im1,im2)
% perform intensity gradient for LK optical flow.
    shiftx = x + flowu;
    shifty = y + flowv;
    
    [xx,yy] = meshgrid(1:size(im1,2),1:size(im1,1));
    
    im_warp = interp2(xx,yy,im2,shiftx,shifty,'linear',NaN);
    mask = isnan(im_warp);
    im_warp(mask) = im1(mask);

    It = im_warp - im1((x-1)*size(im1,1)+y);
end