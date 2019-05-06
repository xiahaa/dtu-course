function im_warp = warpImage(im1, u, v, im2)
    [x,y] = meshgrid(1:size(im1,2),1:size(im1,1));
    shiftx = x + u;
    shifty = y + v;
    % perform warping, fill in nan for missing pixels
    im_warp = interp2(x,y,im2,shiftx,shifty,'linear',NaN);
    mask = isnan(im_warp);
    im_warp(mask) = im1(mask);
end