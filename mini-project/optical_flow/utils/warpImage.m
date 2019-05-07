function im_warp = warpImage(im1, u, v, im2)
% perform warping, fill in nan for missing pixels
    [x,y] = meshgrid(1:size(im1,2),1:size(im1,1));
    shiftx = x + u;
    shifty = y + v;
    if size(im2,3) == 1
        im_warp = warpSingleCh(im1,im2,x,y,shiftx,shifty);
    else
        im_warp = zeros(size(im1,1),size(im1,2),size(im1,3));
        for i = 1:size(im2,3)
            im_warp(:,:,i) = warpSingleCh(im1(:,:,i),im2(:,:,i),x,y,shiftx,shifty);
        end
    end
end

function im_warp = warpSingleCh(im1,im2,x,y,shiftx,shifty)
% perform single channel warping, fill in nan for missing pixels
    im_warp = interp2(x,y,im2,shiftx,shifty,'linear',NaN);
    mask = isnan(im_warp);
    im_warp(mask) = im1(mask);
end