function [Ix, Iy] = grad1(im1)
% perform the simplest gradient calculation for LK optical flow.
    hx = [-1,0,1];
    Ix = imfilter(im1,hx,'symmetric','same');
    Iy = imfilter(im1,hx','symmetric','same');
end