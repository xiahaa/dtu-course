function [Ix, Iy] = grad3(im1)
% third way of computing the gradient.
    h = [-1 9 -45 0 45 -9 1]/60;
    Ix = imfilter(im1,h,'symmetric','same');
    Iy = imfilter(im1,h','symmetric','same');
end