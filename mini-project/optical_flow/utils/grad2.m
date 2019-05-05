function [Ix, Iy] = grad2(im1)
% perform gradient operation with first order gaussian kernel
    hsize = 2;% or 1,3
    sigma = 1;% or 1
    x = -hsize:1:hsize;
    cons1 = sigma*sigma;
    hg = 1/sqrt(2*pi*cons1).*exp(-x.^2./(2*cons1)); hg = hg ./ sum(abs(hg(:)));
    hgx = 1/sqrt(2*pi*cons1).*exp(-x.^2./(2*cons1)).*(-x./cons1); hgx = hgx ./ sum(abs(hgx(:)));
    hgx = fliplr(hgx);% convolution kernel is the flip version.
    Ix = imfilter(im1,hgx,'replicate','same');
%     Ix = imfilter(Ix,hg','replicate','same');
    Iy = imfilter(im1,hgx','replicate','same');
%     Iy = imfilter(Iy,hg,'replicate','same');
end