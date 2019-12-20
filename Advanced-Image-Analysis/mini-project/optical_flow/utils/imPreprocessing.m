function im = imPreprocessing(im)
    if size(im,3) == 3
        im = rgb2gray(im);
    end
    im = im2double(im);
%     ker = fspecial('gaussian',3,1);
%     im = imfilter(im,ker,'replicate','same');
end