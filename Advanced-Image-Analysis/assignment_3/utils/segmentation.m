function imseg = segmentation(im, configuration, f, miuf)
% segmentation an image
    imseg = im;
    for i = 1:length(f)
        imseg(configuration == f(i)) = miuf(f(i));
    end
end