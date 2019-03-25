function imseg = segmentation(im, configuration, f, miuf)
    imseg = im;
    for i = 1:length(f)
        imseg(configuration == f(i)) = miuf(f(i));
    end
end