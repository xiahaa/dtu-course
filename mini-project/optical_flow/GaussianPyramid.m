function pyr = GaussianPyramid(im, layers, verbose)
    pyr = cell(layers,1);
    pyr{1} = im;
    sigma = 1;
    for i = 2:layers
        % gaussian filter
        img = imgaussfilt(pyr{i-1}, sigma);
        % subsampling
        pyr{i} = imresize(img,0.5,'nearest');
    end
    
    if verbose == true
        % debug
        w = size(im,2);
        h = 0;
        for i = 1:layers
            h = h + size(pyr{i},1);
        end
        imgyr1 = zeros(h,w);
        h = 1;
        for i = 1:layers
            imgyr1(h:h+size(pyr{i},1)-1,1:size(pyr{i},2)) = pyr{i};
            h = h + size(pyr{i},1);
        end
        figure;imshow(imgyr1);
    end
end