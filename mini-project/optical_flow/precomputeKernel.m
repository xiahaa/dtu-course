function kernel = precomputeKernel(im, hsize)
% this function precompute the kernel for later optical flow computation.
    width = size(im,2);
    height = size(im,1);
    
    % lagest buffer, will be cut
    m2 = (hsize*2+1)*(hsize*2+1);
    is = zeros(width*height*m2,1);
    js = zeros(width*height*m2,1);
    ss = ones(width*height*m2,1);
    k = 1;
    
    % local kernel id
    [lx,ly] = meshgrid(-hsize:1:hsize,-hsize:1:hsize);
    lx = vec(lx'); ly = vec(ly');
    for j = 1:width
        for i = 1:height
            % local coordinate
            llx = i + lx;
            lly = j + ly;
            % within image region
            valid = llx > 0 & llx <= height & lly > 0 & lly <= width;
            num_valid = sum(valid);
            % to index
            index = (lly-1).*height + llx;
            % assignment
            is(k:k+num_valid) = (j-1)*height + i;
            js(k:k+num_valid) = index(valid);
            % increament k
            k = k + num_valid + 1;
        end
    end
    is(k:end) = []; js(k:end) = []; ss(k:end) = [];
    % sparse kernel from image to LK response
    kernel = sparse(is,js,ss);
end