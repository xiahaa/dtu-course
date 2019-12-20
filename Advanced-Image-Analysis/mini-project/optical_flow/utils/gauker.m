function ker = gauker(hsize,sigma)
    % gaussian filter
    ker = gaussian_kernel_calculator(2, hsize, sigma);
end