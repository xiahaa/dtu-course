function varargout = create_scale_normalized_LoG(Igray, t0, K)
    tsqrt = t0;% prev scale
    LLs = zeros(size(Igray,1),size(Igray,2),K);
    radius = zeros(1,K);
    for k = 1:K
        [LL] = scale_normalized_LoG(Igray, tsqrt);
        LLs(:,:,k) = LL;
        tsqrt = tsqrt * sqrt(2); % scale update
        radius(k) = tsqrt;
    end
    varargout{1} = LLs;
    varargout{2} = radius;
end

function [LL] = scale_normalized_LoG(I, tsqrt)
    sigma = 3;
    g = gassian_fast(tsqrt, sigma);
    dgg = secondorder_der_gaussian_fast(tsqrt, sigma);
    % gaussian filter
    Ixx = imfilter(I, dgg, 'replicate');% along x
    Ixx = imfilter(Ixx, g', 'replicate');% smooth y, since g=gx*gy, dg=dgx gy + gx dgy. partial dev only pick one item
    Iyy = imfilter(I, dgg', 'replicate');% along y
    Iyy = imfilter(Iyy, g, 'replicate');%along x
    % LoG
    LL = Ixx+Iyy;
    % normalized LoG
    LL = LL .* (tsqrt^2);
end
