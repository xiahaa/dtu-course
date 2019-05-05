clc;close all;clear all;

addpath ../../utils;

baseDir = '../../data/optical_flow_data/Basketball';

% load database
buildingScene = imageDatastore(baseDir);
numImages = numel(buildingScene.Files);

% load a certain image
im1 = imread(buildingScene.Files{1});
im1 = imPreprocessing(im1);

im2 = imread(buildingScene.Files{2});
im2 = imPreprocessing(im2);

% image size
[flow_u, flow_v] = denseflowHS(im1, im2);


% display dense flow as an image
img = computeColor(flow_u,flow_v);
figure
if size(img,1)*size(img,2) < 1000
    imshow(img, 'InitialMagnification',1000);hold on;
else
    imshow(img);
end

function showFlowQuiver(im, flow_u, flow_v)
    figure;
    if size(im,1)*size(im,2) < 1000
        imshow(im, 'InitialMagnification',1000);hold on;
    else
        imshow(im);hold on;
    end
    height = size(im,1);
    width = size(im,2);
    % opflow = opticalFlow(flow_u,flow_v);
    % plot(opflow,'DecimationFactor',[1 1],'ScaleFactor',1);
    [yy,xx] = meshgrid(1:height,1:width);
    xx = vec(xx');
    yy = vec(yy');
    quiver(xx,yy,flow_u(:),flow_v(:),'LineWidth',1.5, 'Color','r','MaxHeadSize',1);axis image
end

function [Ix, Iy] = grad2(im1)
% perform gradient operation with first order gaussian kernel
    hsize = 2;% or 1,3
    sigma = 0.5;% or 1
    x = -hsize:1:hsize;
    cons1 = sigma*sigma;
    hg = 1/sqrt(2*pi*cons1).*exp(-x.^2./(2*cons1)); hg = hg ./ sum(abs(hg(:)));
    hgx = 1/sqrt(2*pi*cons1).*exp(-x.^2./(2*cons1)).*(-x./cons1); hgx = hgx ./ sum(abs(hgx(:)));
    hgx = fliplr(hgx);% convolution kernel is the flip version.
    Ix = imfilter(im1,hgx,'replicate','same');
    Ix = imfilter(Ix,hg','replicate','same');
    Iy = imfilter(im1,hgx','replicate','same');
    Iy = imfilter(Iy,hg,'replicate','same');
end

function [Ix, Iy] = grad3(im1)
% third way of computing the gradient.
    h = [-1 9 -45 0 45 -9 1]/60;
    Ix = imfilter(im1,h,'symmetric','same');
    Iy = imfilter(im1,h','symmetric','same');
end

function It = gradIntensity(x,y,flowu,flowv,im1,im2)
% perform intensity gradient for LK optical flow.
    shiftx = x + flowu;
    shifty = y + flowv;

    im_warp = interp2(x,y,im2,shiftx,shifty,'linear',NaN);
    mask = isnan(im_warp);
    im_warp(mask) = im1(mask);

    It = im_warp - im1;
end

function [flow_u, flow_v] = denseflowPyrLK(im1,im2)
    layers = 4;
    pyr1 = GaussianPyramid(im1, layers, 1);
    pyr2 = GaussianPyramid(im2, layers, 1);
   
    hsize = 1;
    
    % last layer
    [pu, pv] = denseflowLK(pyr1{end}, pyr2{end}, [], [], hsize);
    for i = layers-1:-1:1
        pu = pu.*2;
        pv = pv.*2;
        u = imresize(pu,size(pyr1{i}),'bilinear');
        v = imresize(pv,size(pyr2{i}),'bilinear');
        [pu, pv] = denseflowLK(pyr1{i}, pyr2{i}, u, v, hsize);
    end
    flow_u = pu;
    flow_v = pv;
end

function [flow_u, flow_v] = resampleFlow(u, v, sz)
% instead of using hard scale, here I do the flow updating using bilinear
% interpolation and actual ratio between cascaded layers.
    ratio = sz(1) / size(u,1);
    flow_u = imresize(u,sz,'binlinear').*ratio;
    flow_v = imresize(v,sz,'binlinear').*ratio;
end

function im_warp = warpImage(im1, u, v, im2)
    [x,y] = meshgrid(1:size(im1,2),1:size(im1,1));
    shiftx = x + u;
    shifty = y + v;
    % perform warping, fill in nan for missing pixels
    im_warp = interp2(x,y,im2,shiftx,shifty,'linear',NaN);
    mask = isnan(im_warp);
    im_warp(mask) = im1(mask);
end



function ker = boxker(hsize)
    % box filter
    M = 2*hsize+1;
    ker = ones(M,M);
end

function ker = gauker(hsize)
    % gaussian filter
    sigma = 0.5;
    ker = gaussian_kernel_calculator(2, hsize, sigma);
    ker = ker./sum(abs(ker(:)));
end

function [flow_u, flow_v] = denseflowLK(im1, im2, iu, iv, hsize)
% compute optical flow using Lucas-Kanade
    % grad    
    [Ix, Iy] = grad2(im1,2,0.5);
    
    % dense flow
    width = size(im1,2);
    height = size(im1,1);
    [yy,xx] = meshgrid(1:height,1:width);
    xx = vec(xx');
    yy = vec(yy');
    
    if isempty(iu) 
        flow_u = zeros(height, width);
    else
        flow_u = iu;
    end
    
    if isempty(iv) 
    	flow_v = zeros(height, width);
    else
        flow_v = iv;
    end
    
    % precomputing
    Ix2 = Ix.*Ix;
    Iy2 = Iy.*Iy;
    Ixy = Ix.*Iy;
    % filtering
    ker = gauker(1);
    
    sIx2 = imfilter(Ix2,ker,'replicate','same');
    sIy2 = imfilter(Iy2,ker,'replicate','same');
    sIxy = imfilter(Ixy,ker,'replicate','same');
    sIxy2 = sIxy.*sIxy;

    % interpolation
    It = gradIntensity(xx,yy,vec(flow_u),vec(flow_v),im1,im2);
    It = reshape(It, height, width);
%     It = im2 - im1;
    
    Ixt = Ix.*It;
    Iyt = Iy.*It;
        
    sIxt = imfilter(Ixt,ker,'replicate','same');
    sIyt = imfilter(Iyt,ker,'replicate','same');
    
    % LK-flow    
    % opt1: vectorization, faster
    % conditioning: check if the smallest eigenvalue is very close to zero.
    % since it is a 2x2 positive semidefinite matrix, its eigen value has 
    % a analytical solution and must be a real value greater or equal to 0.
%     eig_smallest = 0.5.*(sIx2 + sIy2 - sqrt((sIx2-sIy2).^2+4.*sIxy2));
%     eig_largest = 0.5.*(sIx2 + sIy2 + sqrt((sIx2-sIy2).^2+4.*sIxy2));
%     ratio = eig_smallest ./ eig_largest;
%     invalid = ratio < 0.01 | (eig_smallest < 1e-6);
    
    % second option is to use Harris 
%     score = sIx2.*sIy2 - sIxy2 - 0.01.*(sIx2+sIy2).^2;
    
    % third option is to use the hormonic mean
    score = (sIx2.*sIy2 - sIxy2)./(sIx2+sIy2);
    det1 = sIx2.*sIy2 - sIxy2;
    invalid = score < 0.03*max(score(:)) | abs(det1) < 1e-6;
    
    % add a smaller value to the diagonal elements of those invalid pixels.
    %     sIx2(invalid) = sIx2(invalid) + 0.1;
    %     sIy2(invalid) = sIy2(invalid) + 0.1;

    % analytical inversion
    s = 1./(sIx2.*sIy2 - sIxy2);
    dflow_u = -( sIy2.*sIxt - sIxy.*sIyt).*s;
    dflow_v = -(-sIxy.*sIxt + sIx2.*sIyt).*s;
        
    dflow_u(invalid) = 0;
    dflow_v(invalid) = 0;
    
    flow_u = flow_u + dflow_u;
    flow_v = flow_v + dflow_v;
    
    % debug only
    showFlowQuiver(im1, flow_u, flow_v);
end

function [flow_u, flow_v] = denseflowHS(im1, im2)
% compute optical flow using Horn-Shunck method    
    % grad    
    [Ix, Iy] = grad2(im1,3,1);
    It = im2 - im1;
    
    % precomputing
    Ix2 = Ix.*Ix;
    Iy2 = Iy.*Iy;
    Ixy = Ix.*Iy;    
    Ixt = Ix.*It;
    Iyt = Iy.*It;
    
    flow_u = zeros(size(im1));
    flow_v = zeros(size(im1));
    
    % kernel 1, hard code 3x3 averaging
%     ker_avg = [1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];
    sigma = 0.5;
    hsize = 5;
    ker_avg = gaussian_kernel_calculator(2, hsize, sigma);
    
    max_iter = 100;
    alpha = 1;
    for iter = 1:max_iter
        % arveraging
        ubar = imfilter(flow_u,ker_avg,'replicate','same');
        vbar = imfilter(flow_v,ker_avg,'replicate','same');
        % update
        den = alpha*alpha + Ix2 + Iy2;
        flow_u = ubar - (Ix2.*ubar + Ixy.*vbar + Ixt)./den;
        flow_v = vbar - (Ixy.*ubar + Iy2.*vbar + Iyt)./den;
    end
    
end

