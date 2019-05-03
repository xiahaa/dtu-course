clc;close all;clear all;

addpath ../../utils;

baseDir = '../../data/optical_flow_data';

% load database
buildingScene = imageDatastore(baseDir);
numImages = numel(buildingScene.Files);

% load a certain image
im1 = imread(buildingScene.Files{1});
im1 = imPreprocessing(im1);

im2 = imread(buildingScene.Files{2});
im2 = imPreprocessing(im2);

hsize = 3;
height = size(im1,1);
width = size(im1,2);

% image size
[flow_u, flow_v] = denseflowLK(im1, im2, hsize);
figure;imshow(im1, 'InitialMagnification',1000);hold on;

% opflow = opticalFlow(flow_u,flow_v);
% plot(opflow,'DecimationFactor',[1 1],'ScaleFactor',1);

[yy,xx] = meshgrid(1:height,1:width);
xx = vec(xx');
yy = vec(yy');
quiver(xx,yy,flow_u(:),flow_v(:),'LineWidth',1.5, 'Color','r','MaxHeadSize',1);axis image

% display dense flow as an image
img = computeColor(flow_u,flow_v);
figure
imshow(img, 'InitialMagnification',1000);hold on;

function [Ix, Iy] = grad1(im1)
% perform the simplest gradient calculation for LK optical flow.
    hx = [-1,0,1];
    Ix = imfilter(im1,hx,'replicate','same');
    Iy = imfilter(im1,hx','replicate','same');
end

function [Ix, Iy] = grad2(im1, hsize, sigma)
% perform gradient operation with first order gaussian kernel
    x = -hsize*sigma:1:hsize*sigma;
    cons1 = sigma*sigma;
    hg = 1/sqrt(2*pi*cons1).*exp(-x.^2./(2*cons1));
    hgx = 1/sqrt(2*pi*cons1).*exp(-x.^2./(2*cons1)).*(-x./cons1);
    hgx = fliplr(hgx);% convolution kernel is the flip version.
    Ix = imfilter(im1,hgx,'replicate','same');
    Ix = imfilter(Ix,hg','replicate','same');
    Iy = imfilter(im1,hgx','replicate','same');
    Iy = imfilter(Iy,hg,'replicate','same');
end

function It = gradIntensity(u,v,flowu,flowv,im1,im2)
% perform intensity gradient for LK optical flow.
    % new coorinate 
    xxnew = u + flowu;
    yynew = v + flowv;
    
    width = size(im2,2); height = size(im2,1);
    
    % nearest interpolation
    xx1 = round(xxnew); yy1 = round(yynew);% nearest neighbors
    
    xx1(xx1 < 1) = 1;
    xx1(xx1 > width) = width;
    yy1(yy1 < 1) = 1;
    yy1(yy1 > height) = height;
    
    % to index
    index0 = (u - 1).*size(im1,1) + v;
    index1 = (xx1 - 1).*height + yy1;
    
    It = im2(index1) - im1(index0);
end

function [flow_u, flow_v] = denseflowLK(im1, im2, hsize)
% compute optical flow using Lucas-Kanade
    % grad
    sigma = 1;
    [Ix, Iy] = grad2(im1,hsize,sigma);
    
    % dense flow
    width = size(im1,2);
    height = size(im1,1);
    [yy,xx] = meshgrid(1:height,1:width);
    xx = vec(xx');
    yy = vec(yy');
    
    flow_u = zeros(height, width);
    flow_v = zeros(height, width);
    
    % precomputing
    Ix2 = Ix.*Ix;
    Iy2 = Iy.*Iy;
    Ixy = Ix.*Iy;
    % filtering
    % box filter
%     M = 2*hsize+1;
%     ker = ones(M,M);
    
    % gaussian filter
    ker = gaussian_kernel_calculator(2, hsize, sigma);
    ker = ker./sum(ker(:));
    sIx2 = imfilter(Ix2,ker,'replicate','same');
    sIy2 = imfilter(Iy2,ker,'replicate','same');
    sIxy = imfilter(Ixy,ker,'replicate','same');
    sIxy2 = sIxy.*sIxy;
%     while true
        % interpolation
        It = gradIntensity(xx,yy,vec(flow_u),vec(flow_v),im1,im2);
        It = reshape(It, height, width);
    
        Ixt = Ix.*It;
        Iyt = Iy.*It;
        
        sIxt = imfilter(Ixt,ker,'replicate','same');
        sIyt = imfilter(Iyt,ker,'replicate','same');
        % LK-flow    
        % opt1: vectorization, faster
        % conditioning: check if the smallest eigenvalue is very close to zero.
        % since it is a 2x2 positive semidefinite matrix, its eigen value has 
        % a analytical solution and must be a real value greater or equal to 0.
        eig_smallest = 0.5.*(sIx2 + sIy2 - sqrt((sIx2-sIy2).^2+4.*sIxy2));
        eig_largest = 0.5.*(sIx2 + sIy2 + sqrt((sIx2-sIy2).^2+4.*sIxy2));
        ratio = eig_smallest ./ eig_largest;
        invalid = ratio < 0.01 | (eig_smallest < 1e-6);
    
        % add a smaller value to the diagonal elements of those invalid pixels.
    %     sIx2(invalid) = sIx2(invalid) + 0.1;
    %     sIy2(invalid) = sIy2(invalid) + 0.1;

        % analytical inversion
        s = 1./(sIx2.*sIy2 - sIxy2);
        dflow_u = -( sIy2.*sIxt - sIxy.*sIyt).*s;
        dflow_v = -(-sIxy.*sIxt + sIx2.*sIyt).*s;
        
        dflow_u(invalid) = 0;
        dflow_v(invalid) = 0;
        
        % TODO: add subpixel refinement
%         if max(abs(dflow_u(:))) < 0.1 && max(abs(dflow_v(:))) < 0.1 
%             break;
%         end
%         flow_u = flow_u + dflow_u;
%         flow_v = flow_v + dflow_v;
%     end
    flow_u = flow_u + dflow_u;
    flow_v = flow_v + dflow_v;
end

function pyr = GaussianPyramid(im, layers)
    pyr{1} = im;
    sigma = 1;
    for i = 2:layers
        % gaussian filter
        img = imgaussfilt(pyr{i-1}, sigma);
        % subsampling
        pyr{i-1} = imresize(img,0.5,'nearest');
    end
end

function [flow_u, flow_v] = denseflowHS(im1, im2, hsize)
% compute optical flow using Horn-Shunck method
    % grad
    sigma = 1;
%     [Ix, Iy] = grad2(im1,hsize,sigma);
    [Ix, Iy] = grad1(im1);
    
    % dense flow
    width = size(im1,2);
    height = size(im1,1);
    [yy,xx] = meshgrid(1:height,1:width);
    xx = vec(xx');
    yy = vec(yy');
    
    flow_u = zeros(height, width);
    flow_v = zeros(height, width);
    
    % precomputing
    Ix2 = Ix.*Ix;
    Iy2 = Iy.*Iy;
    Ixy = Ix.*Iy;
    % filtering
    % box filter
%     M = 2*hsize+1;
%     ker = ones(M,M);
    
    % gaussian filter, symmertric, no need to flip
    ker = gaussian_kernel_calculator(2, hsize, sigma);
    sIx2 = imfilter(Ix2,ker,'replicate','same');
    sIy2 = imfilter(Iy2,ker,'replicate','same');
    sIxy = imfilter(Ixy,ker,'replicate','same');
    sIxy2 = sIxy.*sIxy;

    % interpolation
    It = gradIntensity(xx,yy,vec(flow_u),vec(flow_v),im1,im2);
    It = reshape(It, height, width);
    
    Ixt = Ix.*It;
    Iyt = Iy.*It;
        
    sIxt = imfilter(Ixt,ker,'replicate','same');
    sIyt = imfilter(Iyt,ker,'replicate','same');
    % LK-flow    
    % opt1: vectorization, faster
    % conditioning: check if the smallest eigenvalue is very close to zero.
    % since it is a 2x2 positive semidefinite matrix, its eigen value has 
    % a analytical solution and must be a real value greater or equal to 0.
    eig_smallest = 0.5.*(sIx2 + sIy2 - sqrt((sIx2-sIy2).^2+4.*sIxy2));
    eig_largest = 0.5.*(sIx2 + sIy2 + sqrt((sIx2-sIy2).^2+4.*sIxy2));
    ratio = eig_smallest ./ eig_largest;
    invalid = ratio < 0.001 | (eig_smallest < 1e-6);
    
    % second option is to use Harris 
    % score = sIx2.*sIy2 - sIxy2 - 0.01.*(sIx2+sIy2).^2
    
    % third option is to use the hormonic mean
    % score = (sIx2.*sIy2 - sIxy2)./(sIx2+sIy2)

    % analytical inversion
    s = 1./(sIx2.*sIy2 - sIxy2);
    dflow_u = -( sIy2.*sIxt - sIxy.*sIyt).*s;
    dflow_v = -(-sIxy.*sIxt + sIx2.*sIyt).*s;

    dflow_u(invalid) = 0;
    dflow_v(invalid) = 0;
        
    flow_u = flow_u + dflow_u;
    flow_v = flow_v + dflow_v;
    %-------- this ends the initialization, then start HS iteration --------%
    
    % kernel 1, hard code 3x3 averaging
%     ker_avg = [1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];
    ker_avg = gaussian_kernel_calculator(2, hsize, sigma);
    
    max_iter = 50;
    alpha = 0.5;
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

