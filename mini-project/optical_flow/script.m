clc;close all;clear all;

addpath ../../utils;

baseDir = '../../data/optical_flow_data/';

% load database
buildingScene = imageDatastore(baseDir);
numImages = numel(buildingScene.Files);

% load a certain image
im1 = imread(buildingScene.Files{1});
im1 = imPreprocessing(im1);

im2 = imread(buildingScene.Files{2});
im2 = imPreprocessing(im2);

hsize = 1;
height = size(im1,1);
width = size(im1,2);

% image size
[flow_u, flow_v] = denseflowLK(im1, im2, hsize);
figure;imshow(uint8(im1), 'InitialMagnification',1500);hold on;
[yy,xx] = meshgrid(1:height,1:width);
xx = vec(xx');
yy = vec(yy');
quiver(xx,yy,flow_u(:),flow_v(:),'LineWidth',1.5, 'Color','r','MaxHeadSize',1);axis image


function [Ix, Iy] = grad1(im1)
% perform the simplest gradient calculation for LK optical flow.
    hx = [-1,0,1];
    Ix = imfilter(im1,hx,'replicate','same');
    Iy = imfilter(im1,hx','replicate','same');
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
    
    It = -im2(index1) + im1(index0);% for later computation
end

function [flow_u, flow_v] = denseflowLK(im1, im2, hsize)
% compute optical flow using Lucas-Kanade
    % grad
    [Ix, Iy] = grad1(im1);
    
    % dense flow
    width = size(im1,2);
    height = size(im1,1);
    [yy,xx] = meshgrid(1:width,1:height);
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
    M = 2*hsize+1;
    ker = ones(M,M)./(M*M);
    sIx2 = imfilter(Ix2,ker,'replicate','same');
    sIy2 = imfilter(Iy2,ker,'replicate','same');
    sIxy = imfilter(Ixy,ker,'replicate','same');
    
%     while true
        % interpolation
        It = gradIntensity(xx,yy,vec(flow_u),vec(flow_v),im1,im2);
        It = reshape(It,height, width);
    
        Ixt = Ix.*It;
        Iyt = Iy.*It;
        
        sIxt = imfilter(Ixt,ker,'replicate','same');
        sIyt = imfilter(Iyt,ker,'replicate','same');
        % LK-flow    
        % opt1: vectorization, faster
        % conditioning: check if the smallest eigenvalue is very close to zero.
        % since it is a 2x2 positive semidefinite matrix, its eigen value has 
        % a analytical solution and must be a real value greater or equal to 0.
        eig_smallest = 0.5.*(sIx2 + sIy2 - sqrt((sIx2-sIy2).^2+4.*sIxy.^2));
        eig_largest = 0.5.*(sIx2 + sIy2 + sqrt((sIx2-sIy2).^2+4.*sIxy.^2));
        ratio = eig_smallest ./ eig_largest;
        invalid = ratio < 0.05;
    
        % add a smaller value to the diagonal elements of those invalid pixels.
    %     sIx2(invalid) = sIx2(invalid) + 0.1;
    %     sIy2(invalid) = sIy2(invalid) + 0.1;

        % analytical inversion
        s = 1./(sIx2.*sIy2 - sIxy.^2 + 1e-15);
        dflow_u = ( sIy2.*sIxt - sIxy.*sIyt).*s;
        dflow_v = (-sIxy.*sIxt + sIx2.*sIyt).*s;
        
        dflow_u(invalid) = 0;
        dflow_v(invalid) = 0;
        
%         if max(abs(dflow_u(:))) < 0.1 && max(abs(dflow_v(:))) < 0.1 
%             break;
%         end
        flow_u = flow_u + dflow_u;
        flow_v = flow_v + dflow_v;
%     end
    flow_u = flow_u + dflow_u;
    flow_v = flow_v + dflow_v;
end


function im = imPreprocessing(im)
    if size(im,3) == 3
        im = rgb2gray(im);
    end
    im = single(im);
end
