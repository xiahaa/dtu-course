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
[flow_u, flow_v] = flowLK(im1, im2, hsize);
figure;imshow(uint8(im1), 'InitialMagnification',1500);hold on;
[yy,xx] = meshgrid(1:height,1:width);
xx = vec(xx');
yy = vec(yy');
quiver(xx,yy,flow_u(:),flow_v(:),'LineWidth',1.5, 'Color','r','MaxHeadSize',1);axis image


function [Ix, Iy, It] = grad1(im1,im2)
% perform the simplest gradient calculation for LK optical flow.
    hx = [-1,0,1];
    Ix = imfilter(im1,hx,'replicate','same');
    Iy = imfilter(im1,hx','replicate','same');
    It = -im2 + im1;% for later computation
end


function [flow_u, flow_v] = flowLK(im1, im2, hsize)
% compute optical flow using Lucas-Kanade
    % grad
    [Ix, Iy, It] = grad1(im1,im2);
    Ix2 = Ix.*Ix;
    Iy2 = Iy.*Iy;
    Ixy = Ix.*Iy;
    Ixt = Ix.*It;
    Iyt = Iy.*It;
    % box filter
    M = 2*hsize+1;
    ker = ones(M,M)./(M*M);
    % filtering
    sIx2 = imfilter(Ix2,ker,'replicate','same');
    sIy2 = imfilter(Iy2,ker,'replicate','same');
    sIxy = imfilter(Ixy,ker,'replicate','same');
    sIxt = imfilter(Ixt,ker,'replicate','same');
    sIyt = imfilter(Iyt,ker,'replicate','same');
    % LK-flow    
    % opt1: vectorization, faster
    
    
    
    s = 1./(sIx2.*sIy2 - sIxy.^2 + 1e-15);
    flow_u = ( sIy2.*sIxt - sIxy.*sIyt).*s;
    flow_v = (-sIxy.*sIxt + sIx2.*sIyt).*s;
    
    % opt2: use loop, 3-4 times lower than vectorization
%     flow_u = zeros(size(im1));
%     flow_v = zeros(size(im1));
%     tic
%     for i = 1:size(im1,1)
%         for j = 1:size(im1,2)
%             A = [sIx2(i,j) sIxy(i,j);sIxy(i,j) sIy2(i,j)];
%             b = [sIxt(i,j);sIyt(i,j)];
%             res = A\b;
%             flow_u(i,j) = res(1);
%             flow_v(i,j) = res(2);
%         end
%     end
%     toc
end


function im = imPreprocessing(im)
    if size(im,3) == 3
        im = rgb2gray(im);
    end
    im = single(im);
end
