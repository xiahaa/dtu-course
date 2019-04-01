clc;close all;clear all;

addpath ./gui

%% loading
im1 = rgb2gray(imread('000002.bmp'));
im2 = rgb2gray(imread('000003.bmp'));

load('gopro_27.mat');

[im1,~] = undistortImage(im1,cameraParams);
[im2,~] = undistortImage(im2,cameraParams);

% im1 = imresize(im1,0.5);
% im2 = imresize(im2,0.5);

estimateF = 0;

if estimateF == 1
    load('Fdata.mat');
    
    im3 = cat(2,im1,im2);
    figure;imshow(im3);hold on;

    plot(x1(:,1),x1(:,2),'ro');
    plot(x2(:,1)+size(im1,2),x2(:,2),'go');

    shift = size(im1,2);
    cmap = lines(5);
    k = 1;
    for i = 1:size(x1,1)
        ptdraw = [x1(i,1), x1(i,2);
                  x2(i,1)+shift, x2(i,2)];
        plot(ptdraw(:,1),ptdraw(:,2),'LineStyle','-','LineWidth',1,'Color',cmap(k,:));
        k = mod(k+1,5);if k == 0 k = 1;end
    end
    F = DLT8pt(x1',x2');
    save('Fest.mat','F');
else
    load('Fest.mat');
    vgg_gui_F(im1,im2,F');
end

function xh = tohomo(x)
%From inhomogeneous coordinate to homogeneous coordiante.
    xh = [x;ones(1,size(x,2))];
end

function F = DLT8pt(x1,x2)
%Implementation of 8 point algorithm for fundamental matrix estimation.
    if size(x1,1) ~= 3
        x1h = tohomo(x1);
    end
    if size(x2,1) ~= 3
        x2h = tohomo(x2);
    end
    
    A = zeros(size(x1,2),9);
    for i = 1:size(x1,2)
        A(i,:) = kron(x1h(:,i)', x2h(:,i)');
    end
    
    [~,~,V] = svd(A);
    F = V(:,end);
    F = F([1 4 7;2 5 8;3 6 9]);
    F = F./F(3,3);
end





