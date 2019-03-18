clc;close all;clear all;

I = imread('../Ex1/Tiles_perspective_distort.png');
addpath ../Ex1/

I = im2double(I);
% I = rgb2gray(I);
h = figure;
imshow(I);
[x,y] = ginput(4);
q = round([x y]');
q = [q;ones(1,4)];
I = drawlines(I,q,[[1 2];[3 4];[1 4];[2 3]]);
imshow(I);

% 
l1 = cross(q(:,1),q(:,2));
l2 = cross(q(:,3),q(:,4));
l3 = cross(q(:,1),q(:,4));
l4 = cross(q(:,2),q(:,3));

%
p1 = cross(l1,l2);p1 = p1./p1(3);
p2 = cross(l3,l4);p2 = p2./p2(3);

%
linf = cross(p1,p2);
linf = linf ./ linf(3);

H = [1 0 0;0 1 0;linf(1) linf(2) linf(3)];

% H = [1.1 -1.1 0;1.1 1.1 0;0 0 1];
% e
% H = [cos(-5/57.3) sin(-5/57.3) 40;-sin(-5/57.3) cos(-5/57.3) 100;0 0 1];
% s
s=1;
%H = [s*cos(-5/57.3) s*sin(-5/57.3) 0;-s*sin(-5/57.3) s*cos(-5/57.3) 0;0 0 1];
% affine
% H = [s*cos(-5/57.3) s*sin(-5/57.3) 2;-s*sin(-5/57.3) s*cos(-5/57.3) 2;0.001 0.001 1];
tic
Iw = warpping(I,H);
toc
tic
tf = projective2d(H);
Iw1 = imwarp(I,tf');
toc
if s == 1
    if size(Iw,3) == 3
        Iwr = Iw(:,:,1);Iwr = imresize(Iwr,[size(I,1),size(I,2)]);
        Iwg = Iw(:,:,2);Iwg = imresize(Iwg,[size(I,1),size(I,2)]);
        Iwb = Iw(:,:,3);Iwb = imresize(Iwb,[size(I,1),size(I,2)]);
        Iw = cat(3,Iwr,Iwg,Iwb);
    else
        Iw = imresize(Iw,size(I));
    end
    Ig = cat(2,I,Iw);
end
imshow(Ig)