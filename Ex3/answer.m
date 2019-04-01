clc;close all;clear all;

addpath ./gui
addpath ./utils

skip = [1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       first assignment: play with epipolar geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum((skip == 1)) == 0 
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       sec assignment: triangulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100;
p = rand([3,N]) * 10 - 5;% from near to far
p(1,:) = p(1,:) + 15;
C1 = [0;0;0];
C2 = [0;0.2;0];% baseline: 10cm

rad1 = -10/57.3;
R1 = [cos(rad1) 0 -sin(rad1);0 1 0;sin(rad1) 0 cos(rad1)]*[0 -1 0;0 0 -1;1 0 0];
R2 = [cos(-rad1) 0 -sin(-rad1);0 1 0;sin(-rad1) 0 cos(-rad1)]*[0 -1 0;0 0 -1;1 0 0];

t1 = -R1*C1;
t2 = -R2*C2;

figure
h1 = plot3(p(1,:),p(2,:),p(3,:),'g*');hold on;
cam1 = plotCamera('Location',C1,'Orientation',R1,'Opacity',0,'Color',[1 0 0],'Size',0.4,'Label','Camera1');
cam2 = plotCamera('Location',C2,'Orientation',R2,'Opacity',0,'Color',[0 1 0],'Size',0.4,'Label','Camera2');
axis equal
xlabel('x:(m)');
ylabel('y:(m)');
zlabel('z:(m)');
title('Triangulation Simulation');
set(gca,'FontName','Arial','FontSize',20);

K = [500 0 320;0 500 240;0 0 1];

[uv1, in1] = proj(R1,t1, p, K);
[uv2, in2] = proj(R2,t2, p, K);

in = in1 & in2;
% for i = 
q1 = uv1(:,in);
ptrue = p(:,in);
q2 = uv2(:,in);

q1(1:2,:) = q1(1:2,:) + rand([2,size(q1,2)]);%*0.5-0.25;
q2(1:2,:) = q2(1:2,:);% + rand([2,size(q2,2)]);

% im = zeros(K(2,3)*2,K(1,3)*2);
% figure;imshow(im);hold on;
% plot(q1(2,:),q1(1,:),'ro');
% figure;imshow(im);hold on;
% plot(q2(2,:),q2(1,:),'go');

P1 = K*[R1 t1];
P2 = K*[R2 t2];

precons = zeros(3,size(q1,2));
for i = 1:size(q1,2)
    precons(:,i) = triangulationMidpoint(q1(:,i),P1,q2(:,i),P2);
end
figure
plot3(ptrue(1,:),ptrue(2,:),ptrue(3,:),'ro','MarkerSize',8);hold on;
plot3(precons(1,:),precons(2,:),precons(3,:),'g+','MarkerSize',8);hold on;
xlabel('x:(m)');
ylabel('y:(m)');
zlabel('z:(m)');
title('Triangulation Simulation');
legend({'Truth','Reconstruction'});
set(gca,'FontName','Arial','FontSize',20);
grid on;

function [uv1, in] = proj(R, t, p, K)
% function [uv1, in] = proj(R, t, p, K)
%
% Projection of 3D points.
%   Inputs:
%       R,t: extrinsics.
%       K: camera intrinsics.
%       p: 3D points.
%   Outputs:
%       uv1: projected image pixels.
%       in??? array of indicators. Each indicator indicate whether this point is inside the image: 1 inside, 0 otherwise.
% Author: xiahaa@space.dtu.dk
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
    P1 = K*([R t]);
    phomo = tohomogeneous(p);
    uv1 = P1*phomo;
    uv1 = uv1./uv1(3,:);
    in = uv1(1,:) > 0 & uv1(1,:) < K(1,3)*2 & uv1(2,:) > 0 & uv1(2,:) < K(2,3)*2;
end


function P = triangulationMidpoint(x1,P1,x2,P2)
%Implementation of mid-point triangulation method.
    M1 = P1(1:3,1:3);
    c1 = -M1\P1(1:3,4);
    
    M2 = P2(1:3,1:3);
    c2 = -M2\P2(1:3,4);
    
    % ray
    r1 = M1\x1;
    r2 = M2\x2;
    
    A = [r1 -r2];
    b = c2 - c1;
    
    x = A\b;
    
    P = (c1 + x(1).*r1 + c2 + x(2).*r2)*0.5;
end








