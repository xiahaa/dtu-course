clc;close all;clear all;

addpath ./gui
addpath ./utils
addpath ./upload

skip = [1 3 4];

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
if sum((skip == 2)) == 0 
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

    q1(1:2,:) = q1(1:2,:);% + rand([2,size(q1,2)]);%*0.5-0.25;
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
        precons(:,i) = triangulationPoly(q1(:,i),P1,q2(:,i),P2);
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       sec assignment: triangulation with real images "dinosaur"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum((skip == 3)) == 0
    baseDir = './images/';
    buildingScene = imageDatastore(baseDir);
    numImages = numel(buildingScene.Files);

    load(strcat(baseDir,'viff.xy'));
    x = viff(1:end,1:2:72)';  % pull out x coord of all tracks
    y = viff(1:end,2:2:72)';  % pull out y coord of all tracks

    num = 0;
    for n = 1:numImages-1
        im1 = readimage(buildingScene, n);
        imshow(im1); hold on;
        id = x(n,:) ~= -1 & y(n,:) ~= -1;
        plot(x(n,id),y(n,id),'go'); 
        num = num + sum(id);
        hold off;
        pause(0.1);
    end

    load(strcat(baseDir,'dino_Ps.mat'));

    ptcloud = zeros(3,num);
    k = 1;
    for i = 1:size(x,1)-1
        % tracked features
        id = x(i,:) ~= -1 & y(i,:) ~= -1 & x(i+1,:) ~= -1 & y(i+1,:) ~= -1;
        q1 = [x(i,id);y(i,id);];
        q2 = [x(i+1,id);y(i+1,id);];
        q1 = tohomo(q1);
        q2 = tohomo(q2);
        P1 = P{i};
        P2 = P{i+1};
        % triangulation
        precons = zeros(3,size(q1,2));
        for j = 1:size(q1,2)
            precons(:,j) = triangulationNITER2(q1(:,j),P1,q2(:,j),P2);
        end
        ptcloud(:,k:k+size(q1,2)-1) = precons;
        k = k + size(q1,2);
    end
    figure
    plot3(ptcloud(1,:),ptcloud(2,:),ptcloud(3,:),'k.','MarkerSize',10);hold on;
    grid on;
    axis equal;
    view(3);
    for i = 1:size(x,1)
        P1 = P{i};
        [K, R, t, c] = decomposeP(P1);
        plotCamera('Location',c,'Orientation',R,'Opacity',0,'Color',[0 1 0],'Size',0.05);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       sec assignment: triangulation with stereo images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    I1 = imread('im0.png');
    I2 = imread('im1.png');
    J1 = rgb2gray(I1);
    J2 = rgb2gray(I2);
    disparityRange = [0 480];
    disparityMap = disparity(J1,J2,'DisparityRange',disparityRange,'UniquenessThreshold',20);
    figure
    imshow(disparityMap,disparityRange)
    title('Disparity Map')
    colormap jet
    colorbar
    save('disp.mat', 'disparityMap')
else
    I1 = imread('im0.png');
    I2 = imread('im1.png');
    figure;imshow(I1);
    
    load('disp.mat');
    disparityRange = [0 480];
    figure
    imshow(disparityMap,disparityRange)
    title('Disparity Map')
    colormap jet
    colorbar
    f = 4841.191;
    cx = 1441.607;
    cy = 969.311;
    b = 170.458*1e-3;
    id = disparityMap >= 100 & disparityMap <= 480;
    
    ptcloud = stereoTriangulate(disparityMap, f, cx, cy, b,id);
    
    colorcloud = uint8(zeros(size(ptcloud,1),3));
    r = I1(:,:,1); g = I1(:,:,2); b = I1(:,:,3);
    colorcloud(:,1) = r(id);
    colorcloud(:,2) = g(id);
    colorcloud(:,3) = b(id);
    ptCloud = pointCloud(ptcloud,'Color',colorcloud);
    pcshow(ptCloud);
end

function xyz = stereoTriangulate(disp, f, cx, cy, b,id)
    z = f * b ./ disp(id);
    [yy,xx] = find(id);
    x = (xx - cx).*z / f;
    y = (yy - cy).*z / f;
    xyz = [x y z];
end


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
    phomo = tohomo(p);
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








