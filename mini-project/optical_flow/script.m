clc;close all;clear all;

addpath utils;
addpath flowcore;

baseDir = '../../data/optical_flow_data/';

filename = {'Army','Backyard','Basketball','Dumptruck','Evergreen','Grove','Mequon', ...
            'Schefflera','Teddy','Urban','Wooden','Yosemite'};

fileid = 12;

dir = strcat(baseDir, filename{fileid});
% load database
buildingScene = imageDatastore(dir);
numImages = numel(buildingScene.Files);

% load a certain image
im1 = imread(buildingScene.Files{1});
im1 = imPreprocessing(im1);

im2 = imread(buildingScene.Files{2});
im2 = imPreprocessing(im2);

save = 1;

% flow
tic
[flow_u, flow_v] = denseflowPyrLK(im1, im2,0);
time1 = toc;

showFlowQuiver(im1, flow_u, flow_v);
title('Optical Flow LK: Quiver plot');
set(gca,'FontName','Arial','FontSize',20);
if save == true
    fig=gcf;                                     % your figure
    fig.PaperPositionMode='auto';
    print(strcat('./output/',filename{fileid},'_LK_quiver.pdf'),'-dpdf','-fillpage');
end

showFlowImage(flow_u, flow_v);
title('Optical Flow LK: Image plot');
set(gca,'FontName','Arial','FontSize',20);
if save == true
    fig=gcf;                                     % your figure
    fig.PaperPositionMode='auto';
    print(strcat('./output/',filename{fileid},'_LK_rgb.pdf'),'-dpdf','-fillpage');
end

tic
[flow_u, flow_v] = denseflowPyrHS(im1, im2,0);
time2 = toc;

showFlowQuiver(im1, flow_u, flow_v);
title('Optical Flow HS: Quiver plot');
set(gca,'FontName','Arial','FontSize',20);
if save == true
    fig=gcf;                                     % your figure
    fig.PaperPositionMode='auto';
    print(strcat('./output/',filename{fileid},'_HS_quiver.pdf'),'-dpdf','-fillpage');
end

showFlowImage(flow_u, flow_v);
title('Optical Flow HS: Image plot');
set(gca,'FontName','Arial','FontSize',20);
if save == true
    fig=gcf;                                     % your figure
    fig.PaperPositionMode='auto';
    print(strcat('./output/',filename{fileid},'_HS_rgb.pdf'),'-dpdf','-fillpage');
end

fid = fopen(strcat(dir,'/time.txt'),'w');
fprintf(fid, 'LK: %10.6f\n', time1);
fprintf(fid, 'LK: %10.6f\n', time2);
fclose(fid);

function showFlowImage(flow_u, flow_v)
    % display dense flow as an image
    img = computeColor(flow_u,flow_v);
    figure
    if size(img,1)*size(img,2) < 10000
        imshow(img, 'InitialMagnification',10000);hold on;
    else
        imshow(img);
    end
end

function showFlowQuiver(im, flow_u, flow_v)
    figure;
    if size(im,1)*size(im,2) < 10000
        imshow(im, 'InitialMagnification',1000);hold on;
    else
        imshow(im);hold on;
    end
    height = size(im,1);
    width = size(im,2);
    
    [xx1,yy1] = meshgrid(1:width,1:height);
    
    display_grid_height = round(height*0.01);
    display_grid_width = round(width*0.01);
    
    [xx2,yy2] = meshgrid(1:display_grid_width:width, 1:display_grid_height:height);
    uu = interp2(xx1,yy1,flow_u, xx2, yy2);
    vv = interp2(xx1,yy1,flow_v, xx2, yy2);
    
%     opflow = opticalFlow(uu,vv);
%     plot(opflow,'DecimationFactor',[1 1],'ScaleFactor',1);
%     [yy,xx] = meshgrid(1:height,1:width);
%     xx = vec(xx');
%     yy = vec(yy');
    quiver(xx2(:),yy2(:),uu(:),vv(:),4,'LineWidth',1, 'Color','r','MarkerSize',50);axis image
end

function [flow_u, flow_v] = denseflowPyrLK(im1,im2, verbose)
    layers = 4;
    ker = gaussian_kernel_calculator(2, 2, 1);
    
    pyr1 = GaussianPyramid(im1, layers, ker, verbose);
    pyr2 = GaussianPyramid(im2, layers, ker, verbose);
       
    % last layer
    [pu, pv] = denseflowLK(pyr1{end}, pyr2{end}, [], [], verbose);
    for i = layers-1:-1:1
        [u,v]=resampleFlow(pu,pv,size(pyr1{i}));
        [pu, pv] = denseflowLK(pyr1{i}, pyr2{i}, u, v, verbose);
    end
    flow_u = pu;
    flow_v = pv;
end

function [flow_u, flow_v] = denseflowPyrHS(im1,im2,verbose)
    layers = 4;
    ker = gaussian_kernel_calculator(2, 2, 1);% hsize x sigma
    
    pyr1 = GaussianPyramid(im1, layers, ker, verbose);
    pyr2 = GaussianPyramid(im2, layers, ker, verbose);
   
    % last layer
    [pu, pv] = denseflowHS(pyr1{end}, pyr2{end}, [], [], verbose);
    for i = layers-1:-1:1
        [u,v]=resampleFlow(pu,pv,size(pyr1{i}));
        [pu, pv] = denseflowHS(pyr1{i}, pyr2{i}, u, v, verbose);
        pu = medfilt2(pu,[5,5]);
        pv = medfilt2(pv,[5,5]);
    end
    flow_u = pu;
    flow_v = pv;
end


