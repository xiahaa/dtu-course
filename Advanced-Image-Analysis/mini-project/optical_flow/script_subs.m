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
[flow_u1, flow_v1] = denseflowPyrLK(im1, im2, 0);

[flow_u2, flow_v2] = denseflowPyrHS(im1, im2,0);

im1 = im2double(imread(buildingScene.Files{1}));
im2 = im2double(imread(buildingScene.Files{2}));

im2_warp_lk = warpImage(im1, flow_u1, flow_v1, im2);
im2_warp_hs = warpImage(im1, flow_u2, flow_v2, im2);

% imb1 = im2_warp_lk.*0.8 + im1.*0.2;
% imb2 = im2_warp_hs.*0.8 + im1.*0.2;
% figure;imshow(imb1);
% figure;imshow(imb2);

imsubs1 = im2_warp_lk - im1;
imsubs2 = im2_warp_hs - im1;
imsubs1 = abs(imsubs1);
imsubs2 = abs(imsubs2);
for i = 1:3
    ims1 = imsubs1(:,:,i);
    ims2 = imsubs2(:,:,i);
    ims1 = (ims1 - min(ims1(:)))./(max(ims1(:))-min(ims1(:)));
    ims2 = (ims2 - min(ims2(:)))./(max(ims2(:))-min(ims2(:)));
    imsubs1(:,:,i) = ims1;
    imsubs2(:,:,i) = ims2;
end

figure;imshow(imsubs1);colorbar;
title('Image Diff: LK');
set(gca,'FontName','Arial','FontSize',20);
if save == true
    fig=gcf;                                     % your figure
%     fig.PaperPositionMode='auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperPosition = [0 0 20 20];
    print(strcat('./diff/',filename{fileid},'_LK_diff.png'),'-dpng','-r300');
end

figure;imshow(imsubs2);colorbar;
title('Image Diff: HS');
set(gca,'FontName','Arial','FontSize',20);
if save == true
    fig=gcf;                                     % your figure
%     fig.PaperPositionMode='auto';
    fig.PaperUnits = 'centimeters';
    fig.PaperPosition = [0 0 20 20];
    print(strcat('./diff/',filename{fileid},'_HS_diff.png'),'-dpng','-r300');
end



function [flow_u, flow_v] = denseflowPyrLK(im1,im2, verbose)
    layers = 4;
    ker = gaussian_kernel_calculator(2, 2, 1);
    
    pyr1 = GaussianPyramid(im1, layers, ker, verbose);
    pyr2 = GaussianPyramid(im2, layers, ker, verbose);
    
    if verbose == 1
        pyrflowu = cell(layers,1);
        pyrflowv = cell(layers,1);
    end
       
    % last layer
    [pu, pv] = denseflowLK(pyr1{end}, pyr2{end}, [], [], verbose);
    if verbose == 1
        pyrflowu{end} = pu;
        pyrflowv{end} = pu;
    end
    
    for i = layers-1:-1:1
        [u,v]=resampleFlow(pu,pv,size(pyr1{i}));
        [pu, pv] = denseflowLK(pyr1{i}, pyr2{i}, u, v, verbose);
        if verbose == 1
            pyrflowu{i} = pu;
            pyrflowv{i} = pu;
        end
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


