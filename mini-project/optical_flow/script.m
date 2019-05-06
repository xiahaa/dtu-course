clc;close all;clear all;

addpath utils;
addpath flowcore;

baseDir = '../../data/optical_flow_data/';

filename = {'Army','Backyard','Basketball','Dumptruck','Evergreen','Grove','Mequon', ...
            'Schefflera','Teddy','Urban','Wooden','Yosemite'};

fileid = 1;

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
[flow_u, flow_v] = denseflowPyrLK(im1, im2,1);
time1 = toc;

showFlowQuiver(im1, flow_u, flow_v);
title('Optical Flow LK: Quiver plot');
set(gca,'FontName','Arial','FontSize',20);
if save == true
    fig=gcf;                                     % your figure
    fig.PaperPositionMode='auto';
%     fig.PaperSize = [50 50];
%     fig.PaperUnits = 'centimeters';
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


