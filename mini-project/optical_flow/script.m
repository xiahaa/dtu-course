clc;close all;clear all;

addpath utils;

baseDir = '../../data/optical_flow_data/Basketball';

% load database
buildingScene = imageDatastore(baseDir);
numImages = numel(buildingScene.Files);

% load a certain image
im1 = imread(buildingScene.Files{1});
im1 = imPreprocessing(im1);

im2 = imread(buildingScene.Files{2});
im2 = imPreprocessing(im2);

% image size
[flow_u, flow_v] = denseflowPyrLK(im1, im2);

% display dense flow as an image
img = computeColor(flow_u,flow_v);
figure
if size(img,1)*size(img,2) < 1000
    imshow(img, 'InitialMagnification',1000);hold on;
else
    imshow(img);
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

function It = gradIntensity(x,y,flowu,flowv,im1,im2)
% perform intensity gradient for LK optical flow.
    shiftx = x + flowu;
    shifty = y + flowv;
    
    [xx,yy] = meshgrid(1:size(im1,2),1:size(im1,1));
    
    im_warp = interp2(xx,yy,im2,shiftx,shifty,'linear',NaN);
    mask = isnan(im_warp);
    im_warp(mask) = im1(mask);

    It = im_warp - im1((x-1)*size(im1,1)+y);
end

function im_warp = warpImage(im1, u, v, im2)
    [x,y] = meshgrid(1:size(im1,2),1:size(im1,1));
    shiftx = x + u;
    shifty = y + v;
    % perform warping, fill in nan for missing pixels
    im_warp = interp2(x,y,im2,shiftx,shifty,'linear',NaN);
    mask = isnan(im_warp);
    im_warp(mask) = im1(mask);
end

function [flow_u, flow_v] = denseflowPyrLK(im1,im2)
    layers = 4;
    ker = gaussian_kernel_calculator(2, 2, 1);
    
    pyr1 = GaussianPyramid(im1, layers, ker, 1);
    pyr2 = GaussianPyramid(im2, layers, ker, 1);
       
    % last layer
    [pu, pv] = denseflowLK(pyr1{end}, pyr2{end}, [], []);
    for i = layers-1:-1:1
        [u,v]=resampleFlow(pu,pv,size(pyr1{i}));
        [pu, pv] = denseflowLK(pyr1{i}, pyr2{i}, u, v);
    end
    flow_u = pu;
    flow_v = pv;
end

function [flow_u, flow_v] = denseflowLK(im1, im2, iu, iv)
% compute optical flow using Lucas-Kanade
    % dense flow
    width = size(im2,2);
    height = size(im2,1);
    if isempty(iu) 
        flow_u = zeros(height, width);
    else
        flow_u = iu;
    end
    
    if isempty(iv) 
    	flow_v = zeros(height, width);
    else
        flow_v = iv;
    end
    
    im2_warp = warpImage(im1, flow_u, flow_v, im2);
    % grad 
    [Ix, Iy] = grad2(im2_warp);
    It = im2_warp - im1;
    
    % precomputing
    Ix2 = Ix.*Ix;
    Iy2 = Iy.*Iy;
    Ixy = Ix.*Iy;
    % filtering
    ker = gauker(2, 1);% hsize + sigma
    
    sIx2 = imfilter(Ix2,ker,'symmetric','same');
    sIy2 = imfilter(Iy2,ker,'symmetric','same');
    sIxy = imfilter(Ixy,ker,'symmetric','same');
    sIxy2 = sIxy.*sIxy;
    
    Ixt = Ix.*It;
    Iyt = Iy.*It;

    sIxt = imfilter(Ixt,ker,'symmetric','same');
    sIyt = imfilter(Iyt,ker,'symmetric','same');
    
    % LK-flow    
    % opt1: vectorization, faster
    % conditioning: check if the smallest eigenvalue is very close to zero.
    % since it is a 2x2 positive semidefinite matrix, its eigen value has 
    % a analytical solution and must be a real value greater or equal to 0.
%     eig_smallest = 0.5.*(sIx2 + sIy2 - sqrt((sIx2-sIy2).^2+4.*sIxy2));
%     eig_largest = 0.5.*(sIx2 + sIy2 + sqrt((sIx2-sIy2).^2+4.*sIxy2));
%     ratio = eig_smallest ./ eig_largest;
%     invalid = ratio < 0.01 | (eig_smallest < 1e-6);
    
    % second option is to use Harris 
%     score = sIx2.*sIy2 - sIxy2 - 0.01.*(sIx2+sIy2).^2;
    
    % third option is to use the hormonic mean
    score = (sIx2.*sIy2 - sIxy2)./(sIx2+sIy2);
    det1 = sIx2.*sIy2 - sIxy2;
    invalid =  score < 0.05*max(score(:)) | abs(det1) < 1e-6;%
    % add a smaller value to the diagonal elements of those invalid pixels.
    %     sIx2(invalid) = sIx2(invalid) + 0.1;
    %     sIy2(invalid) = sIy2(invalid) + 0.1;

    % analytical inversion
    s = 1./(sIx2.*sIy2 - sIxy2);
    dflow_u = -( sIy2.*sIxt - sIxy.*sIyt).*s;
    dflow_v = -(-sIxy.*sIxt + sIx2.*sIyt).*s;
        
    dflow_u(invalid) = 0;
    dflow_v(invalid) = 0;
    
    flow_u = flow_u + dflow_u;
    flow_v = flow_v + dflow_v;
    
    % debug only
    showFlowQuiver(im1, flow_u, flow_v);
end

function [flow_u, flow_v] = denseflowPyrHS(im1,im2)
    layers = 4;
    ker = gaussian_kernel_calculator(2, 2, 1);% hsize x sigma
    
    pyr1 = GaussianPyramid(im1, layers, ker, 1);
    pyr2 = GaussianPyramid(im2, layers, ker, 1);
   
    % last layer
    [pu, pv] = denseflowHS(pyr1{end}, pyr2{end}, [], []);
    for i = layers-1:-1:1
        [u,v]=resampleFlow(pu,pv,size(pyr1{i}));
        [pu, pv] = denseflowHS(pyr1{i}, pyr2{i}, u, v);
        pu = medfilt2(pu,[5,5]);
        pv = medfilt2(pv,[5,5]);
    end
    flow_u = pu;
    flow_v = pv;
end

function [flow_u, flow_v] = denseflowHS(im1, im2, iu, iv)
% compute optical flow using Horn-Shunck method  
    height = size(im2,1);
    width = size(im2,2);
    if isempty(iu) 
        flow_u = zeros(height, width);
    else
        flow_u = iu;
    end
    
    if isempty(iv) 
    	flow_v = zeros(height, width);
    else
        flow_v = iv;
    end

    maxiter = 3;
    for iter = 1:maxiter
        im2_warp = warpImage(im1, flow_u, flow_v, im2);
        % grad    
        [Ix, Iy] = grad3(im2_warp);
        It = im2_warp - im1;
    
        % precomputing
        Ix2 = Ix.*Ix;
        Iy2 = Iy.*Iy;
        Ixy = Ix.*Iy;    
        Ixt = Ix.*It;
        Iyt = Iy.*It;
    
        % kernel 1, hard code 3x3 averaging
    %     ker_avg = [1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];
        sigma = 1;
        hsize = 2;
        ker_avg = gauker(2, 1);% hsize + sigma
    
        max_iter = 3;
        alpha = 1;
        
        uu = zeros(size(im2));
        vv = zeros(size(im2));
        
        for iter = 1:max_iter
            % arveraging
            ubar = imfilter(uu,ker_avg,'replicate','same');
            vbar = imfilter(vv,ker_avg,'replicate','same');
            % update
            den = alpha*alpha + Ix2 + Iy2;
            
            du = (Ix2.*ubar + Ixy.*vbar + Ixt)./den;
            dv = (Ixy.*ubar + Iy2.*vbar + Iyt)./den;
            
            uu = ubar - du;
            vv = vbar - dv;
            
            if max(abs(du(:))) < 1e-3 && max(abs(dv(:))) < 1e-3
                disp('Exit 1~~~~~~');
                break;
            end
        end
        flow_u = flow_u + uu;
        flow_v = flow_v + vv;
        if max(abs(uu(:))) < 1e-3 && max(abs(vv(:))) < 1e-3
            disp('Exit 2~~~~~~');
            break;
        end
    end
    showFlowQuiver(im1, flow_u, flow_v);
end

