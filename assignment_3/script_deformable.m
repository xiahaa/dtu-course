clc;close all;clear all;

addpath ../data/EX_4_data;
addpath ../utils;

videoname = 'crawling_amoeba.mov';
% videoname = 'echiniscus.mp4';
data_dir = '/../data/EX_4_data/';

debug = 1;

if debug == 0
    % load video
    vidobj = VideoReader(videoname);
    % Create an axes
    currAxes = axes;
    % Read video frames until available
    id = 0;
    while hasFrame(vidobj)
        vidFrame = readFrame(vidobj);
        image(vidFrame, 'Parent', currAxes);
        currAxes.Visible = 'off';
        pause(1/vidobj.FrameRate);
        imwrite(vidFrame,sprintf('../data/EX_4_data/crawling_amoeba/im_%03d.png',id));% save one frame as
        %     the sample data for tuning
        id = id + 1;
    end
else
    typename = {'crawling_amoeba', 'echiniscus'};
    type = 1;
    
    imgpath = strcat(fullfile(pwd),strcat(data_dir,typename{type},'/'));  
    imgFiles = dir(fullfile(imgpath,'*.png'));
    imgFileNames = {imgFiles.name}';
    im = imread(strcat(imgpath,imgFileNames{1}));
    images = zeros(size(imgFileNames,1), size(im,1),size(im,2));    
    m = size(im,1);
    n = size(im,2);
    
    % parameters
    Num = 500;
    stepSize = 30;
    a = 0.5;
    b = 0.5;
    
    Bint = regularization(a, b, Num);
    
    % curve initialization
    alpha = linspace(0,2*pi, Num);
    r = max(m,n) / 5;
    curve(1,:) = m / 2 + r*cos(alpha);
    curve(2,:) = n / 2 + r*sin(alpha);
    
    
    for i = 1:size(imgFileNames,1)
        im = imread(strcat(imgpath,imgFileNames{i}));
        im = imPreprocessing(im, type);
        % show
        imshow(im);hold on;
        % draw
        plot(curve(2,[1:end,1]),curve(1,[1:end,1]),'r-','LineWidth',2);
        
        mask = poly2mask(curve(2,:), curve(1,:), m, n);
        boundary = findBoundary(mask);
        se = strel('disk',5);
        boundary = imdilate(boundary, se);
            
        % find mean intensities inside and outside
        [cin, cout] = meanIntensity(im, curve, boundary);
        
        % compute the displacement along the normal direction
        displacement = computeDisplacement(im, curve, cin, cout);
        
        % draw normals
        quiver(curve(2,:),curve(1,:),displacement(2,:),displacement(1,:));
        
        % regularization
        Bint = regularization(a, b, size(curve,2));

        % update
        curve = (Bint\(curve+displacement.*stepSize)')';
        
        % reinterpolation
        curve = reInterpolate(curve);
        curve = suppressSelfIntersection(curve);
        curve = constraintCurve(curve, m, n);
        
        pause(0.1);
    end
    
    
end

% poisson blending
function boundary = findBoundary(bw)
    boundary = zeros(size(bw));
    for i = 2:size(bw,1)-1
        for j = 2:size(bw,2)-1
            block = bw(i-1:i+1,j-1:j+1);
            if block(5) == 1 && sum(block([1,2,3,4,6,7,8,9])==0)~=0
                boundary(i,j) = 1;
            end
        end
    end
end


function displacement = computeDisplacement(im, curve, cin, cout)
    Icurve = biInterIntensity(im, curve);
    % scalar
    s = (cin - cout).*(2*Icurve - cin - cout);
    % project to normal direction
    normals = computeNormal(curve);
    % debug;
%     imshow(im);hold on;plot(curve(2,:),curve(1,:),'r-');
%     quiver(curve(2,:),curve(1,:),normals(2,:),normals(1,:));
    displacement = s.*normals;
end

function normals = computeNormal(curve)
    tangent = zeros(2,size(curve,2));
    tangent(:,1) = curve(:,2) - curve(:,end);
    tangent(:,2:end-1) = curve(:,3:end) - curve(:,1:end-2);
    tangent(:,end) = curve(:,1) - curve(:,end-1);
    % normalize
    tangent = tangent ./ sqrt(tangent(1,:).^2+tangent(2,:).^2);
    % rotate to outward normal direction
    normals = [tangent(2,:);-tangent(1,:)];
end

function [cin, cout] = meanIntensity(im, curve, boundary)
    polygon = curve;
    polygon(:,end+1) = curve(:,1);
    mask = poly2mask(polygon(2,:),polygon(1,:),size(im,1),size(im,2));
    
    if isempty(boundary)
        maskin = mask;
        maskout = ~mask;
    else
        maskin = mask & boundary;
        maskout = ~mask & boundary;
    end
    
    % debug
%     im(mask) = 1;
%     imshow(im);
    cin = mean(vec(im(maskin)));
    cout = mean(vec(im(maskout)));
end

function im = imPreprocessing(im, type)
    if type == 1
        if size(im,3) == 3
            im = rgb2gray(im);
        end
        im = im2double(im);
    elseif type == 2
        r = im2double(im(:,:,1));
        g = im2double(im(:,:,2));
        b = im2double(im(:,:,3));
        imres = (2*b-(r+g)+2)/4;
        im = imres;
    end
end
