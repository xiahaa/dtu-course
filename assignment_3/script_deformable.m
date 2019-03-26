clc;close all;clear all;

addpath ./utils;

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
    stepSize = 20;
    a = 0.5;
    b = 0.5;
    
    % regularization matrix
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
            
        % find mean intensities inside and outside
        [cin, cout] = meanIntensity(im, curve, []);
        
        % compute the displacement along the normal direction
        displacement = computeDisplacement(im, curve, cin, cout);
        
        % regularization
%         Bint = regularization(a, b, size(curve,2));

        % find final
        displacement = displacement.*stepSize;
%         displacement = findDisplacement(im, curve, displacement);
        
        % draw normals
        quiver(curve(2,:),curve(1,:),displacement(2,:),displacement(1,:));
        
        % update
        curve = (Bint\(curve+displacement)')';
        
        % reinterpolation
        curve = suppressSelfIntersection(curve,m,n);
        curve = constraintCurve(curve, m, n);
        curve = reInterpolate(curve,Num);
        pause(0.1);
        % if write
        %F = getframe;
        %imwrite(F.cdata,strcat('../data/Ex_4_data/output_2/',num2str(i,'%d'),'.png'));
    end 
end

function displacement = findDisplacement(im, curve, displacement)
    % assume normals have already been signed
    [m,n] = size(im);
    
    finalPoint = (curve + displacement);
    startPoint = (curve);
    
    for i = 1:size(curve,2)
        % warp from normals to img
        x = [startPoint(1,i) finalPoint(1,i)];
        y = [startPoint(2,i) finalPoint(2,i)];
        nPoints = max(abs(diff(x)), abs(diff(y)))+1;  % Number of points in line
        
        xraw = linspace(x(1), x(2), nPoints);
        yraw = linspace(y(1), y(2), nPoints)
        cIndex = round(xraw);  % column indices
        rIndex = round(yraw);  % row indices
        
        id = cIndex > 0 & cIndex <= n & rIndex > 0 & rIndex <= m;
        rIndex = rIndex(id);
        cIndex = cIndex(id);
        xraw = xraw(id);
        yraw = yraw(id);
        
        try
            val = im((cIndex-1).*m+rIndex);
        catch
            error('');
        end
        % gradient and find the one with maximum gradient
        gval = gradient(val);
        [~,id] = max((gval));
        displacement(:,i) = [xraw(id)-curve(1,i);yraw(id)-curve(2,i)];
    end
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
