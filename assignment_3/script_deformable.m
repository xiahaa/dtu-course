clc;close all;clear all;

addpath ../data/EX_4_data;
addpath ../utils;

% videoname = 'crawling_amoeba.mov';
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
        imwrite(vidFrame,sprintf('../data/EX_4_data/echiniscus/im_%03d.png',id));% save one frame as
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
        plot(curve(2,[1:end,1]),curve(1,[1:end,1]),'r-');
        
        % find mean intensities inside and outside
        [cin, cout] = meanIntensity(im, curve);
        
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

function curve = suppressSelfIntersection(curve)
    [~,~,segments]=selfintersect(curve(2,:),curve(1,:)); 
    id = ones(1,size(curve,2));
    for i = 1:size(segments,1)
        id(segments(i,1)+1:segments(i,2)-1) = 0;
    end
    id = id == 1;
    curve = curve(:,id);
    curve = reInterpolate(curve);
end

function curve = constraintCurve(curve, m, n)
    id1 = curve(1,:) < 2;
    id2 = curve(1,:) > m-1;
    id3 = curve(2,:) < 2;
    id4 = curve(2,:) > n-1;
    curve(1,id1) = 2;
    curve(1,id2) = m - 1;
    curve(2,id3) = 2;
    curve(2,id4) = n - 1;
end

function curve = reInterpolate(curve)
    accumulateDist = zeros(1,size(curve,2));
    err = curve(:,[2:end]) - curve(:,1:end-1);
    dist = diag(err'*err)';
    for i = 2:size(curve,2)
        accumulateDist(i) = accumulateDist(i-1) + cumsum(dist(i-1));
    end
    accumulateDistNew = linspace(0,accumulateDist(end),size(curve,2));
    try
        xnew = interp1(accumulateDist,curve(1,:),accumulateDistNew);
    catch
        error('s');
    end
    ynew = interp1(accumulateDist,curve(2,:),accumulateDistNew);
    curve = [xnew;ynew];
end

function Bint = regularization(a, b, N)
    L1 = spdiags([-2.*ones(N,1) ones(N,1) ones(N,1)],[0,1,-1],N,N);
    L1(1,N) = 1;
    L1(N,1) = 1;
    L1 = L1.*a;
    L2 = spdiags([-1.*ones(N,1) 4.*ones(N,1) -6.*ones(N,1) 4.*ones(N,1) -1.*ones(N,1)], ...
                 [-2,-1,0,1,2],N,N);
    L2(1,N) = 4;L2(1,N-1)=-1;
    L2(2,N) = -1;
    L2(N,1) = 4;L2(N,2) = -1;
    L2(N-1,1) = -1;
    L2 = L2.*b;
    Bint = (eye(N) - (L1+L2));
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

function Icurve = biInterIntensity(im, curve)
    point1 = floor(curve);point2 = point1;point3 = point1;
    point2(2,:) = point2(2,:) + 1;
    point3(1,:) = point3(1,:) + 1;
    point4 = point2;
    point4(1,:) = point4(1,:) + 1;
    
    d1 = curve(1,:) - point1(1,:);
    d2 = curve(2,:) - point1(2,:);
    
    s1 = (1-d1).*(1-d2);
    s2 = (1-d1).*d2;
    s3 = d1.*(1-d2);
    s4 = (d1).*(d2);
    
    index1 = sub2ind(size(im), point1(1,:),point1(2,:));
    index2 = sub2ind(size(im), point2(1,:),point2(2,:));
    index3 = sub2ind(size(im), point3(1,:),point3(2,:));
    index4 = sub2ind(size(im), point4(1,:),point4(2,:));
    
    Icurve = im(index1).*s1 + ...
             im(index2).*s2 + ...
             im(index3).*s3 + ...
             im(index4).*s4;
end

function [cin, cout] = meanIntensity(im, curve)
    polygon = curve;
    polygon(:,end+1) = curve(:,1);
    mask = poly2mask(polygon(2,:),polygon(1,:),size(im,1),size(im,2));
    % debug
%     im(mask) = 1;
%     imshow(im);
    cin = mean(vec(im(mask)));
    cout = mean(vec(im(~mask)));
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
