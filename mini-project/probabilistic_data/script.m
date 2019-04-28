clc;close all;clear all;

addpath ../../utils;

baseDir = '../../data/probabilistic_data/';

imgid = 9;

% load database
buildingScene = imageDatastore(baseDir);
numImages = numel(buildingScene.Files);

type = 1;

% load a certain image
im = imread(buildingScene.Files{imgid});
im = imPreprocessing(im, type);

% image size
m = size(im,1);
n = size(im,2);
    
% parameters
Num = 500;
stepSize = 5;
a = 0;
b = 5;
    
K = 256;

% regularization matrix
Bint = regularization(a, b, Num);
[Bt,B] = calcB(im, K);

% curve initialization
alpha = linspace(0,2*pi, Num);
r = max(m,n) / 5;
curve(1,:) = m / 2 + r*cos(alpha);
curve(2,:) = n / 2 + r*sin(alpha);


%% start probabilistic deforming
    im = imread(buildingScene.Files{imgid});
    im = imPreprocessing(im, type);
    
    for i = 1:200
        % show
        imshow(im);hold on;
        % draw
        plot(curve(2,[1:end,1]),curve(1,[1:end,1]),'r-','LineWidth',2);
        % compute probability map
        [Pin,maskin,maskout] = calcProbabilityMap(im, curve, [], B, Bt);  
        % compute the displacement along the normal direction
        displacement = computeDisplacement(curve, Pin, maskin, maskout);
        % find final
        displacement = displacement.*stepSize;
        % draw normals
%         quiver(curve(2,:),curve(1,:),displacement(2,:),displacement(1,:));
        % update with regularization
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

function varargout = calcB(im, K)
% compute linear mapping matrix
    M = size(im,1);
    N = size(im,2);
    B = zeros(M*N,K);
    shift = M*N;
    for i = 1:K
        id = find(im == (i - 1));
        B(id + (i-1)*shift) = 1;
    end
    % make B a sparse matrix
    B = sparse(B);
    varargout{1} = B';
    varargout{2} = B;
end

function [Pin, maskin, maskout] = calcProbabilityMap(im, curve, boundary, B, Bt)
%Compute mean intensity for deformable models.
    polygon = curve;
    polygon(:,end+1) = curve(:,1);
    % region inside curve
    mask = poly2mask(polygon(2,:),polygon(1,:),size(im,1),size(im,2));
    % if boundary is empty consider all pixels inside region as in, otherwise as out
    if isempty(boundary)
        maskin = mask;
        maskout = ~mask;
    else
        % otherwise, find relevent in and out
        maskin = mask & boundary;
        maskout = xor(maskin, boundary);
    end
    
    % debug, visualization
%     im(mask) = 1;
%     imshow(im);

    fin = Bt*maskin(:) ./ sum(maskin(:));
    fout = Bt*maskout(:)./ sum(maskout(:));
    
    pin = fin ./ (fin + fout);
    Pin = B * pin;
    Pin = reshape(Pin,size(im,1),size(im,2));
end

function Icurve = biInterIntensity(im, curve)
%Performs bilinear interpolation of intensity for curve points.
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

function displacement = computeDisplacement(curve, Pin, maskin, maskout)
%Compute displacements for curves.
    pins = biInterIntensity(Pin, curve);
%     Pout = biInterIntensity(pout, curve);
%     pins = mean(Pin(maskin));
%     pouts = 1-mean(Pin(maskout));
    % scalar
    s = pins*2 - 1;
    % project to normal direction
    normals = computeNormal(curve);
    % compute displacement
    displacement = s.*normals;
end

function normals = computeNormal(curve)
%Naive method of computing normals
	% firstly, computing tangent direction
    tangent = zeros(2,size(curve,2));
    tangent(:,1) = curve(:,2) - curve(:,end);
    tangent(:,2:end-1) = curve(:,3:end) - curve(:,1:end-2);
    tangent(:,end) = curve(:,1) - curve(:,end-1);
    % normalize
    tangent = tangent ./ sqrt(tangent(1,:).^2+tangent(2,:).^2);
    % rotate to outward normal direction
    normals = [tangent(2,:);-tangent(1,:)];
end

function im = imPreprocessing(im, type)
    if type == 1
        if size(im,3) == 3
            im = rgb2gray(im);
        end
%         im = im2double(im);
    elseif type == 2
        r = im2double(im(:,:,1));
        g = im2double(im(:,:,2));
        b = im2double(im(:,:,3));
        imres = (2*b-(r+g)+2)/4;
        im = imres;
    end
end
