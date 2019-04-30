clc;close all;clear all;

addpath ../../utils;

baseDir = '../../data/probabilistic_data/';

imgid = 8;

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
stepSize = 2;
a = 2;
b = 5;
    
K = 256;

% regularization matrix
Bint = regularization(a, b, Num);

% curve initialization
alpha = linspace(0,2*pi, Num);
r = max(m,n) / 20;

imshow(im);
[xx,yy] = ginput(1);

curve(1,:) = yy + r*cos(alpha);
curve(2,:) = xx + r*sin(alpha);

usePatchProbability = true;
if usePatchProbability == false
    [Bt,B] = calcB(im, K);
else
    M = 9;
    L = 1554;
    [Bt,B] = calcPatchB(im, M, L,imgid);
%     for i = 1:size(B,1)
        Cinv = 1./(sum(B,2)+1e-15);
%     end
end

%% start probabilistic deforming
    im = imread(buildingScene.Files{imgid});
    im = imPreprocessing(im, type);
    
    for i = 1:200
        % show
        figure(1);imshow(im);hold on;
        % draw
        plot(curve(2,[1:end,1]),curve(1,[1:end,1]),'r-','LineWidth',2);
        hold off;
        % compute probability map
        [Pin,maskin,maskout] = calcProbabilityMap(im, curve, [], B, Bt, Cinv);  
        figure(2);imshow(Pin,[]);
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
%         curve = reInterpolate(curve,Num);
        pause(0.1);
        % if write
        %F = getframe;
        %imwrite(F.cdata,strcat('../data/Ex_4_data/output_2/',num2str(i,'%d'),'.png'));
    end
    
function varargout = calcPatchB(im, M, L, imgid)
    W = size(im,2);
    H = size(im,1);
    % form pathches
    numPatches = ((W - M) + 1) * ((H-M) + 1);
    patches = zeros(numPatches, M*M);
    patchesMap = zeros(numPatches, M*M);
    hsize = floor(M*0.5);
    
    [lx,ly] = meshgrid(-hsize:1:hsize,-hsize:1:hsize);
    lx = vec(lx'); ly = vec(ly');
    k = 1;
    for i = hsize+1:1:H-hsize
        for j = hsize+1:1:W-hsize
            % local coordinate
            llx = i + lx;
            lly = j + ly;
            % to index
            index = (lly-1).*H + llx;
            patches(k,:) = im(index)';
            % form map
            patchesMap(k,:) = index';
            k = k + 1;
        end
    end
    
    % kmeans
    if 0
        opts = statset('MaxIter',100);        
        [patchClusterID, clusterMean] = kmeans(patches, L, 'Distance', 'sqeuclidean','Options', opts);
        
        % form B
        M2 = M*M;
        is = zeros(size(patchClusterID,1)*M2,1);
        js = zeros(size(patchClusterID,1)*M2,1);
        ss = ones(size(patchClusterID,1)*M2,1);
        k = 1;
        for i = 1:size(patchClusterID,1)
            % which cluster
            ii = (patchClusterID(i)-1)*M2;
            for j = 1:M2
                % which entry
                jj = patchesMap(i,j);
                is(k) = ii + j;
                js(k) = jj;
                k = k + 1;
            end
        end
        Bt = sparse(is,js,ss);
        
        save(strcat('res',num2str(imgid),'.mat'),'patchClusterID','clusterMean','Bt');
    else
        load(strcat('res',num2str(imgid),'.mat'));
    end
    
    varargout{1} = Bt;
    varargout{2} = Bt';
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

function [Pin, maskin, maskout] = calcProbabilityMap(im, curve, boundary, B, Bt, Cinv)
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
    
    pin = fin ./ (fin + fout + 1e-12);
    if ~isempty(Cinv)
        Pin = Cinv.*(B * pin);
    else
        Pin = B * pin;
    end
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

function Icurve = nearestInterIntensity(im, curve)
    point1 = round(curve);
    id1 = point1(1,:) < 0;
    id2 = point1(1,:) > size(im,1);
    id3 = point1(2,:) < 0;
    id4 = point1(2,:) > size(im,2);
    point1(1,id1) = 1;
    point1(1,id2) = size(im,1);
    point1(2,id3) = 1;
    point1(2,id4) = size(im,2);
    Icurve = im((point1(2,:)-1).*size(im,1)+point1(1,:));
end

function displacement = computeDisplacement(curve, Pin, maskin, maskout)
%Compute displacements for curves.
    pins = nearestInterIntensity(Pin, curve);
    
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
