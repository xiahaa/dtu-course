clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

% add path to the data directory
addpath ../data/EX_2_data;
addpath ../utils/

skip = [1 2 4 5 6 7];
drawGif = 0;

% since we will work with second order derivation, we need this template
sigma = 3;% window size, 3t ~ 99%
t = 0.5;% standard deviation

% dgg = derivative_gaussian1d_generator(t, sigma, 2);% 2 order 
% g = derivative_gaussian1d_generator(t, sigma, 0);% 0 order

% fast solution
g = gassian_fast(t, sigma);
dgg = secondorder_der_gaussian_fast(t, sigma);

fontSize = 15;

if isempty(find(skip == 1,1))
    %% ex1: image smoothing
    I0 = imread('test_blob_varying.png');
    I0 = im2double(I0);
    if size(I0,3) == 3
        Igray = rgb2gray(I0);
    else
        Igray = I0;
    end
    imshow(I0,[]);
    % gaussian filter
    Ixx = imfilter(Igray, dgg, 'replicate');% along x
    Ixx = imfilter(Ixx, g', 'replicate');% smooth y, since g=gx*gy, dg=dgx gy + gx dgy. partial dev only pick one item
    Iyy = imfilter(Igray, dgg', 'replicate');% along y
    Iyy = imfilter(Iyy, g, 'replicate');%along x

    figure
    Ixxyy = img_merge(Ixx,Iyy);
    imshow(Ixxyy,[]);
    title('$2^{nd}\ derv\ of\ Gaussian\ conv\ with\ I$','FontSize',fontSize,'Interpreter','latex');

    % 
    LL = Ixx+Iyy;
    figure
    imshow(LL,[])
    title('$LOG\ of\ I$','FontSize',fontSize,'Interpreter','latex');

    half_win_size = 1;
    [blobs_center] = detect_blobs(LL, half_win_size);
    % debug only
    num_blobs = size(blobs_center, 1);
    ccmap = jet(num_blobs);
    [M,N] = size(Igray);
    Ic = cat(3,Igray,Igray,Igray);
    radius = sqrt(2)*t;% since in my code, t is standard deviation
    w = linspace(0,2*pi,100);
    rcw = radius*cos(w);
    rsw = radius*sin(w);
    for i = 1:num_blobs
        xx = round(blobs_center(i,1) + rcw);
        yy = round(blobs_center(i,2) + rsw);
        valid = xx > 0 & xx <= M & yy > 0 & yy <= N;
        indices = sub2ind([M,N], xx(valid), yy(valid));
        nn = [-M +M -1 +1];
        indices = [indices indices+nn(1) indices+nn(2) indices+nn(3) indices+nn(4)];
        indices = indices(indices > 0 & indices < M*N);
        Ic(indices) = ccmap(i,1);
        Ic(indices+M*N) = ccmap(i,2);
        Ic(indices+M*N*2) = ccmap(i,3);
    end
    figure
    imshow(Ic,[]);
    title('$Single-scale\ blob\ detection$','FontSize',fontSize,'Interpreter','latex');

    % 1st imlementation, always start from 0
    K = 6;
    tsqrt = t;% prev scale
    prevLayer = Igray;
    blobs = [];
    for k = 1:K
        % do layer from the previsou layer
        % gaussian filter
        Ixx = imfilter(prevLayer, dgg, 'replicate');% along x
        Ixx = imfilter(Ixx, g', 'replicate');% smooth y, since g=gx*gy, dg=dgx gy + gx dgy. partial dev only pick one item
        Iyy = imfilter(prevLayer, dgg', 'replicate');% along y
        Iyy = imfilter(Iyy, g, 'replicate');%along x
        LL = Ixx+Iyy;
        LL = LL * tsqrt.^2;
        half_win_size = 1;
        [blobs_center] = detect_blobs(LL, half_win_size);
        blobs = [blobs_center repmat(sqrt(2)*tsqrt,size(blobs_center,1),1)];
        tsqrt = tsqrt * sqrt(2); % scale update
    end
    %% drawing
    w = linspace(0,2*pi,100);
    dtured = [153/255 0 0];
    for i = 1:size(blobs,1)
        radius = blobs(i,3);
        rcw = radius*cos(w);
        rsw = radius*sin(w);
        xx = round(blobs(i,1) + rcw);
        yy = round(blobs(i,2) + rsw);
        valid = xx > 0 & xx <= M & yy > 0 & yy <= N;
        indices = sub2ind([M,N], xx(valid), yy(valid));
        nn = [-M +M -1 +1];
        indices = [indices indices+nn(1) indices+nn(2) indices+nn(3) indices+nn(4)];
        indices = indices(indices > 0 & indices < M*N);
        Ic(indices) = dtured(1);%ccmap(i,1);
        Ic(indices+M*N) = dtured(2);%ccmap(i,2);
        Ic(indices+M*N*2) = dtured(3);%ccmap(i,3);
    end
    figure
    imshow(Ic,[]);
    title('$Multi-scale\ blob\ detection$','FontSize',fontSize,'Interpreter','latex');
end

%% ex2: TODO
if isempty(find(skip == 2,1))
    close all;
    imgnames = {'CT_lab_high_res.png', ...
             'CT_lab_med_res.png', ...
             'CT_lab_low_res.png', ...
             'CT_synchrotron.png', ...
             'Optical.png', ...
             'SEM.png'};
    % 4 5, need tuning
    for kk = 1%size(imgnames,1) 
        I0 = imread(imgnames{kk});
        I0 = im2double(I0);
        if size(I0,3) == 3
            Igray = rgb2gray(I0);
        else
            Igray = I0;
        end
%         Igray = imresize(Igray,[480,640]);
        [M,N] = size(Igray);
        K = 6;
        t0 = 2^(1/2);
        [LLNs, radius] = create_scale_normalized_LoG(Igray, t0, K);
%         half_win_size = 1;
%         blobs = detect_blobs(LLNs, half_win_size);
        LLMax = max(abs(LLNs),[],3);

        blobs = [];
        for i = 1:1:1
            half_win_size = 1;
            [blobs_center] = detect_blobs(LLMax, half_win_size);
            blobs = [blobs;blobs_center repmat(radius(i),size(blobs_center,1),1)];
        end
        %% drawing
        Ic2 = cat(3,Igray,Igray,Igray);
        w = linspace(0,2*pi,100);
        dtured = [153/255 0 0];
        for i = 1:size(blobs,1)
            val = LLNs(blobs(i,1),blobs(i,2),:);
            [~,maxid] = max(abs(val));
            r = radius(maxid);
            rcw = r*cos(w);
            rsw = r*sin(w);
            xx = round(blobs(i,1) + rcw);
            yy = round(blobs(i,2) + rsw);
            valid = xx > 0 & xx <= M & yy > 0 & yy <= N;
            indices = sub2ind([M,N], xx(valid), yy(valid));
            nn = [-M +M -1 +1];
            indices = [indices indices+nn(1) indices+nn(2) indices+nn(3) indices+nn(4)];
            indices = indices(indices > 0 & indices < M*N);
            Ic2(indices) = dtured(1);%ccmap(i,1);
            Ic2(indices+M*N) = dtured(2);%ccmap(i,2);
            Ic2(indices+M*N*2) = dtured(3);%ccmap(i,3);
        end
        figure
        imshow(Ic2,[]);
        title('$Multi-scale\ blob\ detection$','FontSize',fontSize,'Interpreter','latex');
    end
end
   
if isempty(find(skip == 3,1))
    close all;
    imgnames = {'CT_lab_high_res.png', ...
             'CT_lab_med_res.png', ...
             'CT_lab_low_res.png', ...
             'CT_synchrotron.png', ...
             'Optical.png', ...
             'SEM.png'};
    % 4 5, need tuning
    for kk = 6%size(imgnames,1) 
        I0 = imread(imgnames{kk});
        I0 = im2double(I0);
        if size(I0,3) == 3
            Igray = rgb2gray(I0);
        else
            Igray = I0;
        end
%         Igray = imresize(Igray,[480,640]);
        [M,N] = size(Igray);
        
        t = 4;%1 4 2 3 3 1 4 4, 5 4,6 4
        sigma = 3;
        g = gassian_fast(t, sigma);
        Iblur = imfilter(Igray, g, 'replicate');
        Iblur = imfilter(Iblur, g', 'replicate');
        
%         threshold = 0.65*max(Iblur(:));
%         Ib = imbinarize(Iblur, threshold);
        imshow(Iblur,[]);
        
        half_win_size = 1;
        [local_maxima, Idebug] = find_maxima(Iblur, half_win_size);
        figure
        imshow(Idebug,[]);
        
        [M,N] = size(Igray);
        K = 10;
        t0 = 2^(1/5);
        [LLNs, radius] = create_scale_normalized_LoG(Igray, t0, K);
%         blobs = [];
%         for i = 1:1:K
%             half_win_size = 1;
%             [blobs_center] = detect_blobs(LLNs(:,:,i), half_win_size);
%             blobs = [blobs_center repmat(radius(i),size(blobs_center,1),1)];
%         end
        LLMax = max(abs(LLNs),[],3);
        %% drawing
        Ic2 = cat(3,Igray,Igray,Igray);
        w = linspace(0,2*pi,100);
        dtured = [153/255 0 0];
        for i = 1:size(local_maxima,1)
            val = LLNs(local_maxima(i,1),local_maxima(i,2),:);
            [~,maxid] = max(abs(val));
%             maxid
            r = radius(maxid);
            rcw = r*cos(w);
            rsw = r*sin(w);
            xx = round(local_maxima(i,1) + rcw);
            yy = round(local_maxima(i,2) + rsw);
            valid = xx > 0 & xx <= M & yy > 0 & yy <= N;
            indices = sub2ind([M,N], xx(valid), yy(valid));
%             nn = [-M +M -1 +1];
%             indices = [indices indices+nn(1) indices+nn(2) indices+nn(3) indices+nn(4)];
            indices = indices(indices > 0 & indices < M*N);
            Ic2(indices) = dtured(1);%ccmap(i,1);
            Ic2(indices+M*N) = dtured(2);%ccmap(i,2);
            Ic2(indices+M*N*2) = dtured(3);%ccmap(i,3);
        end
        
        figure
        imshow(Ic2,[]);
        title('$Multi-scale\ blob\ detection$','FontSize',fontSize,'Interpreter','latex');
    end
end
    


function res = gassian_fast(t, sigma)
    x = round(-sigma*t):1:round(sigma*t);
    res = (2^(1/2).*exp(-x.^2./(2*t^2)))./(2*pi^(1/2)*(t^2)^(1/2));
end

function res = secondorder_der_gaussian_fast(t, sigma)
    x = round(-sigma*t):1:round(sigma*t);
    res = [(2^(1/2).*x.^2.*exp(-x.^2./(2*t^2)))./(2*t^4*pi^(1/2)*(t^2)^(1/2)) - ...
        (2^(1/2).*exp(-x.^2./(2*t^2)))./(2*t^2*pi^(1/2)*(t^2)^(1/2))];
end

function varargout = create_scale_normalized_LoG(Igray, t0, K)
%     tsqrt = t0;% prev scale
%     prevLayer = Igray;
%     LLs = zeros(size(Igray,1),size(Igray,2),K);
%     radius = zeros(1,K);
% 
%     sigma = 3;
%     g = gassian_fast(tsqrt, sigma);
%     dgg = secondorder_der_gaussian_fast(tsqrt, sigma);
%     
%     for k = 1:K
%         % do layer from the previsou layer
%         % gaussian filter
%         Ixx = imfilter(prevLayer, dgg, 'replicate');% along x
%         Ixx = imfilter(Ixx, g', 'replicate');% smooth y, since g=gx*gy, dg=dgx gy + gx dgy. partial dev only pick one item
%         Iyy = imfilter(prevLayer, dgg', 'replicate');% along y
%         Iyy = imfilter(Iyy, g, 'replicate');%along x
%         LL = Ixx+Iyy;
%         prevLayer = LL;
%         LL = LL * (tsqrt.^2);
%         LLs(:,:,k) = LL;
%         tsqrt = tsqrt * sqrt(2); % scale update
%         radius(k) = tsqrt;
%     end
    tsqrt = t0;% prev scale
    LLs = zeros(size(Igray,1),size(Igray,2),K);
    radius = zeros(1,K);
    for k = 1:K
        [LL] = scale_normalized_LoG(Igray, tsqrt);
        LLs(:,:,k) = LL;
        tsqrt = tsqrt * sqrt(2); % scale update
        radius(k) = tsqrt;
    end

    varargout{1} = LLs;
    varargout{2} = radius;
end

function [LL] = scale_normalized_LoG(I, tsqrt)
    sigma = 3;
    g = gassian_fast(tsqrt, sigma);
    dgg = secondorder_der_gaussian_fast(tsqrt, sigma);
    % gaussian filter
    Ixx = imfilter(I, dgg, 'replicate');% along x
    Ixx = imfilter(Ixx, g', 'replicate');% smooth y, since g=gx*gy, dg=dgx gy + gx dgy. partial dev only pick one item
    Iyy = imfilter(I, dgg', 'replicate');% along y
    Iyy = imfilter(Iyy, g, 'replicate');%along x
    % LoG
    LL = Ixx+Iyy;
    % normalized LoG
    LL = LL .* (tsqrt^2);
end

function [blobs_center] = detect_blobs(I, half_win_size)
%     I = abs(I);
%     thresh = 0.2*max(I(:));
    I = abs(I);
    val = I(:);
    meanval = mean(val);
    stdv = var(val);
    thresh = meanval+1*stdv;%0.3*max(I(:));
    mx = imdilate(I, strel('square',2*half_win_size+1));

    % Make mask to exclude points within radius of the image boundary. 
    bordermask = zeros(size(I));
    bordermask(half_win_size+1:end-half_win_size, half_win_size+1:end-half_win_size) = 1;
    
    % Find maxima, threshold, and apply bordermask
    cimmx = (I==mx) & (I>thresh) & bordermask;
    [r,c] = find(cimmx);        % Find row,col coords.
    blobs_center = [r c];
end

function [varargout] = find_maxima(I, half_win_size)
%     I = abs(I);
    val = I(:);
    meanval = mean(val);
    stdv = var(val);
    thresh = meanval+2*stdv;%0.3*max(I(:));
    mx = imdilate(I, strel('square',2*half_win_size+1));

    % Make mask to exclude points within radius of the image boundary. 
    bordermask = zeros(size(I));
    bordermask(half_win_size+1:end-half_win_size, half_win_size+1:end-half_win_size) = 1;
    
    % Find maxima, threshold, and apply bordermask
    cimmx = (I==mx) & (I>thresh) & bordermask;
    [r,c] = find(cimmx);        % Find row,col coords.
    blobs_center = [r c];
    
    % visulization
    cc = blobs_center;
    num_segments = size(cc, 1);
    [M,N] = size(I);
    Idebug = cat(3,I,I,I);
    for i = 1:num_segments
        indices = sub2ind(size(I), cc(i,1), cc(i,2));
        Idebug(indices) = 1;
        Idebug(indices+M*N) = 0;
        Idebug(indices+M*N*2) = 0;
    end
    varargout{1} = blobs_center;
    if nargout == 2
        varargout{2} = Idebug;
    end
end
