clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

% add path to the data directory
addpath ../data/EX_2_data;
addpath ../utils/

skip = [1 2 3 ];
drawGif = 0;

% since we will work with second order derivation, we need this template
sigma = 3;% window size, 3t ~ 99%
t = 1;% standard deviation

% dgg = derivative_gaussian1d_generator(t, sigma, 2);% 2 order 
% g = derivative_gaussian1d_generator(t, sigma, 0);% 0 order

% fast solution
g = gassian_fast(t, sigma);
dgg = secondorder_der_gaussian_fast(t, sigma);

fontSize = 15;

if isempty(find(skip == 1,1))
    %% ex1: single scale detection
    I0 = imread('test_blob_uniform.png');
    % to double
    I0 = im2double(I0);
    % if color, to gray
    if size(I0,3) == 3
        Igray = rgb2gray(I0);
    else
        Igray = I0;
    end
    % viz
    imshow(I0,[]);
    % gaussian filter
    Ixx = imfilter(Igray, dgg, 'replicate');% along x
    Ixx = imfilter(Ixx, g', 'replicate');% smooth y, since g=gx*gy, dg=dgx gy + gx dgy. partial dev only pick one item
    Iyy = imfilter(Igray, dgg', 'replicate');% along y
    Iyy = imfilter(Iyy, g, 'replicate');%along x
    % combine, show
    figure
    Ixxyy = img_merge(Ixx,Iyy);
    imshow(Ixxyy,[]);
    title('$2^{nd}\ derv\ of\ Gaussian\ conv\ with\ I$','FontSize',fontSize,'Interpreter','latex');

    % LoG
    LL = Ixx+Iyy;
    figure
    imshow(LL,[])
    title('$LOG\ of\ I$','FontSize',fontSize,'Interpreter','latex');
    
    % blob detection: local minima and maxima with a threshold that the
    % response should be at least 0.1*maximum response
    half_win_size = 1;
    [blobs_center] = detect_blobs(LL, half_win_size);
    % debug only
    num_blobs = size(blobs_center, 1);
    % color map
    ccmap = jet(num_blobs);
    [M,N] = size(Igray);
    % color image for drawing
    Ic = cat(3,Igray,Igray,Igray);
    % radius
    radius = sqrt(2)*t;% since in my code, t is standard deviation
    w = linspace(0,2*pi,100);
    % precompute of the drawing pixels locally
    rcw = radius*cos(w);
    rsw = radius*sin(w);
    for i = 1:num_blobs
        xx = round(blobs_center(i,1) + rcw);
        yy = round(blobs_center(i,2) + rsw);
        valid = xx > 0 & xx <= M & yy > 0 & yy <= N;
        indices = sub2ind([M,N], xx(valid), yy(valid));
        % if too thin, you can draw the 4 neighbors
%         nn = [-M +M -1 +1];% draw also 4 neighbors
%         indices = [indices];% indices+nn(1) indices+nn(2) indices+nn(3) indices+nn(4)];
%         indices = indices(indices > 0 & indices < M*N);
        Ic(indices) = ccmap(i,1);
        Ic(indices+M*N) = ccmap(i,2);
        Ic(indices+M*N*2) = ccmap(i,3);
    end
    figure
    imshow(Ic,[]);
    title('$Single-scale\ blob\ detection$','FontSize',fontSize,'Interpreter','latex');
end


if isempty(find(skip == 2,1))
    %% ex2: multiple scales detections
    I0 = imread('test_blob_varying.png');
    % to double
    I0 = im2double(I0);
    % if color, to gray
    if size(I0,3) == 3
        Igray = rgb2gray(I0);
    else
        Igray = I0;
    end
    
    K = 10;
    t0 = 2^(1/2);
    [LLNs, radius] = create_scale_normalized_LoG(Igray, t0, K);
%         half_win_size = 1;
%         blobs = detect_blobs(LLNs, half_win_size);
    LLMax = max(abs(LLNs),[],3);
    blobs = [];
    
    for i = 1:1:1
        half_win_size = 9;
        [blobs_center] = detect_blobs(LLMax, half_win_size);
        blobs = [blobs;blobs_center repmat(radius(i),size(blobs_center,1),1)];
    end
    % drawing
    [M,N] = size(Igray);
    Ic = cat(3,Igray,Igray,Igray);
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
        Ic(indices) = dtured(1);%ccmap(i,1);
        Ic(indices+M*N) = dtured(2);%ccmap(i,2);
        Ic(indices+M*N*2) = dtured(3);%ccmap(i,3);
    end
    figure
    imshow(Ic,[]);
    title('$Multi-scale\ blob\ detection$','FontSize',fontSize,'Interpreter','latex');
end

%% ex3: detecting blobs in real data
if isempty(find(skip == 3,1))
    close all;
    imgnames = {'CT_lab_high_res.png', ...
             'CT_lab_med_res.png', ...
             'CT_lab_low_res.png', ...
             'CT_synchrotron.png', ...
             'Optical.png', ...
             'SEM.png'};
    for kk = 5%size(imgnames,1) 
        I0 = imread(imgnames{kk});
        I0 = im2double(I0);
        if size(I0,3) == 3
            Igray = rgb2gray(I0);
        else
            Igray = I0;
        end
        % if the image is too large, then resize a little bit
%         Igray = imresize(Igray,[480,640]);
        [M,N] = size(Igray);
        
        K = 8;
        t0 = 2^(1/2);
        [LLNs, radius] = create_scale_normalized_LoG(Igray, t0, K);
%         half_win_size = 1;
%         blobs = detect_blobs(LLNs, half_win_size);
        LLMax = min(abs(LLNs),[],3);

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
            Ic2(indices) = dtured(1);%ccmap(i,1);
            Ic2(indices+M*N) = dtured(2);%ccmap(i,2);
            Ic2(indices+M*N*2) = dtured(3);%ccmap(i,3);
        end
        figure
        imshow(Ic2,[]);
        title('$Multi-scale\ blob\ detection$','FontSize',fontSize,'Interpreter','latex');
    end
end
   
if isempty(find(skip == 4,1))
    close all;
    imgnames = {'CT_lab_high_res.png', ...
             'CT_lab_med_res.png', ...
             'CT_lab_low_res.png', ...
             'CT_synchrotron.png', ...
             'Optical.png', ...
             'SEM.png'};
    smooth_scales = [4 3 1 4 4 4];
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
        
        t = smooth_scales(kk);%1 4 2 3 3 1 4 4, 5 4,6 4
        sigma = 3;
        g = gassian_fast(t, sigma);
        Iblur = imfilter(Igray, g, 'replicate');
        Iblur = imfilter(Iblur, g', 'replicate');
        
%         threshold = 0.65*max(Iblur(:));
%         Ib = imbinarize(Iblur, threshold);
        imshow(Iblur,[]);
        
        half_win_size = 2;
        [local_maxima, Idebug] = find_maxima(Iblur, half_win_size);
        figure
        imshow(Idebug,[]);
        
        [M,N] = size(Igray);
        K = 6;
        t0 = 2^(1/5);
        [LLNs, radius] = create_scale_normalized_LoG(Igray, t0, K);
%         blobs = [];
%         for i = 1:1:K
%             half_win_size = 1;
%             [blobs_center] = detect_blobs(LLNs(:,:,i), half_win_size);
%             blobs = [blobs_center repmat(radius(i),size(blobs_center,1),1)];
%         end
        LLMax = min((LLNs),[],3);
        %% drawing
        Ic2 = cat(3,Igray,Igray,Igray);
        w = linspace(0,2*pi,100);
        dtured = [153/255 0 0];
        for i = 1:size(local_maxima,1)
            val = LLNs(local_maxima(i,1),local_maxima(i,2),:);
            [~,maxid] = min((val));
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
        title('Multi-scale blob detection','FontSize',fontSize,'FontName','Arial');
    end
end


