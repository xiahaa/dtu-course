clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

% add path to the data directory
addpath ../data/EX_2_data;
addpath ../utils/

skip = [ 3 4 5 6 7];
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
        tsqrt = tsqrt * sqrt(2); % scale update
        half_win_size = 1;
        [blobs_center] = detect_blobs(LL, half_win_size);
        blobs = [blobs_center repmat(sqrt(2)*tsqrt,size(blobs_center,1),1)];
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

%% ex2: computing contour length
if isempty(find(skip == 2,1))
    close all;
    imgnames = {'CT_lab_high_res.png', ...
             'CT_lab_low_res.png', ...
             'CT_lab_low_res.png', ...
             'CT_synchrotron.png', ...
             'Optical.png', ...
             'SEM.png'};
    % 4 5, need tuning
    for kk = 2%size(imgnames,1) 
        I0 = imread(imgnames{kk});
        I0 = im2double(I0);
        if size(I0,3) == 3
            Igray = rgb2gray(I0);
        else
            Igray = I0;
        end
%         Igray = imresize(Igray,[480,640]);
        [M,N] = size(Igray);
        K = 2;
        t0 = 1;
        LLNs = zeros(M,N,K);
        Ic2 = cat(3,Igray,Igray,Igray);
        ccmap = jet(128);
        w = linspace(0,2*pi,100);
        dtured = [153/255 0 0];
        for k = 1:K
            t = sqrt(2^(k-1)*t0);
            [LLN] = scale_normalized_LoG(Igray, t);
            LLNs(:,:,k) = LLN;
            half_win_size = 1;
            [blobs_center] = detect_blobs(LLNs(:,:,k), half_win_size);

            radius = sqrt(2)*t;
            rcw = radius*cos(w);
            rsw = radius*sin(w);

            num_blobs = size(blobs_center, 1);
            for i = 1:num_blobs
                xx = round(blobs_center(i,1) + rcw);
                yy = round(blobs_center(i,2) + rsw);
                valid = xx > 0 & xx <= M & yy > 0 & yy <= N;
                indices = sub2ind([M,N], xx(valid), yy(valid));
%                 nn = [-M +M -1 +1];
%                 indices = [indices indices+nn(1) indices+nn(2) indices+nn(3) indices+nn(4)];
                indices = indices(indices > 0 & indices < M*N);
                Ic2(indices) = dtured(1);%ccmap(i,1);
                Ic2(indices+M*N) = dtured(2);%ccmap(i,2);
                Ic2(indices+M*N*2) = dtured(3);%ccmap(i,3);
            end
        end
        figure
        imshow(Ic2,[]);
        title('$Multi-scale\ blob\ detection$','FontSize',fontSize,'Interpreter','latex');
    end
end
   
if isempty(find(skip == 3,1))
    close all;
    imgnames = {'CT_lab_high_res.png', ...
             'CT_lab_low_res.png', ...
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
        
        t = sqrt(2);%1
        sigma = 3;
        g = gassian_fast(t, sigma);
        
        Iblur = imfilter(Igray, g, 'replicate');
        Iblur = imfilter(Iblur, g', 'replicate');
        
%         threshold = 0.65*max(Iblur(:));
%         Ib = imbinarize(Iblur, threshold);
        imshow(Iblur,[]);
        
        half_win_size = 3;
        [local_maxima, Idebug] = find_maxima(Iblur, half_win_size);
        figure
        imshow(Idebug,[]);
        
        [M,N] = size(Igray);
        K = 8;
        t0 = 1;
        LLNs = zeros(M,N,K);
        scales = zeros(1,K);
        for k = 1:K
            t = sqrt(2)^(k-1)*t0;
            [LLN] = scale_normalized_LoG(Igray, t);
            LLNs(:,:,k) = LLN;
            scales(k) = t;
        end
        
        Ic2 = cat(3,Igray,Igray,Igray);
        ccmap = jet(128);
        w = linspace(0,2*pi,100);
        dtured = [153/255 0 0];
        
        for k = 1:size(local_maxima,1)
            val = LLNs(local_maxima(k,1),local_maxima(k,2),:);
            [~,minid] = min(val);
            scale = scales(minid);
            radius = sqrt(2)*scale;
            rcw = radius*cos(w);
            rsw = radius*sin(w);
            xx = round(local_maxima(k,1) + rcw);
            yy = round(local_maxima(k,2) + rsw);
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
    
function res = gassian_fast(t, sigma)
    x = round(-sigma*t):1:round(sigma*t);
    res = (2^(1/2).*exp(-x.^2./(2*t^2)))./(2*pi^(1/2)*(t^2)^(1/2));
end

function res = secondorder_der_gaussian_fast(t, sigma)
    x = round(-sigma*t):1:round(sigma*t);
    res = [(2^(1/2).*x.^2.*exp(-x.^2./(2*t^2)))./(2*t^4*pi^(1/2)*(t^2)^(1/2)) - ...
        (2^(1/2).*exp(-x.^2./(2*t^2)))./(2*t^2*pi^(1/2)*(t^2)^(1/2))];
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
    [xx,yy]=meshgrid([-half_win_size:1:half_win_size], ...
                     [-half_win_size:1:half_win_size]);
    xx = xx(:);
    yy = yy(:);
    
    idc = xx == 0 & yy == 0;
    xx(idc) = []; yy(idc) = [];
    
    [M,N] = size(I);
    
    maxcc = zeros(M*N,2);
    mincc = zeros(M*N,2);
    k1 = 0;
    k2 = 0;
    
    for i = 1:M
        subx = xx + i;
        for j = 1:N
            suby = yy + j;
            valid = subx>0 & subx < M & suby>0 & suby < N;
            validx = subx(valid);validy = suby(valid);
            indxy = sub2ind([M,N],validx,validy);
%             indxy = indxy(indxy > 0 & indxy <= M*N);
            maxval = max(I(indxy));
            if I(i,j) > maxval
                % maxima
                k1 = k1+1;
                maxcc(k1,:) = [i,j];
            else
                minval = min(I(indxy));
                if I(i,j) < minval
                    % minima
                    k2 = k2 + 1;
                    mincc(k2,:) = [i,j];
                end
            end
        end
    end
    maxcc(k1+1:end,:) = [];
    mincc(k2+1:end,:) = [];
    blobs_center = [maxcc;mincc];
end

function [varargout] = find_maxima(I, half_win_size)
    [xx,yy]=meshgrid([-half_win_size:1:-1,1:1:half_win_size], ...
                     [-half_win_size:1:-1,1:1:half_win_size]);
    xx = xx(:);
    yy = yy(:);
    [M,N] = size(I);
    
    maxcc = zeros(M*N,2);
    mincc = zeros(M*N,2);
    k1 = 0;
    k2 = 0;
    
    for i = 1:M
        subx = xx + i;
        for j = 1:N
            suby = yy + j;
            valid = subx>0 & subx < M & suby>0 & suby < N;
            validx = subx(valid);validy = suby(valid);
            indxy = sub2ind([M,N],validx,validy);
%             indxy = indxy(indxy > 0 & indxy <= M*N);
            maxval = max(I(indxy));
            if I(i,j) > maxval
                % maxima
                k1 = k1+1;
                maxcc(k1,:) = [i,j];
%             else
%                 minval = min(I(indxy));
%                 if I(i,j) < minval
%                     % minima
%                     k2 = k2 + 1;
%                     mincc(k2,:) = [i,j];
%                 end
            end
        end
    end
    maxcc(k1+1:end,:) = [];
    mincc(k2+1:end,:) = [];
    
    % further refinement
    indices = sub2ind(size(I), maxcc(:,1), maxcc(:,2));
    val = I(indices);
    meanval = mean(val);
    stdv = var(val);
    maxcc = maxcc(abs(val - meanval) < 5*stdv, :);
    
    blobs_center = maxcc;
%     % visulization
    cc = [maxcc];
    num_segments = size(cc, 1);
    ccmap = jet(num_segments);
    [M,N] = size(I);
    Idebug = cat(3,I,I,I);
    for i = 1:num_segments
        indices = sub2ind(size(I), cc(i,1), cc(i,2));
        Idebug(indices) = 1;%ccmap(i,1);
        Idebug(indices+M*N) = 0;%ccmap(i,2);
        Idebug(indices+M*N*2) = 0;%ccmap(i,3);
    end
%     imshow(Ic,[]);
    varargout{1} = blobs_center;
    if nargout == 2
        varargout{2} = Idebug;
    end
end
