clc;close all;clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

% add path to the data directory
% addpath ./EX_1_data;
addpath ../utils/;
data_dir = '../data/EX_1_data/';

skip = [1 2  4 5 6 7];
drawGif = 0;

if isempty(find(skip == 1,1))
    %% ex1: image smoothing
    I0 = imread(strcat(data_dir,'fibres_xcth.png'));
    I0 = im2double(I0);
    imshow(I0,[]);
    % gaussian filter
    % 1st try with 2D Gaussian filter
    gxy = gaussian_kernel_calculator(2, 2, 3);
    Ixy = imfilter(I0,gxy,'replicate');
    % 2nd try with two seperate 1D Gaussian filter
    gx = gaussian_kernel_calculator(1, 2, 3);
    Ix = imfilter(I0,gx,'replicate');% along x direction
    
    If = convImg(I0, gx);

    
    Ix = imfilter(Ix,gx','replicate');% again along y direction
    disp(max(abs(Ixy(:)-Ix(:))));
    If = img_merge(Ixy,Ix);

    
    imshow(If,[]);
    title('Images: left (raw) and right (Gaussian Smoothed)', 'FontSize', 15, 'Interpreter','latex');
    TV0 = total_variantion_calculator(I0);
    TV1 = total_variantion_calculator(Ixy);
    disp(['TV original: ', num2str(TV0), '; TV filtered: ', num2str(TV1)]);
end

%% ex2: computing contour length
if isempty(find(skip == 2,1))
    close all;
    I0 = imread(strcat(data_dir,'fuel_cell_3.tif'));
    segments = img_segmentation(I0);
    % visulization
    num_segments = size(segments, 2);
    ccmap = jet(num_segments);
    Ic = cat(3, I0, I0, I0);
    Ic = im2double(Ic);
    % imshow(Ic,[]);
    Icc = Ic;
    [M,N] = size(I0);
    for i = 1:num_segments
        segment = segments{i};
        indices = sub2ind(size(I0), segment.point_set(1,:), segment.point_set(2,:));
        Icc(indices) = ccmap(i,1);
        Icc(indices+M*N) = ccmap(i,2);
        Icc(indices+M*N*2) = ccmap(i,3);
        if segment.boundary_length > 100
            center = round(mean(segment.point_set,2));
            center = flipud(center);
            Icc = insertText(Icc,center',num2str(segment.boundary_length),'FontSize',20);
        end
    end
    Im = img_merge(Ic,Icc);
    imshow(Im,[]);
    title('Computing boundary length', 'FontSize', 15, 'Interpreter','latex');
end

%% ex3: curve smoothing
if isempty(find(skip == 3,1))
    dino = load(strcat(data_dir,'dino.txt'));
    dino_noisy = load(strcat(data_dir,'dino_noisy.txt'));
    figure(1)
    subplot(2,2,1);
    plot(dino(:,1),dino(:,2),'-', 'LineWidth', 3);
    legend({'Truth'},'FontSize', 15,'Interpreter','latex');
    subplot(2,2,2);
    plot(dino_noisy(:,1),dino_noisy(:,2),'r-', 'LineWidth', 2);hold on;
    ps1 = curve_smoothing(dino_noisy, 1);
    plot(ps1(:,1),ps1(:,2),'b-', 'LineWidth', 2);hold on;
    legend({'Curve noise','Curve snake-smoothing'}, 'FontSize', 15,'Interpreter','latex');
    subplot(2,2,3);
    plot(dino_noisy(:,1),dino_noisy(:,2),'r-', 'LineWidth', 2);hold on;
    ps2 = curve_smoothing(dino_noisy, 2);
    plot(ps2(:,1),ps2(:,2),'g-', 'LineWidth', 2);hold on;    
    legend({'Curve noise','Curve implicit approach'}, 'FontSize', 15, 'Interpreter','latex');
    subplot(2,2,4);
    plot(dino_noisy(:,1),dino_noisy(:,2),'r-', 'LineWidth', 2);hold on;
    ps3 = curve_smoothing(dino_noisy, 3);
    plot(ps3(:,1),ps3(:,2),'c-', 'LineWidth', 2);hold on;
    legend({'Curve noise','Curve length\&curvature min'}, 'FontSize', 15, 'Interpreter','latex');
    title('Curve smoothing', 'FontSize', 15, 'Interpreter','latex');
end

%% ex4: image unwarping
if isempty(find(skip == 4,1))
    imgpath = strcat(fullfile(pwd),strcat(data_dir,'dental/'));
    imgFiles = dir(fullfile(imgpath,'*.png'));
    imgFileNames = {imgFiles.name}';
    im = imread(strcat(imgpath,imgFileNames{1}));
    images = zeros(size(imgFileNames,1), size(im,1),size(im,2));

    axis tight manual % this ensures that getframe() returns a consistent size
    giffilename = 'unwarping.gif';
    for j = 1:1:size(imgFileNames,1)
        im = imread(strcat(imgpath,imgFileNames{j}));
        im = im2double(im);
        cc = cc_detect(im);
    %     [M,N] = size(im);
    %     cc = [M/2,N/2];
        Iw = polar_unwarping(im,[cc(1),cc(2)]);    
%         imshow(Iw,[]);
        Ic = img_merge(im,Iw);
        [imind,cm] = gray2ind(Ic,256); 
        imshow(Ic,[]);
        title('Image Warpping','FontSize',15,'Interpreter','latex');
        if drawGif == 1
            if j == 1
                imwrite(imind,cm,giffilename,'gif', 'Loopcount',inf); 
            else
                imwrite(imind,cm,giffilename,'gif','WriteMode','append'); 
            end
        end
    end
end

%% ex5: 3D images
if isempty(find(skip == 5,1))
    imgpath = strcat(fullfile(pwd),strcat(data_dir,'dental/'));
    imgFiles = dir(fullfile(imgpath,'*.png'));
    imgFileNames = {imgFiles.name}';
    im = imread(strcat(imgpath,imgFileNames{1}));
    images = zeros(size(imgFileNames,1), size(im,1),size(im,2));    
    [m,n] = size(im);
    p = size(imgFileNames,1);
    [X,Y,Z] = meshgrid(1:m,1:n,1:p);
    V = zeros(m,n,p);
    for j = 1:1:size(imgFileNames,1)
        im = imread(strcat(imgpath,imgFileNames{j}));
        if size(im,3)~=1
            im = rgb2gray(im);
        end
        threshold = graythresh(im);
        imbw = imbinarize(im,threshold);
        im(~imbw(:)) = 0;
        imshow(im);
        V(:,:,j) = im2double(im);
    end
    pp = patch(isosurface(X,Y,Z,V,0));
    isonormals(X,Y,Z,V,pp)
    pp.FaceColor = [153/255,153/255,153/255];
    pp.EdgeColor = 'none';
    daspect([1 1 1]);
    view(3);
    axis tight;
    camlight;
    lighting gouraud
    title('Volumetric Image','Interpreter','latex','FontSize',15);
end

%% ex6: PCA
if isempty(find(skip == 6,1))
    imgpath = strcat(fullfile(pwd),strcat(data_dir,'mixed_green/'));
    imgFiles = dir(fullfile(imgpath,'*.png'));
    imgFileNames = {imgFiles.name}';
    im = imread(strcat(imgpath,imgFileNames{1}));
    images = zeros(size(imgFileNames,1), size(im,1),size(im,2));    
    [M,N] = size(im);
    channels = size(imgFileNames,1);
    X = zeros(M*N,channels);
    ims = zeros(M,N,channels);
    for j = 1:1:size(imgFileNames,1)
        im = imread(strcat(imgpath,imgFileNames{j}));
        im = im2double(im);
        imshow(im);
        % just store
        ims(:,:,j) = im;
        X(:,j) = im(:);
    end
    % PCA
    Xmean = mean(X,1);
    Xbar = (X-ones(M*N,1)*Xmean);
    Cov = 1/(M*N).*(Xbar'*Xbar);
    [evec, eval] = eig(Cov);
    % MATLAB routine
    coeff = pca(X,'Algorithm','eig');
    evec = fliplr(evec);
    % compare
    err = zeros(size(evec,2),1);
    for i = 1:size(evec,2)
        eveccorr = evec(:,i).*(sign(evec(1,i))*sign(coeff(1,i)));
        err(i) = max(abs(eveccorr-coeff(:,i)));
    end
    disp(['Maximum error: ',num2str(max(err))]);
    Q = Xbar*evec;
    img = zeros(M,N,size(Q,2));
    for i = 1:size(Q,2)
        Qi = Q(:,i);
        Qi = reshape(Qi,M,N);
        imshow(Qi);
        img(:,:,i) = Qi;
%         imwrite(Qi,sprintf('./EX_1_data/pca/pca_%02d.png',i));
    end
    figure
    imshow(img(:,:,1:3));
    title('PCA result: the first 3 channels','FontSize',15,'Interpreter','latex');
end

% ex7: bacterial growth from movie frames
if isempty(find(skip == 7,1))
    v = VideoReader(strcat(data_dir,'listeria_movie.mp4'));%Read all video frames.
    ncs = [];
    dg = derivative_gaussian1d_generator(1, 5, 1);
    g = derivative_gaussian1d_generator(1, 5, 0);
    k = 1;
    axis tight manual % this ensures that getframe() returns a consistent size
    giffilename = 'bac.gif';
    while hasFrame(v)
        frame = readFrame(v);
        [num,Ic] = cell_counter(frame,dg,g);
        ncs = [ncs;num];
        [imind,cm] = rgb2ind(Ic,256); 
        if drawGif == 1
            if k == 1
                imwrite(imind,cm,giffilename,'gif', 'Loopcount',inf); 
            else
                imwrite(imind,cm,giffilename,'gif','WriteMode','append'); 
            end
        end
        k = k + 1;
    end    
    figure;
    red_color = [153/255,0,0];
    plot(ncs,'-', 'LineWidth', 3, 'Color', red_color);
    xlabel('time','Interpreter','latex','FontSize',15);
    ylabel('Number','Interpreter','latex','FontSize',15)
    title('Bacterial growth plot','Interpreter','latex','FontSize',15);
end

function [num,Imm] = cell_counter(Irgb,dg,g)
    Irgb = im2double(Irgb);
    % rgb2grey
    Igray = rgb2gray(Irgb);
    Igray = im2double(Igray);
    % gradient
%     Ih = edge(Igray,'sobel',[],'horizontal');
%     Iv = edge(Igray,'sobel',[],'vertical');
    Ih = imfilter(Igray, dg, 'replicate');
    Ih = imfilter(Ih, g', 'replicate');% horizontal
    Iv = imfilter(Igray, dg', 'replicate');% 
    Iv = imfilter(Iv, g, 'replicate');% vertical
    M = sqrt(Ih.^2+Iv.^2);% amplitude
    % smoothing
    gx = gaussian_kernel_calculator(1, 3, 5);
    Ms = imfilter(M,gx,'replicate');% along x direction
    Ms = imfilter(Ms,gx','replicate');% again along y direction
    threshold = max(Ms(:))*0.3;
    Mb = imbinarize(Ms,threshold);
    % segmentation
%     tic
    segments = img_segmentation_fast(Mb);
%     toc
    Ic = cat(3, Mb, Mb, Mb);
    % visulization
    num_segments = size(segments, 2);
    ccmap = jet(num_segments);
    Ic = im2double(Ic);
    [M,N] = size(Igray);
    num_cells = 0;
    for i = 1:num_segments
        segment = segments{i};
        indices = sub2ind(size(Igray), segment.point_set(1,:), segment.point_set(2,:));
        Ic(indices) = ccmap(i,1);
        Ic(indices+M*N) = ccmap(i,2);
        Ic(indices+M*N*2) = ccmap(i,3);
        center = round(mean(segment.point_set,2));
        center = flipud(center);
        Ic = insertText(Ic,center',num2str(size(segment.point_set,2)));
        num_cells = num_cells + size(segment.point_set,2);
    end
    Imm = img_merge(Irgb,Ic);
    Imm = imresize(Imm,[400,800]);
    imshow(Imm,[]);
    num = num_cells;
end


