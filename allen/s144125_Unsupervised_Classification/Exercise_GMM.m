% Lecture 3 - Unsupervised Classification - GMM
clc; clear; close all

% Adds a path into our data folder
addpath('..\Spectral Labs')
addpath('..\') % For the MakeFig.m function

fid = fopen('igalmss','r');
X_vec = fread(fid,'uint8');
fclose(fid);

% From igalmss.hdr we get dimensions
nrows = 512;
ncols = 512;
nvars = 4;

% We arrange the matrix so we have each pixels values in all spectra
% side-by-side
% X as a matrix
X_mat = reshape(X_vec, ncols*nrows, nvars);
% For displaying image
X = reshape(X_vec, ncols, nrows, nvars);

% All combinations of of the 4 bands
combs = nchoosek(1:nvars,3);

MakeFig(100, 100, 400, 400)
for i = 1:nvars
    subplot(2, 2, i)
    imshowrgb(X, combs(i, :))
    title(['Image using band ',num2str(combs(i, 1)),', ',num2str(combs(i, 2)),' and ',num2str(combs(i, 1)),''])
end

%% Evaluating our clusters
tic
eva = evalclusters(X_mat,'gmdistribution','CalinskiHarabasz','KList', 1:6)
toc

% Evaluating the pixel data we have, evalclusters tells us that the optimal
% number of classes with GMM is 2. We therefore move forward with 2
%% Classyfying
% Number of clusters
N_k = eva.OptimalK;
N_k = 5;

GMModel = fitgmdist(X_mat, N_k);
%% Inspecting GMModel
clc; close all
% Colormap
c_map = jet;
% Length of jet colormap
N_colormap = size(c_map, 1);
% Index for color of each classes
Idx_class = round(linspace(1, N_colormap, N_k));
% With "posterior" we get the how much each pixel belongs to each class
P = posterior(GMModel,X_mat);
% Number of pixels
N_pix = size(P, 1);

C = zeros(N_pix, 1);

for i = 1:N_pix
    C(i) = find(P(i, :) == max(P(i, :)));
end
C = reshape(C, ncols, nrows);

MakeFig(100, 100, 800, 400)
subplot(1,2,1)
imagesc(C)
colormap jet
title(['Pixels assigned directly to classes (',num2str(N_k),' classes)'])

% Assigning a percentage of each color spectrum, based of proportions (pi)
C_rgb = zeros(N_pix, 3);

for i = 1:N_pix
    for j = 1:N_k
        C_rgb(i, :) = C_rgb(i, :) + P(i, j)*c_map(Idx_class(j), :);
    end
end

C_rgb = reshape(C_rgb, ncols, nrows, 3);

subplot(1,2,2)
imagesc(C_rgb)
colormap jet
title(['Pixels assigned proportionately to classes (',num2str(N_k),' classes)'])
