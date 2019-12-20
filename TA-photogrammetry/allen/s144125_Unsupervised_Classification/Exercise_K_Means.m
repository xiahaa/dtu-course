% Lecture 3 - Unsupervised Classification - K-means
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

% Binary file is loaded as one long vector, so we arrange it according to 
% dimensions. 
X = reshape(X_vec, ncols, nrows, nvars);

% All combinations of of the 4 bands
combs = nchoosek(1:nvars,3);

MakeFig(100, 100, 800, 600)
for i = 1:nvars
    subplot(2, 2, i)
    imshowrgb(X, combs(i, :))
    title(['Image using band ',num2str(combs(i, 1)),', ',num2str(combs(i, 2)),' and ',num2str(combs(i, 1)),''])
end

%% Evaluation data
% We arrange the matrix so we have each pixels values in all spectra
% side-by-side
% X as a matrix
X_mat = reshape(X_vec, ncols*nrows, nvars);

tic
eva = evalclusters(X_mat,'kmeans','CalinskiHarabasz','KList', 1:6)
toc

%% K-means function
clc; close all

% Number of Clusters
N_k = 10;

% Computing k-means
[IDX, ~, ~, ~] = kmeans(X_mat, N_k);

% Arranging back so that we can visualize the results
X_kmean = reshape(IDX, ncols, nrows);

imagesc(X_kmean)
colorbar
colormap jet

%% Elbow function
% We perform the operation for cases of classes 1 - 10. For each sum of
% Euclidean distances, we save the value and plot an elbow-like curve

tic
% Sum of summed distances
sum_sq = zeros(N_k, 1);

for i = 1:N_k
    % Computing k-means
    [~, ~, sumd, ~] = kmeans(X_mat, i);
    sum_sq(i) = sum(sumd);
end
toc
% Elbow function
figure
plot(1:N_k, sum_sq)
title('Elbow curve')
xlabel('Number of classes')
ylabel('Sum of squared distances')
% We see that around 5 classes the minimization of errors start to become
% insignificant, thus we pick 5 classes for future cases
%% Using replicate
close all
% We use replicate to do the k-means clustering with 5 classes several
% times, and then pick the best result

% We pick 5 from elbow curve and evalcluster results
N_k = 5;

MakeFig(10, 10, 800, 300)
for i = 1:2
    tic
    [idx,C,sumd,D] = kmeans(X_mat,N_k,'MaxIter',10000,'Display','final','Replicates',10);
    toc

    % Arranging back so that we can visualize the results
    X_kmean = reshape(idx, ncols, nrows);

    % Plotting the result with smallest summed Euclidean distance
    subplot(1, 2, i)
    imagesc(X_kmean)
    colorbar
    colormap jet
    title(['Image pixels split into ',num2str(N_k),' classes'])
end
%% PCA
% We now repeat the process using PCA

% Subtracting mean
X_means = X_mat - repmat(mean(X_mat), ncols*nrows, 1);
% Variance covariance (dispersion) matrix
% MATLAB function
S = cov(X_means);

% Eigenvalues and eigenvectors
[V,Lambda] = eig(S);
lambda = diag(Lambda);

% Principal components
Q = X_means*V;
% Rearranging Principal Components
Q_im = reshape(Q, nrows, ncols, nvars);

MakeFig(10, 10, 800, 400)
subplot(1, 2, 1)
imagesc(Q_im(:, :, 3))
colormap gray

subplot(1, 2, 2)
imagesc(Q_im(:, :, 4))
colormap gray

%% Elbow function
N_k = 10;

tic
% Sum of summed distances
sum_sq = zeros(N_k, 1);

% We use only the 2 first principal components
for i = 1:N_k
    % Computing k-means
    [~, ~, sumd, ~] = kmeans(Q(:, 1:2), i);
    sum_sq(i) = sum(sumd);
end
toc
% Elbow function
figure
plot(1:N_k, sum_sq)

% We get that 4 clusters would once more suffice

%% Replicate with PCA
close all
% We use replicate to do the k-means clustering with 4 classes several
% times, and then pick the best result

N_k = 4;
tic
[idx_Q, ~, ~, ~] = kmeans(Q(:, 1:2),N_k,'MaxIter',10000,'Display','final','Replicates',10);
toc

% Arranging back so that we can visualize the results
PCA_kmean = reshape(idx_Q, ncols, nrows);

% Plotting the result with smallest summed Euclidean distance
imagesc(PCA_kmean)
colorbar
colormap jet

% We see that the results are quite plagued with noise when using the first
% two principle components

