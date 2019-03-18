clc; clear; close all
% Lecture 5 - Supervised Classification

addpath('..\Spectral Labs')
addpath('..\')

% read uint8 images (unsigned interger 8)
X_rgb = imread('udsnit_2005_08_18_farve.tiff'); % RGB
X_ir = imread('udsnit_2005_08_18_band453.tiff'); % Infrared

% Combined image
X = cat(3,X_rgb,X_ir); % 6-band image

[nrows, ncols, nvars] = size(X);
nobs = nrows*ncols;

%% Draw traing areas (with no overlap), for example these

% % Reading of interest polygon
% t1 = roipoly(X(:,:,1:3)); % water
% t2 = 2*roipoly(X(:,:,1:3)); % forest
% t3 = 3*roipoly(X(:,:,1:3)); % urban
% t4 = 4*roipoly(X(:,:,1:3)); % meadow
% 
% save training_areas.mat t1 t2 t3 t4
% 
% close all
%% Training image
close all; clc
load training_areas.mat
X_train = t1+t2+t3+t4;
clear t1 t2 t3 t4
nclass = length(unique(X_train))-1; % -1 since backround is a unique too

MakeFig(100, 100, 600, 600)

subplot(2,1,1)
imshow(X(:,:,1:3))
title('Original Picture of Copenhagen from 2005')

subplot(2,1,2)
imshow(X_train, [0, 4])
colormap jet
axis equal
xlim([1, ncols])
ylim([1, nrows])
title('Drawn polygons Copenhagen from 2005')

%% Computing Mahalanobis distances and likelihood
close all; clc
% Arranging image into nobs x nvars matrix
X_mat = double(reshape(X, nobs, nvars));
d = nan(nobs, nclass);
t = cell(4, 1);
for i = 1:nclass
    % From the training image, we extract the pixel intensity of each band
    % for each class
    t{i} = X_mat(X_train==i,:);
    d(:,i) = mahal(X_mat,t{i});
%     [d(:,i),S,m] = mahal1(X,t);
%     % Plotting
%     figure, imshow(reshape(d(:,i),nrows,ncols),[])
    % Dispersion of class
    Sigma_c = cov(t{i}); % We will use the determinant of this in the likelihood
    % Prior conditions P(omega_c)
    prior = 1/nclass;
    % Likelihood
    L(:, i) = - log(prior) - 1/2*log(det(Sigma_c)) - 1/2*d(:, i);
end

%% Finding the minimum Mahalanobis distance and assign class
[mind, class_1] = min(d,[],2);
figure, imshow(reshape(class_1,nrows,ncols),[0 nclass])
% Creating colormap
classlut = [0 0 0;0 0 255; 100 255 0; 230 150 0; 150 100 125; 255 255 0]/255;
colormap(gca,classlut)
%% Finding the maximum likelihood and assign class
[maxL, class_2] = max(L, [], 2);
figure, imshow(reshape(class_2,nrows,ncols),[0 nclass])
% Creating colormap
classlut = [0 0 0;0 0 255; 100 255 0; 230 150 0; 150 100 125; 255 255 0]/255;
colormap(gca,classlut)

%% Confusion matrix
% class_sample = randsample

% Known class
g1 = reshape(X_train, [], 1);
g1(g1 == 0) = nan; % Background (non-polygon) set to nan
% Predicted class
g2 = class_2;

C = confusionmat(g1, g2);
% confusionchart(C);

C = [C, sum(C, 2)];
C = [C; sum(C, 1)];
C = [[C, zeros(nclass + 1, 1)]; zeros(1, nclass + 2)];
C(end-1:end, end) = nan;
C(end, end-1) = nan;

% Accuracies
for i = 1:nclass
    C(i, end) = round(100*C(i, i)/C(i, end-1), 1);
    C(end, i) = round(100*C(i, i)/C(end-1, i), 1);
end
C
 
%% PCA

% Column centering -> unmeaning
X_unmean = X_mat - repmat(mean(X_mat), ncols*nrows, 1);
% Variance covariance (dispersion) matrix
S = cov(X_unmean);

% Eigenvalues and eigenvectors
[V, Lambda] = eig(S);
lambda = diag(Lambda);

% Sorting eigen values in descendig order, and using index to sort eigen
% vectors likewise
[lambda, sort_idx] = sort(lambda,'descend');
V = V(:, sort_idx);

% Principal components
Q = X_unmean*V;
% Rearranging Principal Components
Q_im = reshape(Q, nrows, ncols, nvars);

MakeFig(10, 10, 1200, 400)
for i = 1:nvars
    subplot(2, 3, i)
    imagesc(Q_im(:, :, i))
    colormap gray
    axis equal
    xlim([1, ncols])
    ylim([1, nrows])
    title(['PC nr. ',num2str(i),''])
end

% We check if our PCA is correct if cov(Q) = Lambda
isequal(round(diag(cov(Q)), 2), round(lambda, 2))

%% Checking variance proportions
% Proportion
prop = zeros(size(lambda));
for i = 1:length(prop)
    prop(i) = sum(lambda(1:i))/sum(lambda)*100;
end
format long g
disp([lambda,prop])
format short

% We see that including just the first 3 principal components will give us
% almost 99% of the variance in the data
%% Likelihood of class
% Data is the first 3 principal components
PCA_data = Q(:, 1:2);

for i = 1:nclass
    % From the training image, we extract the pixel intensity of each band
    % for each class
    t{i} = PCA_data(X_train==i,:);
    d(:,i) = mahal(PCA_data,t{i});
%     [d(:,i),S,m] = mahal1(X,t);
%     % Plotting
%     figure, imshow(reshape(d(:,i),nrows,ncols),[])
    % Dispersion of class
    Sigma_c = cov(t{i}); % We will use the determinant of this in the likelihood
    % Prior conditions P(omega_c)
    prior = 1/nclass;
    % Likelihood
    L(:, i) = - log(prior) - 1/2*log(det(Sigma_c)) - 1/2*d(:, i);
end

%% Finding the maximum likelihood and assign class for PC's
[maxL, class_3] = max(L, [], 2);
figure, imshow(reshape(class_3,nrows,ncols),[0 nclass])
% Creating colormap
classlut = [0 0 0;0 0 255; 100 255 0; 230 150 0; 150 100 125; 255 255 0]/255;
colormap(gca,classlut)

%% Difference in methods
class_diff = reshape(class_3,nrows,ncols) - reshape(class_2,nrows,ncols);
imshow(class_diff, [min(min(class_diff)) max(max(class_diff))])